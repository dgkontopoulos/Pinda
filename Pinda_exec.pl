#!/usr/bin/perl

####################
#$VERSION = '0.02';#
####################

=head1 NAME

PINDA - Pipeline for INtraspecies Duplication Analysis

=head1 AUTHOR

Pinda has been developed by Dimitrios - Georgios Kontopoulos as his
final year project, under the supervision of Prof. Nicholas M. Glykos
at the Department of Molecular Biology and Genetics of Democritus
University of Thrace, Greece.

=head1 LICENSE

This program is free software: you can redistribute it and/or modify 
it under the terms of the GNU Affero General Public License as 
published by the Free Software Foundation, either version 3 of the 
License, or (at your option) any later version.

For more information, see http://www.gnu.org/licenses/.

=cut

use Bio::AlignIO;
use Bio::Perl;
use Bio::Root::IO;
use Bio::SearchIO;
use Bio::Search::Iteration::GenericIteration;
use Bio::Search::Result::BlastResult;
use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlastPlus;
use Bio::Tools::Run::StandAloneBlastPlus::BlastMethods;
use Bio::TreeIO;
use Bio::Tree::TreeFunctionsI;
use Bio::Tree::TreeI;
use File::stat;
use FreezeThaw qw(freeze thaw);
use List::MoreUtils qw(uniq);
use LWP::Simple qw(get);
use Math::Cephes qw(:all);
use MIME::Lite;
use Statistics::Basic qw(:all nofill);
use Sys::CPU;
use Time::localtime;

use feature qw(say);

use bignum;
use strict;
use warnings;

open STDERR, '>', '/tmp/err' or die $!;

########################################
#On to the serious computation phase...#
########################################
my $seq2counter      = 0;
my $hit_old          = 0;
my $iteration_number = 0;
my $resnum           = 0;
my $resnum2          = 0;
my $sequences_number = 0;
my $tdcounter        = 0;

my (
    $acnumber,   $alignment,  $dbfetch,      $des,
    $hit,        $i,          $line,         $match_line,
    $number,     $one_gos,    $org1,         $out,
    $out_fh,     $sequence,   $sequences_fh, $starting_point,
    $tdbg,       @accession,  @candidate,    @cand_sans,
    @sequences2, @realign,    @reslines,     @seq,
    @seq2,       %common,     %hsp,          %hsp_pos,
    %hsp_seq,    %textcommon, %textncommon,  %texts,
    %texts2
);
my $start_timer = time;

my $email         = $ARGV[0];
my $organism      = $ARGV[1];
my $prid          = $ARGV[2];
my $db            = $ARGV[3];
my $lcr_filtering = $ARGV[4];
my $one           = $ARGV[5];
my $masking       = $ARGV[6];

my $slurm_queue = `squeue`;
my $slurm_id;
my $squeue_line;
{
    local $/ = "\n";
    open( my $squeue_fh, '<', \$slurm_queue );
    while ( $squeue_line = <$squeue_fh> )
    {
        if ( $squeue_line =~ /1 orion/ )
        {
            if ( $squeue_line =~ /\s+(\d+)\s+/ )
            {
                $slurm_id = $1;
                last;
            }
        }
    }
    close $squeue_fh;
}

my $tmp = '../tmps/blast/' . $prid . '.tmp';
my $query_line;
if ( $one eq 'QUERY' )
{
    open my $tmp_fh, '<', $tmp;
    local $/ = undef;
    while ( my $line = <$tmp_fh> )
    {
        if ( $line =~ />.*\n/ )
        {
            $query_line = $';
            $query_line = uc $query_line;
            last;
        }
        elsif ( $line !~ />/ )
        {
            $query_line = $line;
            $query_line = uc $query_line;
            last;
        }
    }
    close $tmp_fh;
}
my $query_length = length $query_line;
my $SWISSPROT    = '/usr/local/databases/Swissprot/uniprot_sprot.fasta';
my $UNIPROT      = '/usr/local/databases/UniProt/UniProt.fasta';
my $NT           = '/usr/local/databases/nt/nt.fasta';
my $database     = $db;

if ( $db eq 'Swiss-Prot' )
{
    $db = $SWISSPROT;
}
elsif ( $db eq 'UniProt' )
{
    $db = $UNIPROT;
}
elsif ( $db eq 'nt' )
{
    $db = $NT;
}

my $input = '/var/www/Pinda/tmps/blast/' . $prid . '.tmp';
open my $input_fh, '<', $input or die $!;
my $line_input;
{
    local $/ = "\n";
    while ( $line = <$input_fh> )
    {
        $line_input .= $line . '<br>';
    }
    if ( $line_input =~ /<br>$/ )
    {
        $line_input = $`;
    }
}
close $input_fh;
chomp $line_input;
$line_input =~ s/\n//g;

my $email_data = <<"EMAIL_END";
<center><br>
<a href='http://orion.mbg.duth.gr/Pinda/cgi/Pinda.cgi'>
<img src='http://orion.mbg.duth.gr/Pinda/pindalogo.png'></a>
</center><br><br><br><hr />
<b>Input Sequence:<br><font face='Courier New'>$line_input
</font><br></b></b></b><hr />
<b>Organism:</b> $organism<br><hr />
<b>Database:</b> $database<hr /><br><center>
EMAIL_END

if ( $db =~ /nt[.]fasta/ )
{
    $email_data .= <<"EMAIL_END";
	<font size='2'>
	<img src='http://orion.mbg.duth.gr/Pinda/caution.png' width=39
	height=35><b><br>
	The nt database is partially non-redundant.<br>
	You are advised to pay extra attention when interpreting the results.
	</b></font><br>
EMAIL_END
}

##############################
#Get the number of cpu cores.#
##############################
my $cpu_n = Sys::CPU::cpu_count();
####################
#E-value threshold.#
####################
if ( $db !~ /nt[.]fasta/ )
{
    my $e_th = '0.00000000001';

    $out = '../outs/psiblast/' . $prid . '.tmp';

    if ( $lcr_filtering == 1 )
    {
        system(
"psiblast -query $tmp -num_alignments 7000 -num_iterations 50 -evalue $e_th -db $db -num_threads $cpu_n -out $out -seg yes"
        );
    }
    else
    {
        system(
"psiblast -query $tmp -num_alignments 7000 -num_iterations 50 -evalue $e_th -db $db -num_threads $cpu_n -out $out -seg no"
        );
    }
    #####################################################################
    #Shorten the PSI-BLAST output file. We only need the last iteration.#
    #####################################################################
    open $out_fh, '<', $out or die $!;
    while ( $line = <$out_fh> )
    {
        if ( $line =~ /Results from round (\d+)/ )
        {
            $resnum = $1;
        }
    }
    close $out_fh or die $!;
    open $out_fh, '<', $out or die $!;
    while ( $line = <$out_fh> )
    {
        if ( $line =~ /Results from round $resnum/ )
        {
            while ( $line = <$out_fh> )
            {
                $reslines[$resnum2] = $line;
                $resnum2++;
            }
        }
    }
    close $out_fh or die $!;
    open $out_fh, '>', $out or die $!;
    foreach my $line (@reslines)
    {
        print {$out_fh} $line;
    }
    close $out_fh or die $!;

    my $blast = Bio::SearchIO->new(
        -format => 'blast',
        -file   => "$out"
    );

    my $query_found = 0;
    #################################
    #Start parsing PSI-BLAST output.#
    #################################
    while ( my $result = $blast->next_result )
    {
        while ( my $it = $result->next_iteration )
        {
            $i = 1000.0;
            my $number = 0;
            while (( $hit = $it->next_hit_new )
                || ( $hit = $it->next_hit_old ) )
            {
                if ( ( $it->number ) > $iteration_number )
                {
                    $iteration_number = $it->number;
                }
                if ( $hit->description =~ /(OS\=\w+)\s+(\w+)/ )
                {
                    if ( $hit->accession =~ /tr[|](\w+)[|]/ )
                    {
                        my $ac = $1;
                        local $/ = undef;
                        open my $out_fh, '<', $out;
                        while ( my $readline = <$out_fh> )
                        {
                            if ( $readline =~ />tr[|]$ac[|]/ )
                            {
                                $readline = $';
                                if ( $readline =~ /OS\=(.+\s+.+)\s+/ )
                                {
                                    $org1 = $1;
                                    if ( $org1 =~ /\n/ )
                                    {
                                        $org1 = $` . $';
                                    }
                                }
                                else
                                {
                                    undef $org1;
                                }
                            }
                        }
                        close $out_fh;
                    }
                    else
                    {
                        my $ac = $hit->accession;
                        local $/ = undef;
                        open my $out_fh, '<', $out;
                        while ( my $readline = <$out_fh> )
                        {
                            if ( $readline =~ />sp[|]$ac[|]/ )
                            {
                                $readline = $';
                                if ( $readline =~ /OS\=(.+\s+.+)\s+/ )
                                {
                                    $org1 = $1;
                                    if ( $org1 =~ /\n/ )
                                    {
                                        $org1 = $` . $';
                                    }
                                }
                                else
                                {
                                    undef $org1;
                                }
                            }
                        }
                    }
                    if ( defined $org1 && $org1 =~ /$organism/i )
                    {
                        if ( $hit->accession =~ /tr[|](\w+)[|]/ )
                        {
                            $accession[$number] = $1;
                            $acnumber = $accession[$number];
                        }
                        else
                        {
                            $accession[$number] = $hit->accession;
                            $acnumber = $accession[$number];
                        }
                        ###################################
                        #Get the full sequence of the hit.#
                        ###################################
                        if ( $hit_old ne $acnumber )
                        {
                            $dbfetch = get(
"http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=uniprotkb&id=$acnumber&format=fasta&style=raw"
                            );
                            if ( !( defined $dbfetch ) )
                            {
                                for ( 0 .. 3 )
                                {
                                    $dbfetch = get(
"http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=uniprotkb&id=$acnumber&format=fasta&style=raw"
                                    );
                                    if ( defined $dbfetch )
                                    {
                                        last;
                                    }
                                }
                            }
                            if ( defined $dbfetch && $dbfetch =~ /\n/ )
                            {
                                $match_line = $';
                                my $seq_length = length $match_line;
                                if (   $seq_length > $query_length / 5
                                    && $seq_length <= 3 * $query_length )
                                {
                                    $seq[$number] =
">$accession[$number] $org1\n$match_line\n\n";
                                    $number++;
                                }
                                elsif ( $seq_length > $query_length / 5 )
                                {
                                    my $hsp_old_length;
                                    while ( my $hsp = $hit->next_hsp )
                                    {
                                        if ( !( defined $hsp_old_length )
                                            || $hsp->length('hit') >
                                            $hsp_old_length )
                                        {
                                            $hsp{ $accession[$number] } =
                                              ">$accession[$number] $org1\n"
                                              . $hsp->hit_string . "\n\n";
                                            $hsp_seq{ $accession[$number] } =
                                              $hsp->hit_string;
                                            $hsp_pos{ $accession[$number] } =
                                                '('
                                              . $hsp->start('hit') . ' -> '
                                              . $hsp->end('hit') . ')';
                                            $hsp_old_length =
                                              $hsp->length('hit');
                                        }
                                    }
                                    $seq[$number] = $hsp{ $accession[$number] };
                                    undef $hsp_old_length;
                                    delete $hsp{ $accession[$number] };
                                    $number++;
                                }
                            }
                            if ( defined $query_line )
                            {
                                $query_line =~ s/\n//g;
                                $query_line =~ s/\s//g;
                                $match_line =~ s/\n//g;
                                $match_line =~ s/\s//g;
                                if ( ( uc $query_line ) eq ( uc $match_line ) )
                                {
                                    $query_found = '1';
                                    $one         = $acnumber;
                                }
                            }
                        }
                        $hit_old = $acnumber;
                    }
                }
            }
        }
    }
    if ( defined $query_line && $query_found == '0' )
    {
        my $number = @seq;
        $seq[$number] = ">QUERY $organism\n" . $query_line . "\n\n";
    }
}
else
{
    $cpu_n--;
    $out = '../outs/blast/' . $prid . '.tmp';
    if ( $lcr_filtering == 1 )
    {
        system(
"blastn -query $tmp -db $db -evalue 0.00000000001 -num_threads $cpu_n -out $out -dust yes"
          ) == 0
          or die $?;
    }
    else
    {
        system(
"blastn -query $tmp -db $db -evalue 0.00000000001 -num_threads $cpu_n -out $out -dust no"
          ) == 0
          or die $?;
    }

    my $blast = Bio::SearchIO->new(
        -format => 'blast',
        -file   => "$out"
    );

    my $hit_old = 0;
    my ( $accession, $input_hit, $org, @organism );
    my $list        = 0;
    my $number      = 0;
    my $query_found = 0;
    #############################
    #Start parsing BLAST output.#
    #############################
    while ( my $result = $blast->next_result )
    {
        while ( my $hit = $result->next_hit )
        {
            if ( $hit->name =~ /ref[|](.+)[|]/ )
            {
                if ( $1 =~ /[.]\d*/ )
                {
                    $accession = $`;
                }
                $dbfetch = get(
"http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=refseqn&id=$accession&style=raw"
                );
                if ( !( defined $dbfetch ) )
                {
                    for ( 0 .. 3 )
                    {
                        $dbfetch = get(
"http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=refseqn&id=$accession&style=raw"
                        );
                        if ( defined $dbfetch )
                        {
                            last;
                        }
                    }
                }

                ######################################
                #Populate the organism dropdown list.#
                ######################################
                open( my $db_fh, '<', \$dbfetch );
                local $/ = "\n";
                while ( my $db_line = <$db_fh> )
                {
                    if ( $db_line =~ /ORGANISM\s+(.+\s+.+)/ )
                    {
                        $org = $1;
                        $organism[$list] = $org;
                        $list++;
                        last;
                    }
                }
                close $db_fh;
            }
            elsif ($hit->name =~ /gb[|](.+)[|]/
                || $hit->name =~ /dbj[|](.+)[|]/
                || $hit->name =~ /emb[|](.+)[|]/
                || $hit->name =~ /tpg[|](.+)[|]/ )
            {
                if ( $1 =~ /[.]\d*/ )
                {
                    $accession = $`;
                }
                $dbfetch = get(
"http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=embl&id=$accession&style=raw"
                );
                if ( !( defined $dbfetch ) )
                {
                    for ( 0 .. 3 )
                    {
                        $dbfetch = get(
"http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=embl&id=$accession&style=raw"
                        );
                        if ( defined $dbfetch )
                        {
                            last;
                        }
                    }
                }

                ######################################
                #Populate the organism dropdown list.#
                ######################################
                open( my $db_fh, '<', \$dbfetch );
                local $/ = "\n";
                while ( my $db_line = <$db_fh> )
                {
                    if ( $db_line =~ /OS\s+(.+\s+.+)/ )
                    {
                        $org = $1;
                        $organism[$list] = $org;
                        $list++;
                        last;
                    }
                }
            }
            ####################
            #Get the sequences.#
            ####################
            if ( defined $accession && $org =~ /$organism/i )
            {
                if ( !( defined $hit_old ) || $hit_old ne $accession )
                {
                    if ( $accession =~ /[_]/ )
                    {
                        $dbfetch = get(
"http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=refseqn&id=$accession&format=fasta&style=raw"
                        );
                        if ( !( defined $dbfetch ) )
                        {
                            for ( 0 .. 3 )
                            {
                                $dbfetch = get(
"http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=refseqn&id=$accession&format=fasta&style=raw"
                                );
                                if ( defined $dbfetch )
                                {
                                    last;
                                }
                            }
                        }
                    }
                    else
                    {
                        $dbfetch = get(
"http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=embl&id=$accession&format=fasta&style=raw"
                        );
                        if ( !( defined $dbfetch ) )
                        {
                            for ( 0 .. 3 )
                            {
                                $dbfetch = get(
"http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=embl&id=$accession&format=fasta&style=raw"
                                );
                                if ( defined $dbfetch )
                                {
                                    last;
                                }
                            }
                        }
                    }
                    if ( $dbfetch =~ /\n/ )
                    {
                        $match_line = $';
                    }
                    if ( defined $query_line )
                    {
                        $query_line =~ s/\n//g;
                        $query_line =~ s/\s//g;
                        $match_line =~ s/\n//g;
                        $match_line =~ s/\s//g;
                        if ( ( uc $query_line ) eq ( uc $match_line ) )
                        {
                            $query_found = '1';
                            $one         = $accession;
                            $accession   = "III" . $accession . "III";
                        }
                    }

                    $match_line = uc $match_line;
                    $match_line =~ tr/\n//d;
                    my $seq_length = length $match_line;
                    if (   $seq_length > $query_length / 5
                        && $seq_length <= 3 * $query_length )
                    {
                        $seq[$number] = ">$accession $org\n$match_line\n\n";
                        $number++;
                    }
                    elsif ( $seq_length > $query_length / 5 )
                    {
                        my $hsp_old_length;
                        while ( my $hsp = $hit->next_hsp )
                        {
                            if ( !( defined $hsp_old_length )
                                || $hsp->length('hit') > $hsp_old_length )
                            {
                                $hsp{$accession} = ">$accession $org\n"
                                  . $hsp->hit_string . "\n\n";
                                $hsp_seq{$accession} = $hsp->hit_string;
                                $hsp_pos{$accession} = '('
                                  . $hsp->start('hit') . ' -> '
                                  . $hsp->end('hit') . ')';
                                $hsp_old_length = $hsp->length('hit');
                            }
                        }
                        $seq[$number] = $hsp{$accession};
                        $number++;
                        undef $hsp_old_length;
                        delete $hsp{$accession};
                    }
                }
            }
            $hit_old = $accession;
        }
    }
    if ( defined $query_line && $query_found == '0' )
    {
        my $number = @seq;
        $seq[$number] = ">IIIQUERYIII $organism\n" . $query_line . "\n";
    }

}

###########################
#Append sequences to file.#
###########################
@seq2 = uniq(@seq);
my $sequences = '/var/www/Pinda/seq/final_seq/' . $prid . '.tmp';
open $sequences_fh, '>', $sequences or die $!;

foreach my $sequence (@seq2)
{
    if ( $sequence =~ /$organism/i )
    {
        if ( $sequence =~ /(\w+)[.]\d*/ )
        {
            $sequence = $` . $1 . $';
        }
        print {$sequences_fh} $sequence;

    }
}
close $sequences_fh or die $!;
if ( $db !~ /nt[.]fasta/ )
{
    open $sequences_fh, '<', $sequences or die $!;
    my $line2;
    {
        local $/ = "\n";
        while ( $line2 = <$sequences_fh> )
        {
            ##################################
            #Add stars to the input sequence.#
            ##################################
            $line2 =~ s/$one/***$one***/;
            $sequences2[$seq2counter] = $line2;
            $seq2counter++;
        }
    }
    close $sequences_fh or die $!;

    open $sequences_fh, '>', $sequences or die $!;
    foreach my $seq2counter (@sequences2)
    {
        print {$sequences_fh} $seq2counter;
    }
    close $sequences_fh or die $!;
}

my $fnaln     = '/var/www/Pinda/results/final_alns/multalign/' . $prid . '.aln';
my $ph        = '/var/www/Pinda/results/trees/phs/' . $prid . '.ph';
my $phb       = '/var/www/Pinda/results/trees/phbs/' . $prid . '.phb';
my $drawntree = '/var/www/Pinda/results/trees/drawn/' . $prid . '.png';
my $zip       = '/var/www/Pinda/results/trees/zips/' . $prid . '.zip';

open $sequences_fh, '<', $sequences or die $!;
my $line3;
{
    local $/ = "\n\n";
    while ( $line3 = <$sequences_fh> )
    {
        if ( $line3 !~ /^\s*$/ )
        {
            ################################
            #Count the resulting sequences.#
            ################################
            $sequences_number++;
        }
    }
}
close $sequences_fh or die $!;

if ( $one ne 'QUERY' )
{
    $email_data .= <<"EMAIL_END";
	<br>This sequence has been identified as "$one" from "$organism".
	<br>
EMAIL_END
}

my $email2 = $email;
$email2 =~ s/([\@])/\\$1/;
if ( $sequences_number > 3 )
{
    #########################################################
    #Alignment, NJ tree plotting, bootstrapping and parsing.#
    #########################################################
    align( $sequences, $fnaln, $db, 1 );
    if ( $db !~ /nt[.]fasta/ && $masking == 1 && $sequences_number <= 150 )
    {
        my $conf_val =
          '/var/www/Pinda/results/final_alns/multalign/conf/' . $prid . '.tmp';
        alignment_masking( $fnaln, $conf_val ) == 0 or die $?;
    }
    elsif ( $db !~ /nt[.]fasta/ && $masking == 1 && $sequences_number > 150 )
    {
        $email_data .= <<"EMAIL_END";
		<br><center>The resulting hits were over 150.
		<br>Due to computational limits, masking was <b>NOT</b> performed.</center>
EMAIL_END
    }
    my $fnaln2 = $fnaln . ".fasta";
    if ( !( -e $fnaln2 ) )
    {
        system("cp $fnaln $fnaln2");
    }
    system("clustalw -INFILE=$fnaln2 -OUTFILE=$ph -tree") == 0 or die $?;
    system(
"clustalw -INFILE=$fnaln2 -OUTFILE=$phb -bootstrap=1000 -bootlabels=node"
      ) == 0
      or die $?;
    rename "../results/final_alns/multalign/$prid.aln.ph",
      "../results/trees/phs/$prid.ph"
      or die $!;
    rename "../results/final_alns/multalign/$prid.aln.phb",
      "../results/trees/phbs/$prid.phb"
      or die $!;
    tree_manipulation1($phb);
    ##########################################################
    #Parsing the phb tree for distances and bootstrap values.#
    ##########################################################
    system("../Pinda.R -parser $phb > /var/www/Pinda/parsing/$prid.tmp") == 0
      or die $!;
    system("../Pinda.R -lengths_1 $phb > /var/www/Pinda/parsing/$prid\_1.tmp")
      == 0
      or die $?;
    system("../Pinda.R -lengths_2 $phb > /var/www/Pinda/parsing/$prid\_2.tmp")
      == 0
      or die $?;
    tree_manipulation2($phb);
    system("zip -j $zip $ph $phb") == 0 or die $?;
    ################
    #Draw the tree.#
    ################
    system("../Pinda.R $drawntree $phb") == 0
      or die $?;
    system("rm $ph $phb") == 0 or die $?;

    parser( "../parsing/$prid.tmp", "../parsing/$prid\_1.tmp",
        "../parsing/$prid\_2.tmp" );

    ###############################################
    #Compute level of confidence for duplications.#
    ###############################################
    compute_probability(@candidate);
    my $realign_num = @realign;
    if ( $realign_num >= 1 )
    {
        open $sequences_fh, '>', $sequences or die $!;
        my $one_dbfetch;
        if ( $db =~ /nt[.]fasta/ )
        {
            if ( $one =~ /^\D{2}\_/ )
            {
                $one_dbfetch = get(
"http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=refseqn&id=$one&format=fasta&style=raw"
                );
                if ( !( defined $one_dbfetch ) )
                {
                    for ( 0 .. 3 )
                    {
                        $one_dbfetch = get(
"http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=refseqn&id=$one&format=fasta&style=raw"
                        );
                        if ( defined $one_dbfetch )
                        {
                            last;
                        }
                    }
                }
                $one_dbfetch = uc $one_dbfetch;
                if ( $one_dbfetch =~ /\n/ )
                {
                    $one_dbfetch = ">III" . $one . "III\n" . $' . "\n\n";
                }
            }
            elsif ( $one eq 'QUERY' )
            {
                $one_dbfetch = ">IIIQUERYIII\n$query_line\n\n";
            }
            else
            {
                $one_dbfetch = get(
"http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=embl&id=$one&format=fasta&style=raw"
                );
                if ( !( defined $one_dbfetch ) )
                {
                    for ( 0 .. 3 )
                    {
                        $one_dbfetch = get(
"http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=embl&id=$one&format=fasta&style=raw"
                        );
                        if ( defined $dbfetch )
                        {
                            last;
                        }
                    }
                }
                if ( $one_dbfetch =~ /\n/ )
                {
                    $one_dbfetch = ">III" . $one . "III\n" . $' . "\n\n";
                }
            }
        }
        else
        {
            if ( $one eq 'QUERY' )
            {
                $one_dbfetch = ">QUERY\n$query_line\n\n";
            }
            else
            {
                $one_dbfetch = get("http://www.uniprot.org/uniprot/$one.fasta");
            }
        }
        if ( $db !~ /nt[.]fasta/ )
        {
            $one_dbfetch =~ s/$one/***$one***/;
        }
        say {$sequences_fh} $one_dbfetch;
        my $dbfetch;
        foreach my $reseq (@realign)
        {
            if ( defined $hsp_seq{$reseq} && defined $hsp_pos{$reseq} )
            {
                $dbfetch =
                  ">$reseq $hsp_pos{$reseq}\n" . $hsp_seq{$reseq} . "\n\n";
            }
            else
            {
                if ( $db =~ /nt[.]fasta/ )
                {
                    if ( $reseq =~ /^\D{2}\_/ )
                    {
                        $dbfetch = get(
"http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=refseqn&id=$reseq&format=fasta&style=raw"
                        );
                        if ( !( defined $dbfetch ) )
                        {
                            for ( 0 .. 3 )
                            {
                                $dbfetch = get(
"http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=refseqn&id=$reseq&format=fasta&style=raw"
                                );
                                if ( defined $dbfetch )
                                {
                                    last;
                                }
                            }
                        }
                        $dbfetch = uc $dbfetch;
                        if ( $dbfetch =~ /\n/ )
                        {
                            $dbfetch = ">$reseq\n" . $' . "\n\n";
                        }
                    }
                    else
                    {
                        $dbfetch = get(
"http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=embl&id=$reseq&format=fasta&style=raw"
                        );
                        if ( !( defined $dbfetch ) )
                        {
                            for ( 0 .. 3 )
                            {
                                $dbfetch = get(
"http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=embl&id=$reseq&format=fasta&style=raw"
                                );
                                if ( defined $dbfetch )
                                {
                                    last;
                                }
                            }
                        }
                        if ( $dbfetch =~ /\n/ )
                        {
                            $dbfetch = ">$reseq\n" . $' . "\n\n";
                        }
                    }
                }
                else
                {
                    $dbfetch =
                      get("http://www.uniprot.org/uniprot/$reseq.fasta");
                    if ( !( defined $dbfetch ) )
                    {
                        for ( 0 .. 3 )
                        {
                            $dbfetch = get(
                                "http://www.uniprot.org/uniprot/$reseq.fasta");
                            if ( defined $dbfetch )
                            {
                                last;
                            }
                        }
                    }
                }
            }
            say {$sequences_fh} $dbfetch;
        }
        close $sequences_fh or die $!;
        $fnaln .= "2";
        align( $sequences, $fnaln, $db );
    }
    if ( $db !~ /nt[.]fasta/ )
    {
        ##############################
        #Find common gene ontologies.#
        ##############################
        GOs_in_common();
    }
    alarm 0;
    #######################
    #Final results output.#
    #######################
    $email_data .= <<"EMAIL_END";
		<center><br><font size='3' face='Georgia' color='330033'>
        <a href=http://orion.mbg.duth.gr/Pinda/results/final_alns/multalign/$prid.aln>Alignment of all hits</a>
EMAIL_END
    my $clu_file = "/var/www/Pinda/results/final_alns/multalign/$prid.aln.clu";
    if ( -e $clu_file )
    {
        $email_data .= <<"EMAIL_END";
		 | <a href=http://orion.mbg.duth.gr/Pinda/results/final_alns/multalign/$prid.aln.clu>Masked Alignment</a>
EMAIL_END
    }
    if ( $realign_num >= 1 )
    {
        $email_data .= <<"EMAIL_END";
         | <a href=http://orion.mbg.duth.gr/Pinda/results/final_alns/multalign/$prid.aln2>Alignment of top hits</a>
EMAIL_END
    }

    $email_data .= <<"EMAIL_END";
        </font>
        <br><br>
        <table border='1'>
        <tr bgcolor=FFFF66><th><center>Possible duplications of
EMAIL_END

    if ( $one ne 'QUERY' )
    {
        if ( $db !~ /nt[.]fasta/ )
        {
            $email_data .= <<"EMAIL_END";
            <a href=http://www.uniprot.org/uniprot/$one>
            $one</a>.</center></th>
            <th><center>Z Value</center></th>
            <th><center>Level Of Confidence</center></th>
            <th>GOs in common (out of $one_gos)</th>
            <th>GOs not found in $one</th></tr>
EMAIL_END
        }
        else
        {
            if ( $one =~ /^\D{2}\_/ )
            {
                $email_data .= <<"EMAIL_END";
                <a href=http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=refseqn&id=$one&style=raw>
                $one</a>.</center></th>
                <th><center>Z Value</center></th>
				<th><center>Level Of Confidence</center></th></tr>
EMAIL_END
            }
            else
            {
                $email_data .= <<"EMAIL_END";
            <a href=http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=embl&id=$one&style=raw>
            $one</a>.</center></th>
            <th><center>Z Value</center></th>
            <th><center>Level Of Confidence</center></th></tr>
EMAIL_END
            }
        }
    }
    else
    {
        if ( $db !~ /nt[.]fasta/ )
        {
            $email_data .= <<"EMAIL_END";
            $one.</center></th>
            <th><center>Z Value</center></th>
            <th><center>Level Of Confidence</center></th>
            <th>GO terms</th></tr>
EMAIL_END
        }
        else
        {
            if ( $one =~ /^\D{2}\_/ )
            {
                $email_data .= <<"EMAIL_END";
                $one.</center></th>
                <th><center>Z Value</center></th>
				<th><center>Level Of Confidence</center></th></tr>
EMAIL_END
            }
            else
            {
                $email_data .= <<"EMAIL_END";
                $one.</center></th>
                <th><center>Z Value</center></th>
				<th><center>Level Of Confidence</center></th></tr>
EMAIL_END
            }
        }
    }
    ##################################
    #Table of Z Values/Probabilities.#
    ##################################
    my $top_can;
    #################
    #Input sequence.#
    #################
    if ( $candidate[0] =~ /(\d?\d?\d?.?\d+e?-?\d*) \w+/ )
    {
        $top_can = $1;
    }

    foreach my $can (@candidate)
    {
        ################################################
        #List every sequence, except for the input one.#
        ################################################
        if (   $can !~ /$starting_point/
            && $can =~
            /(-?\d?\d?\d?.?\d+e?-?\d*) (\d?\d?\d?.?\d+e?-?\d*) (\w+)/ )
        {
            my $z_temp    = $1;
            my $conf_temp = $2;
            my $uni_temp  = $3;
            ###################################
            #Color the html table alternately.#
            ###################################
            if ( $tdcounter == 1 )
            {
                $tdbg      = 'F8FBFE';
                $tdcounter = 0;
            }
            else
            {
                $tdbg = 'EAF1FB';
                $tdcounter++;
            }
            if ( $db !~ /nt[.]fasta/ )
            {
                if ( defined $hsp_pos{$3} )
                {
                    $email_data .= <<"EMAIL_END";
<tr bgcolor=$tdbg><td><center><a href=http://www.uniprot.org/uniprot/$3>$3</a> $hsp_pos{$3}
</center></td>
EMAIL_END
                }
                else
                {
                    $email_data .= <<"EMAIL_END";
<tr bgcolor=$tdbg><td><center><a href=http://www.uniprot.org/uniprot/$3>$3</a>
</center></td>
EMAIL_END
                }
            }
            else
            {
                if ( $uni_temp =~ /^\D{2}\_/ )
                {
                    if ( defined $hsp_pos{$uni_temp} )
                    {
                        $email_data .= <<"EMAIL_END";
<tr bgcolor=$tdbg><td><center>
<a href=http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=refseqn&id=$uni_temp&style=raw>$uni_temp</a> $hsp_pos{$uni_temp}
</center></td>
EMAIL_END
                    }
                    else
                    {
                        $email_data .= <<"EMAIL_END";
<tr bgcolor=$tdbg><td><center>
<a href=http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=refseqn&id=$uni_temp&style=raw>$uni_temp</a>
</center></td>
EMAIL_END
                    }
                }
                else
                {
                    if ( defined $hsp_pos{$uni_temp} )
                    {
                        $email_data .= <<"EMAIL_END";
<tr bgcolor=$tdbg><td><center>
<a href=http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=embl&id=$uni_temp&style=raw>$uni_temp</a> $hsp_pos{$uni_temp}
</center></td>
EMAIL_END
                    }
                    else
                    {
                        $email_data .= <<"EMAIL_END";
<tr bgcolor=$tdbg><td><center>
<a href=http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=embl&id=$uni_temp&style=raw>$uni_temp</a>
</center></td>
EMAIL_END
                    }
                }
            }
            $email_data .=
              sprintf
'<td align=left><center>%5.2f</center></td><td align=left><center>%5.1f%%</center></td>',
              $z_temp, $conf_temp;
            if ( $one ne 'QUERY' )
            {
                if ( $db !~ /nt[.]fasta/ )
                {
                    if ( $common{$uni_temp} =~ /(\w+)\s(\w+)/ )
                    {
                        $email_data .= <<"EMAIL_END";
						<td title="$textcommon{$uni_temp}"><center>$1</center></td>
						<td title="$textncommon{$uni_temp}"><center>$2</center></td></tr>
EMAIL_END
                    }
                }
            }
            else
            {
                if ( $db !~ /nt[.]fasta/ )
                {
                    if ( $texts{$uni_temp} =~ /(\w+)/ )
                    {
                        $email_data .= <<"EMAIL_END";
						<td title="$texts2{$uni_temp}"><center>$1</center></td></tr>
EMAIL_END
                    }
                }
            }
        }
        elsif ($can !~ /$starting_point/
            && $can =~ /(-?\d?\d?\d?.?\d+e?-?\d*) (\w+)/ )
        {
            my $z_temp   = $1;
            my $uni_temp = $2;
            ###################################
            #Color the html table alternately.#
            ###################################
            if ( $tdcounter == 1 )
            {
                $tdbg      = 'F8FBFE';
                $tdcounter = 0;
            }
            else
            {
                $tdbg = 'EAF1FB';
                $tdcounter++;
            }
            if ( $db !~ /nt[.]fasta/ )
            {
                if ( defined $hsp_pos{$2} )
                {
                    $email_data .= <<"EMAIL_END";
<tr bgcolor=$tdbg><td><center>
<a href=http://www.uniprot.org/uniprot/$2>$2</a> $hsp_pos{$2}
</center></td>
EMAIL_END
                }
                else
                {
                    $email_data .= <<"EMAIL_END";
<tr bgcolor=$tdbg><td><center>
<a href=http://www.uniprot.org/uniprot/$2>$2</a>
</center></td>
EMAIL_END
                }
            }
            else
            {
                if ( $uni_temp =~ /^\D{2}\_/ )
                {
                    if ( defined $hsp_pos{$uni_temp} )
                    {
                        $email_data .= <<"EMAIL_END";
<tr bgcolor=$tdbg><td><center>
<a href=http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=refseqn&id=$uni_temp&style=raw>$uni_temp</a> $hsp_pos{$uni_temp}
</center></td>
EMAIL_END
                    }
                    else
                    {
                        $email_data .= <<"EMAIL_END";
<tr bgcolor=$tdbg><td><center>
<a href=http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=refseqn&id=$uni_temp&style=raw>$uni_temp</a>
</center></td>
EMAIL_END
                    }
                }
                else
                {
                    if ( defined $hsp_pos{$uni_temp} )
                    {
                        $email_data .= <<"EMAIL_END";
<tr bgcolor=$tdbg><td><center>
<a href=http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=embl&id=$uni_temp&style=raw>$uni_temp</a> $hsp_pos{$uni_temp}
</center></td>
EMAIL_END
                    }
                    else
                    {
                        $email_data .= <<"EMAIL_END";
<tr bgcolor=$tdbg><td><center>
<a href=http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=embl&id=$uni_temp&style=raw>$uni_temp</a>
</center></td>
EMAIL_END
                    }
                }
            }
            $email_data .=
              sprintf '<td align=left><center>%5.2f</center></td><td></td>',
              $z_temp;
            if ( $one ne 'QUERY' )
            {
                if ( $db !~ /nt[.]fasta/ )
                {
                    if ( $common{$uni_temp} =~ /(\w+)\s(\w+)/ )
                    {
                        $email_data .= <<"EMAIL_END";
						<td title="$textcommon{$uni_temp}"><center>$1</center></td>
						<td title="$textncommon{$uni_temp}"><center>$2</center></td></tr>
EMAIL_END
                    }
                }
            }
            else
            {
                if ( $db !~ /nt[.]fasta/ )
                {
                    if ( $texts{$uni_temp} =~ /(\w+)/ )
                    {
                        $email_data .= <<"EMAIL_END";
						<td title="$texts2{$uni_temp}"><center>$1</center></td></tr>
EMAIL_END
                    }
                }
            }
        }
    }

    #############
    #Drawn tree.#
    #############
    $email_data .= <<"EMAIL_END";
        </table><img src='http://orion.mbg.duth.gr/Pinda/results/trees/drawn/$prid.png'><br>
        <a href=http://orion.mbg.duth.gr/Pinda/results/trees/zips/$prid.zip>
        NJ trees produced in .ph/.phb format
        </a>
EMAIL_END

}
else
{
    if ( $sequences_number == 3 )
    {
        ##########################################
        #Alignment, NJ-tree plotting and parsing.#
        ##########################################
        align( $sequences, $fnaln, $db, 1 );
        if ( $db !~ /nt[.]fasta/ && $masking == 1 )
        {
            my $conf_val = '/var/www/Pinda/results/final_alns/multalign/conf/' 
              . $prid . '.tmp';
            alignment_masking( $fnaln, $conf_val ) == 0 or die $?;
        }
        my $fnaln2 = $fnaln . ".fasta";
        if ( !( -e $fnaln2 ) )
        {
            system("cp $fnaln $fnaln2");
        }
        system("clustalw -INFILE=$fnaln2 -OUTFILE=$ph -tree") == 0 or die $?;
        rename "../results/final_alns/multalign/$prid.aln.ph",
          "../results/trees/phs/$prid.ph"
          or die $!;
        tree_manipulation1($ph);
        ####################################
        #Parsing the ph tree for distances.#
        ####################################
        system("../Pinda.R -parser $ph > /var/www/Pinda/parsing/$prid.tmp") == 0
          or die $?;
        system(
            "../Pinda.R -lengths_1 $ph > /var/www/Pinda/parsing/$prid\_1.tmp")
          == 0
          or die $?;
        system(
            "../Pinda.R -lengths_2 $ph > /var/www/Pinda/parsing/$prid\_2.tmp")
          == 0
          or die $?;

        tree_manipulation2($ph);

        system("zip -j $zip $ph") == 0 or die $?;
        ################
        #Draw the tree.#
        ################
        system("../Pinda.R $drawntree $ph") == 0
          or die $?;
        system("rm $ph") == 0 or die $?;

        $email_data .= <<"EMAIL_END";
			<font size='3' face='Georgia' color='330033'><br><br>
            <center>
            The number of sequences for this particular organism is three.
            <br><b>The dendrogram cannot be bootstrapped.</b></font>
            </center>
EMAIL_END

        parser( "../parsing/$prid.tmp", "../parsing/$prid\_1.tmp",
            "../parsing/$prid\_2.tmp" );
        compute_probability(@candidate);
        my $realign_num = @realign;
        if ( $realign_num >= 1 )
        {
            open $sequences_fh, '>', $sequences or die $!;
            my $one_dbfetch;
            if ( $db =~ /nt[.]fasta/ )
            {
                if ( $one =~ /^\D{2}\_/ )
                {
                    $one_dbfetch = get(
"http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=refseqn&id=$one&format=fasta&style=raw"
                    );
                    if ( !( defined $one_dbfetch ) )
                    {
                        for ( 0 .. 3 )
                        {
                            $one_dbfetch = get(
"http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=refseqn&id=$one&format=fasta&style=raw"
                            );
                            if ( defined $one_dbfetch )
                            {
                                last;
                            }
                        }
                    }
                }
                elsif ( $one eq 'QUERY' )
                {
                    $one_dbfetch = ">IIIQUERYIII\n$query_line\n\n";
                }
                else
                {
                    $one_dbfetch = get(
"http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=embl&id=$one&format=fasta&style=raw"
                    );
                    if ( !( defined $one_dbfetch ) )
                    {
                        for ( 0 .. 3 )
                        {
                            $one_dbfetch = get(
"http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=embl&id=$one&format=fasta&style=raw"
                            );
                            if ( defined $one_dbfetch )
                            {
                                last;
                            }
                        }
                    }
                    if ( !( defined $one_dbfetch ) )
                    {
                        for ( 0 .. 3 )
                        {
                            $one_dbfetch = get(
"http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=embl&id=$one&format=fasta&style=raw"
                            );
                            if ( defined $one_dbfetch )
                            {
                                last;
                            }
                        }
                    }
                }
                if ( $one_dbfetch =~ /\n/ )
                {
                    $one_dbfetch = ">III" . $one . "III\n" . $' . "\n\n";
                }
            }
            else
            {
                if ( $one eq 'QUERY' )
                {
                    $one_dbfetch = ">QUERY\n$query_line\n\n";
                }
                else
                {
                    $one_dbfetch =
                      get("http://www.uniprot.org/uniprot/$one.fasta");
                }
            }
            $one_dbfetch = uc $one_dbfetch;
            if ( $db !~ /nt[.]fasta/ )
            {
                $one_dbfetch =~ s/$one/***$one***/;
            }
            say {$sequences_fh} $one_dbfetch;
            my $dbfetch;
            foreach my $reseq (@realign)
            {
                if ( defined $hsp_seq{$reseq} && defined $hsp_pos{$reseq} )
                {
                    $dbfetch =
                      ">$reseq $hsp_pos{$reseq}\n" . $hsp_seq{$reseq} . "\n\n";
                }
                else
                {
                    if ( $db =~ /nt[.]fasta/ )
                    {
                        if ( $reseq =~ /^\D{2}\_/ )
                        {
                            $dbfetch = get(
"http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=refseqn&id=$reseq&format=fasta&style=raw"
                            );
                            if ( !( defined $dbfetch ) )
                            {
                                for ( 0 .. 3 )
                                {
                                    $dbfetch = get(
"http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=refseqn&id=$reseq&format=fasta&style=raw"
                                    );
                                    if ( defined $dbfetch )
                                    {
                                        last;
                                    }
                                }
                            }
                            if ( $dbfetch =~ /\n/ )
                            {
                                $dbfetch = ">$reseq\n" . $' . "\n\n";
                            }
                        }
                        else
                        {
                            $dbfetch = get(
"http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=embl&id=$reseq&format=fasta&style=raw"
                            );
                            if ( !( defined $dbfetch ) )
                            {
                                for ( 0 .. 3 )
                                {
                                    $dbfetch = get(
"http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=embl&id=$reseq&format=fasta&style=raw"
                                    );
                                    if ( defined $dbfetch )
                                    {
                                        last;
                                    }
                                }
                            }
                            if ( $dbfetch =~ /\n/ )
                            {
                                $dbfetch = ">$reseq\n" . $' . "\n\n";
                            }
                        }
                    }
                    else
                    {
                        $dbfetch =
                          get("http://www.uniprot.org/uniprot/$reseq.fasta");
                    }
                    $dbfetch = uc $dbfetch;
                }
                say {$sequences_fh} $dbfetch;
            }
            close $sequences_fh or die $!;
            $fnaln .= "2";
            align( $sequences, $fnaln, $db );
        }
        if ( $db !~ /nt[.]fasta/ )
        {
            ##############################
            #Find common gene ontologies.#
            ##############################
            GOs_in_common();
        }

        $email_data .= <<"EMAIL_END";
			<center><br><font size='3' face='Georgia' color='330033'>
			<a href=http://orion.mbg.duth.gr/Pinda/results/final_alns/multalign/$prid.aln>Alignment of all hits</a>
EMAIL_END
        my $clu_file =
          "/var/www/Pinda/results/final_alns/multalign/$prid.aln.clu";
        if ( -e $clu_file )
        {
            $email_data .= <<"EMAIL_END";
			| <a href=http://orion.mbg.duth.gr/Pinda/results/final_alns/multalign/$prid.aln.clu>Masked Alignment</a>
EMAIL_END
        }
        if ( $realign_num >= 1 )
        {
            $email_data .= <<"EMAIL_END";
         | <a href=http://orion.mbg.duth.gr/Pinda/results/final_alns/multalign/$prid.aln2>Alignment of top hits</a>
EMAIL_END
        }

        $email_data .= <<"EMAIL_END";
            </font>
            <br><br>
            <table border='1'>
            <tr bgcolor=FFFF66><th><center>Possible duplications of
EMAIL_END

        if ( $one ne 'QUERY' )
        {
            if ( $db !~ /nt[.]fasta/ )
            {
                $email_data .= <<"EMAIL_END";
            <a href=http://www.uniprot.org/uniprot/$one>
            $one</a>.</center></th>
            <th><center>Z Value</center></th>
            <th><center>Level Of Confidence</center></th>
            <th>GOs in common (out of $one_gos)</th>
            <th>GOs not found in $one</th></tr>
EMAIL_END
            }
            else
            {
                if ( $one =~ /^\D{2}\_/ )
                {
                    $email_data .= <<"EMAIL_END";
                <a href=http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=refseqn&id=$one&style=raw>
                $one</a>.</center></th>
                <th><center>Z Value</center></th>
				<th><center>Level Of Confidence</center></th></tr>
EMAIL_END
                }
                else
                {
                    $email_data .= <<"EMAIL_END";
            <a href=http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=embl&id=$one&style=raw>
            $one</a>.</center></th>
            <th><center>Z Value</center></th>
            <th><center>Level Of Confidence</center></th></tr>
EMAIL_END
                }
            }
        }
        else
        {
            if ( $db !~ /nt[.]fasta/ )
            {
                $email_data .= <<"EMAIL_END";
            $one.</center></th>
            <th><center>Z Value</center></th>
            <th><center>Level Of Confidence</center></th>
            <th>GO terms</th></tr>
EMAIL_END
            }
            else
            {
                if ( $one =~ /^\D{2}\_/ )
                {
                    $email_data .= <<"EMAIL_END";
                $one.</center></th>
                <th><center>Z Value</center></th>
				<th><center>Level Of Confidence</center></th></tr>
EMAIL_END
                }
                else
                {
                    $email_data .= <<"EMAIL_END";
                $one.</center></th>
                <th><center>Z Value</center></th>
				<th><center>Level Of Confidence</center></th></tr>
EMAIL_END
                }
            }
        }

        my $top_can;
        #################
        #Input sequence.#
        #################
        if ( $candidate[0] =~ /(\d?\d?\d?.?\d+e?-?\d*) \w+/ )
        {
            $top_can = $1;
        }

        foreach my $can (@candidate)
        {
            if (   $can !~ /$starting_point/
                && $can =~
                /(-?\d?\d?\d?.?\d+e?-?\d*) (\d?\d?\d?.?\d+e?-?\d*) (\w+)/ )
            {
                my $z_temp    = $1;
                my $conf_temp = $2;
                my $uni_temp  = $3;
                ###################################
                #Color the html table alternately.#
                ###################################
                if ( $tdcounter == 1 )
                {
                    $tdbg      = 'F8FBFE';
                    $tdcounter = 0;
                }
                else
                {
                    $tdbg = 'EAF1FB';
                    $tdcounter++;
                }
                if ( $db !~ /nt[.]fasta/ )
                {
                    if ( defined $hsp_pos{$3} )
                    {
                        $email_data .= <<"EMAIL_END";
<tr bgcolor=$tdbg><td><center>
<a href=http://www.uniprot.org/uniprot/$3>$3</a> $hsp_pos{$3}
</center></td>
EMAIL_END
                    }
                    else
                    {
                        $email_data .= <<"EMAIL_END";
<tr bgcolor=$tdbg><td><center>
<a href=http://www.uniprot.org/uniprot/$3>$3</a>
</center></td>
EMAIL_END
                    }
                }
                else
                {
                    if ( $uni_temp =~ /^\D{2}\_/ )
                    {
                        if ( defined $hsp_pos{$uni_temp} )
                        {
                            $email_data .= <<"EMAIL_END";
<tr bgcolor=$tdbg><td><center>
<a href=http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=refseqn&id=$uni_temp&style=raw>$uni_temp</a> $hsp_pos{$uni_temp}
</center></td>
EMAIL_END
                        }
                        else
                        {
                            $email_data .= <<"EMAIL_END";
<tr bgcolor=$tdbg><td><center>
<a href=http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=refseqn&id=$uni_temp&style=raw>$uni_temp</a>
</center></td>
EMAIL_END
                        }
                    }
                    else
                    {
                        if ( defined $hsp_pos{$uni_temp} )
                        {
                            $email_data .= <<"EMAIL_END";
<tr bgcolor=$tdbg><td><center>
<a href=http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=embl&id=$uni_temp&style=raw>$uni_temp</a> $hsp_pos{$uni_temp}
</center></td>
EMAIL_END
                        }
                        else
                        {
                            $email_data .= <<"EMAIL_END";
<tr bgcolor=$tdbg><td><center>
<a href=http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=embl&id=$uni_temp&style=raw>$uni_temp</a>
</center></td>
EMAIL_END
                        }
                    }
                }
                $email_data .=
                  sprintf
'<td align=left><center>%5.2f</center></td><td align=left><center>%5.1f%%</center></td>',
                  $z_temp, $conf_temp;
                if ( $one ne 'QUERY' )
                {
                    if ( $db !~ /nt[.]fasta/ )
                    {
                        if ( $common{$uni_temp} =~ /(\w+)\s(\w+)/ )
                        {
                            $email_data .= <<"EMAIL_END";
						<td title="$textcommon{$uni_temp}"><center>$1</center></td>
						<td title="$textncommon{$uni_temp}"><center>$2</center></td></tr>
EMAIL_END
                        }
                    }
                }
                else
                {
                    if ( $db !~ /nt[.]fasta/ )
                    {
                        if ( $texts{$uni_temp} =~ /(\w+)/ )
                        {
                            $email_data .= <<"EMAIL_END";
						<td title="$texts2{$uni_temp}"><center>$1</center></td></tr>
EMAIL_END
                        }
                    }
                }
            }
            elsif ($can !~ /$starting_point/
                && $can =~ /(-?\d?\d?\d?.?\d+e?-?\d*) (\w+)/ )
            {
                my $z_temp   = $1;
                my $uni_temp = $2;
                ###################################
                #Color the html table alternately.#
                ###################################
                if ( $tdcounter == 1 )
                {
                    $tdbg      = 'F8FBFE';
                    $tdcounter = 0;
                }
                else
                {
                    $tdbg = 'EAF1FB';
                    $tdcounter++;
                }
                if ( $db !~ /nt[.]fasta/ )
                {
                    if ( defined $hsp_pos{$2} )
                    {
                        $email_data .= <<"EMAIL_END";
<tr bgcolor=$tdbg><td><center>
<a href=http://www.uniprot.org/uniprot/$2>$2</a> $hsp_pos{$2}
</center></td>
EMAIL_END
                    }
                    else
                    {
                        $email_data .= <<"EMAIL_END";
<tr bgcolor=$tdbg><td><center>
<a href=http://www.uniprot.org/uniprot/$2>$2</a>
</center></td>
EMAIL_END
                    }
                }
                else
                {
                    if ( $uni_temp =~ /^\D{2}\_/ )
                    {
                        if ( defined $hsp_pos{$uni_temp} )
                        {
                            $email_data .= <<"EMAIL_END";
<tr bgcolor=$tdbg><td><center>
<a href=http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=refseqn&id=$uni_temp&style=raw>$uni_temp</a> $hsp_pos{$uni_temp}
</center></td>
EMAIL_END
                        }
                        else
                        {
                            $email_data .= <<"EMAIL_END";
<tr bgcolor=$tdbg><td><center>
<a href=http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=refseqn&id=$uni_temp&style=raw>$uni_temp</a>
</center></td>
EMAIL_END
                        }
                    }
                    else
                    {
                        if ( defined $hsp_pos{$uni_temp} )
                        {
                            $email_data .= <<"EMAIL_END";
<tr bgcolor=$tdbg><td><center>
<a href=http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=embl&id=$uni_temp&style=raw>$uni_temp</a> $hsp_pos{$uni_temp}
</center></td>
EMAIL_END
                        }
                        else
                        {
                            $email_data .= <<"EMAIL_END";
<tr bgcolor=$tdbg><td><center>
<a href=http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=embl&id=$uni_temp&style=raw>$uni_temp</a>
</center></td>
EMAIL_END
                        }
                    }
                }
                $email_data .=
                  sprintf '<td align=left><center>%5.2f</center></td><td></td>',
                  $z_temp;
                if ( $one ne 'QUERY' )
                {
                    if ( $db !~ /nt[.]fasta/ )
                    {
                        if ( $common{$uni_temp} =~ /(\w+)\s(\w+)/ )
                        {
                            $email_data .= <<"EMAIL_END";
						<td title="$textcommon{$uni_temp}"><center>$1</center></td>
						<td title="$textncommon{$uni_temp}"><center>$2</center></td></tr>
EMAIL_END
                        }
                    }
                }
                else
                {
                    if ( $db !~ /nt[.]fasta/ )
                    {
                        if ( $texts{$uni_temp} =~ /(\w+)/ )
                        {
                            $email_data .= <<"EMAIL_END";
						<td title="$texts2{$uni_temp}"><center>$1</center></td></tr>
EMAIL_END
                        }
                    }
                }
            }
        }

        $email_data .= <<"EMAIL_END";
            </table><center>
            <img src='http://orion.mbg.duth.gr/Pinda/results/trees/drawn/$prid.png'>
            <br>
            <a href=http://orion.mbg.duth.gr/Pinda/results/trees/zips/$prid.zip>
            NJ trees produced in .ph/.phb format
            </a></center>
EMAIL_END

    }
    if ( $sequences_number == 2 )
    {
        #################
        #Alignment only.#
        #################
        align( $sequences, $fnaln, $db );
        $email_data .= <<"EMAIL_END";
            <center>
            <font size='3' face='Georgia' color='330033'><br><br>
            Only <b>two</b> similar sequences from <b>$organism</b>
            have been identified.
            <br><b>Phylogenetic analysis is meaningless therefore.</b>
            </font>
            </center>
EMAIL_END

        $email_data .= <<"EMAIL_END";
            <center><br><font size='3' face='Georgia' color='330033'><br>
                <a href=http://orion.mbg.duth.gr/Pinda/results/final_alns/multalign/$prid.aln>
                Alignment</a></font></center>
EMAIL_END
    }
    if ( $sequences_number <= 1 )
    {
        alarm 0;
        #####################
        #Nothing to do here.#
        #####################
        if ( $db =~ /nt[.]fasta/ )
        {
            $email_data .= <<"EMAIL_END";
            <center>
            <font size='3' face='Georgia' color='330033'><br><br>
            No similar sequences from <b>$organism</b>
            have been identified.
            <br><b>Phylogenetic analysis is meaningless therefore.</b>
            </font>
            </center>
EMAIL_END
        }
        else
        {
            $email_data .= <<"EMAIL_END";
            <center>
            <font size='3' face='Georgia' color='330033'><br><br>
            No similar sequences from <b>$organism</b>
            have been identified.
            <br><b>Phylogenetic analysis is meaningless therefore.</b>
            </font>
            </center>
EMAIL_END
        }
    }
}
##############################################
#Calculate how much time it took for the job.#
##############################################
my $end_timer   = time;
my $run_time    = $end_timer - $start_timer;
my $job_average = "/var/www/Pinda/job_times";
my ( $pr_jobs, $pr_time, $dn_jobs, $dn_time );
open my $job_average_fh, '<', $job_average;
local $/ = "\n";
while ( my $line = <$job_average_fh> )
{

    if ( $line =~ /Protein Jobs[:] (\d+) Average Time[:] (\d+)/ )
    {
        $pr_jobs = $1;
        $pr_time = $2;
    }
    elsif ( $line =~ /DNA Jobs[:] (\d+) Average Time[:] (\d+)/ )
    {
        $dn_jobs = $1;
        $dn_time = $2;
    }
}
close $job_average_fh;

if ( $db =~ /nt[.]fasta/ )
{
    $dn_time = ( ( $dn_jobs * $dn_time ) + $run_time ) / ( $dn_jobs + 1 );
    $dn_jobs++;
}
else
{
    $pr_time = ( ( $pr_jobs * $pr_time ) + $run_time ) / ( $pr_jobs + 1 );
    $pr_jobs++;
}
open $job_average_fh, '>', $job_average;
say {$job_average_fh} "Protein Jobs: $pr_jobs Average Time: $pr_time";
say {$job_average_fh} "DNA Jobs: $dn_jobs Average Time: $dn_time";
close $job_average_fh;

$email_data .= <<"ENDHTML";
<br><br>
<font size='1'>Temporary files from this job, including alignments, .ph/.phb
trees and plotted trees, will be <b>DELETED</b> after ten days.</font>
</center>
</body>
</html>
ENDHTML

##############
#Send e-mail.#
##############
send_email( $one, $email );

my $job_counting = "/var/www/Pinda/running_jobs";
my $protein_jobs;
my $dna_jobs;
open my $job_counting_fh, '<', $job_counting;
local $/ = "\n";
while ( my $line = <$job_counting_fh> )
{

    if ( $line =~ /Protein[:] (\d+)/ )
    {
        $protein_jobs = $1;
    }
    elsif ( $line =~ /DNA[:] (\d+)/ )
    {
        $dna_jobs = $1;
    }
}
close $job_counting_fh;
if ( $db =~ /nt[.]fasta/ )
{
    $dna_jobs--;
}
else
{
    $protein_jobs--;
}
open $job_counting_fh, '>', $job_counting;
say {$job_counting_fh} "Protein: $protein_jobs";
say {$job_counting_fh} "DNA: $dna_jobs";
close $job_counting_fh;
system("rm /var/www/Pinda/slurm_errors/slurm-$slurm_id.out");

#############################
#Multiple Sequence Alignment#
#############################
sub align
{
    my $out = $_[1] . ".fasta";
    my $db  = $_[2];
    my $again;
    if ( defined $_[3] )
    {
        $again = $_[3];
    }
    if ( $db !~ /nt[.]fasta/ && defined $again && $masking == 1 )
    {
        system(
"/usr/local/bin/clustalo -i $_[0] -o $_[1] --outfmt=clu --threads=4 -v --force"
        );
        system(
"/usr/local/bin/clustalo -i $_[0] -o $out --outfmt=fasta --threads=4 -v --force"
        );
    }
    elsif (( $db !~ /nt[.]fasta/ && !( defined $again ) )
        || ( $db !~ /nt[.]fasta/ && $masking == 0 ) )
    {
        system(
"/usr/local/bin/clustalo -i $_[0] -o $_[1] --outfmt=clu --threads=4 -v --force"
        );
    }
    else
    {
        system("/usr/local/bin/kalign -i $_[0] -o $_[1] -f clu -q");
        my @sequences2  = ();
        my $seq2counter = 0;
        open my $fnaln_fh, '<', $fnaln or die $!;
        my $line2;
        {
            local $/ = "\n";
            while ( $line2 = <$fnaln_fh> )
            {
                ##################################
                #Add stars to the input sequence.#
                ##################################
                my $one_temp = "III" . $one . "III";
                $line2 =~ s/$one_temp/***$one***/;

                $sequences2[$seq2counter] = $line2;
                $seq2counter++;
            }
        }
        close $fnaln_fh or die $!;

        open $fnaln_fh, '>', $fnaln or die $!;
        foreach my $seq2counter (@sequences2)
        {
            print {$fnaln_fh} $seq2counter;
        }
        close $fnaln_fh or die $!;

        open $fnaln_fh, '<', $fnaln;
        local $/ = undef;
        my $line = q{};
        while ( $line = <$fnaln_fh> )
        {
            if ( $line =~ /in ClustalW format/ )
            {
                $line = $';
                last;
            }
        }
        close $fnaln_fh;
        if ( defined $line )
        {
            open $fnaln_fh, '>', $fnaln;
            print {$fnaln_fh} "CLUSTAL W (1.83) multiple sequence alignment";
            print {$fnaln_fh} $line;
            close $fnaln_fh;
        }
    }
    return 0;
}

############################################################
#Mask alignments using confidence values generated by ZORRO#
############################################################
sub alignment_masking
{
    my $input  = $_[0] . ".fasta";
    my $output = $_[1];
    system("/usr/local/bin/zorro $input > $output");

    my $counter  = 0;
    my $counter2 = 0;
    my @low_conf;

    open my $output_fh, '<', $output or die $!;
    local $/ = "\n";
    while ( my $line = <$output_fh> )
    {
        if ( $line =~ /(\d+[.]\d+)/ )
        {
            if ( $1 <= 0.4 )
            {
                $low_conf[$counter2] = $counter + 1;
                $counter2++;
            }
        }
        $counter++;
    }
    close $output_fh;

    if ( $counter2 >= 1 && ( $counter2 / $counter ) < 0.75 )
    {
        open my $input_fh, '<', $input or die $!;
        local $/ = undef;
        my $newline = <$input_fh>;
        close $input_fh;
        my $seqc = 0;
        my $resc = 0;
        my @sequence;

        while ( $newline =~ />[*]?[*]?[*]?\w+/ )
        {
            if ( $newline =~ /(>.+\s)/ )
            {
                $sequence[$seqc][0] = $1;
                if ( $' =~ />/ || $' =~ /\s$/ )
                {
                    $resc = 0;
                    my $temp = $`;
                    $newline = '>' . $';
                    $temp =~ s/\n//g;
                    foreach my $residue ( split //, $temp )
                    {
                        $resc++;
                        $sequence[$seqc][$resc] = $residue;
                    }
                }
                else
                {
                    last;
                }
                $seqc++;
            }
        }

        foreach my $position (@low_conf)
        {
            for my $seqtmp ( 0 .. $seqc - 1 )
            {
                undef $sequence[$seqtmp][$position];
            }
        }

        open $input_fh, '>', $input or die $!;
        for my $seqtmp ( 0 .. $seqc - 1 )
        {
            print {$input_fh} $sequence[$seqtmp][0];
            for ( 1 .. $resc )
            {
                if ( defined $sequence[$seqtmp][$_] )
                {
                    print {$input_fh} $sequence[$seqtmp][$_];
                }
            }
            say {$input_fh} q{};
        }
        close $input_fh;
        my $clu_file = $_[0] . ".clu";
        system("sreformat clustal $input > $clu_file");
    }
    elsif ( $counter2 == 0 )
    {
        $email_data .= <<"ENDHTML";
		<br><center>All this alignment's columns have high confidence values.
		<br>Masking is not needed.</center>
ENDHTML
    }
    else
    {
        $email_data .= <<"ENDHTML";
		<br><center>Most of this alignment's columns have poor confidence values.
		<br>Masking was <b>NOT</b> performed.</center>
ENDHTML
        system("cp $_[0] $input");
    }
    return 0;
}

################################################################
#Removing the TRICHOTOMY word, setting negative values to zero.#
################################################################
sub tree_manipulation1
{
    my $yacounter = 0;
    my ( $tree_fh, @tree );
    open $tree_fh, '<', $_[0] or die $!;
    my $line;
    {
        local $/ = "\n";
        while ( $line = <$tree_fh> )
        {
            $line =~ s/TRICHOTOMY//;
            $line =~ s/-\d[.]\d*/0/;
 ###############################################################################
 #Fix the NJ negative branch length artifact, by setting those numbers to zero.#
 ###############################################################################
            $tree[$yacounter] = $line;
            $yacounter++;
        }
    }
    close $tree_fh or die $!;
    open $tree_fh, '>', $_[0] or die $!;
    foreach my $yacounter (@tree)
    {
        print {$tree_fh} $yacounter;
    }
    close $tree_fh or die $!;
    return 0;
}
#######################################################################
#Dividing the bootstrap values by 10, to reach a maximum value of 100.#
#######################################################################
sub tree_manipulation2
{
    my $yacounter = 0;
    my ( $bvdiv, $tree_fh );
    my @tree = ();
    open $tree_fh, '<', $_[0] or die $!;
    my $line;
    {
        local $/ = "\n";
        while ( $line = <$tree_fh> )
        {
            if ( $line =~ /^(\d+):/ )
            {
                $bvdiv = $1 / 10.0;
                $line  = $` . $bvdiv . q{:} . $';
            }
            $tree[$yacounter] = $line;
            $yacounter++;
        }
    }
    close $tree_fh or die $!;
    open $tree_fh, '>', $_[0] or die $!;
    foreach my $yacounter (@tree)
    {
        print {$tree_fh} $yacounter;
    }
    close $tree_fh or die $!;
    return 0;
}

##########################################################
#Parse the tree to exract distances and bootstrap values.#
##########################################################
sub parser
{
    my $cancounter    = 0;
    my $ginomenon     = 1;
    my $linocounter   = 0;
    my $node_distance = 0;
    my $plc           = 0;
    my $plc2          = 0;
    my $plcn          = 0;
    my $plcn2         = 0;
    my $sscounter     = 0;
    my $unicounter    = 0;

    my (
        $continue,             $degree_of_confidence, $last_yes,
        $search_id,            $tree_for_parsingfh,   $tree_node_distancesfh,
        $tree_tip_distancesfh, @compare_seq,          @parsed_lines,
        @parsed_linesnew,      @parsing_lines,        @star_seq,
        @uni_ids,              %distanceshash,        %distanceshash2
    );

    open $tree_for_parsingfh, '<', $_[0] or die $!;
    my $line;
    {
        local $/ = "\n\n";
        ########################
        #Get all the sequences.#
        ########################
        while ( $line = <$tree_for_parsingfh> )
        {
            ##########################
            #Find the input sequence.#
            ##########################
            if ( $line =~ /\s"(\w\_\_\_\w+\_\_\_)"/ )
            {
                $starting_point = $1;
                $uni_ids[$unicounter] = $1;
                $unicounter++;
            }
            ###############################
            #Find all the other sequences.#
            ###############################
            if ( $line =~ /"(\w{6,})"/ )
            {
                if ( $1 !~ /\D\d{1,4}\_/ )
                {
                    $uni_ids[$unicounter] = $1;
                    $unicounter++;
                }
            }
            if ( $' =~ /"(\w{6,})"/ && $1 !~ /\D\d{1,4}\_/ )
            {
                $uni_ids[$unicounter] = $1;
                $unicounter++;
            }
            $parsing_lines[$linocounter] = $line;
            $linocounter++;
        }
    }

    #################################
    #Parsing distance between nodes.#
    #################################
    open $tree_node_distancesfh, '<', $_[1]
      or die $!;
    my $neoline;
    {
        local $/ = "\n";
        while ( $neoline = <$tree_node_distancesfh> )
        {
            if ( $neoline !~ /\$/ && $neoline =~ /\w/ )
            {
                do
                {
                    ####################################
                    #This line contains sequence names.#
                    ####################################
                    if ( $neoline !~ /[.]/ )
                    {
                        if (   $neoline =~ /\s?(\w+)/
                            || $neoline =~ /(Root)/ )
                        {
                            $parsed_lines[$plc][$plc2] = $1;
                            $plc2++;
                            $continue = $';
                            if ( $continue =~ /\w/ )
                            {
                                $neoline = $continue;
                            }
                            else
                            {
                                $plc++;
                            }
                        }
                    }
                    elsif (
                        ############################
                        #This line contains values.#
                        ############################
                        $neoline =~ /(\d+[.]\d+)/
                      )
                    {
                        $distanceshash{ $parsed_lines[ $plc - 1 ][$plc2] } = $1;
                        $plc2++;
                        $continue = $';
                        if ( $continue =~ /\w/ )
                        {
                            $neoline = $continue;
                        }
                    }
                } while ( $continue =~ /\w/ );
                $plc2 = 0;
            }
        }
    }
    close $tree_node_distancesfh or die $!;
    ######################################################
    #Parsing distance from every tip to its closest node.#
    ######################################################
    open $tree_tip_distancesfh, '<', $_[2]
      or die $!;
    my $neoline2;
    {
        local $/ = "\n";
        while ( $neoline2 = <$tree_tip_distancesfh> )
        {
            if ( $neoline2 !~ /\$/ && $neoline2 =~ /\w/ )
            {
                do
                {
                    ####################################
                    #This line contains sequence names.#
                    ####################################
                    if ( $neoline2 !~ /[.]/ )
                    {
                        if ( $neoline2 =~ /\_?\_?\_?(\w{6,})\_?\_?\_?\s/ )
                        {
                            $parsed_linesnew[$plcn][$plcn2] = $1;
                            $plcn2++;
                            $continue = $';
                            if ( $continue =~ /\w/ )
                            {
                                $neoline2 = $continue;
                            }
                            else
                            {
                                $plcn++;
                            }
                        }
                    }
                    ############################
                    #This line contains values.#
                    ############################
                    elsif ( $neoline2 =~ /(\d+[.]\d+)/ )
                    {
                        $distanceshash2{ $parsed_linesnew[ $plcn - 1 ][$plcn2] }
                          = $1;
                        $plcn2++;
                        $continue = $';
                        if ( $continue =~ /\w/ )
                        {
                            $neoline2 = $continue;
                        }
                    }
                } while ( $continue =~ /\w/ );
                $plcn2 = 0;
            }
        }
    }
    close $tree_tip_distancesfh or die $!;

    my $p_l_s = @parsing_lines;
    foreach my $line (@parsing_lines)
    {
        if ( $line =~ /$starting_point/ )
        {
            #####################################
            #Get the node of the input sequence.#
            #####################################
            if ( $line =~ /parts\$(\w)(\w+)/ )
            {
                $star_seq[$sscounter] = $1 . $2;
                $search_id = $1 . $2;
                $sscounter++;
            }
            for ( 0 .. $p_l_s - 1 )
            {
                if ( $parsing_lines[$_] =~ /$search_id"/ )
                {
                    #####################################
                    #Get the node leading to the node...#
                    #####################################
                    if ( $parsing_lines[$_] =~ /parts\$(\w)(\w+)/ )
                    {
                        $star_seq[$sscounter] = $1 . $2;
                        $search_id = $1 . $2;
                        $sscounter++;
                        $_--;
                    }
                }
            }
        }
    }
    $last_yes = 0;
    foreach my $uni (@uni_ids)
    {
        $sscounter = 0;
        foreach my $line (@parsing_lines)
        {
            if ( $line =~ /$uni/ )
            {
                ####################################################
                #Get the closest node for  all the other sequences.#
                ####################################################
                if ( $line =~ /parts\$(\w)(\w+)/ )
                {
                    $compare_seq[$sscounter] = $1 . $2;
                    $search_id = $1 . $2;
                    $sscounter++;
                }
                for ( 0 .. $p_l_s - 1 )
                {
                    if ( $parsing_lines[$_] =~ /$search_id"/ )
                    {
                        #####################################
                        #Get the node leading to the node...#
                        #####################################
                        if ( $parsing_lines[$_] =~ /parts\$(\w)(\w+)/ )
                        {
                            $compare_seq[$sscounter] = $1 . $2;
                            $search_id = $1 . $2;
                            $sscounter++;
                            $_--;
                        }
                    }
                }
            }
        }
        foreach my $node (@compare_seq)
        {
            if ( $node =~ /(\d+)\_?/ )
            {
                #################################
                #Add the distance for this node.#
                #################################
                $node_distance += $distanceshash{"$node"};
                $ginomenon *= $1 / 1000.0;
            }
            foreach my $star_node (@star_seq)
            {
              ##################################################################
              #Find the common node for the input sequence and every other one.#
              ##################################################################
                if ( $node eq $star_node )
                {
                    foreach my $star_node (@star_seq)
                    {
                        ######################################
                        #Add distances until the common node.#
                        ######################################
                        if (   $star_node ne $node
                            && $star_node =~ /(\d+)\_?/ )
                        {
                            $node_distance += $distanceshash{"$star_node"};
                            $ginomenon *= $1 / 1000.0;
                        }
                        ####################################
                        #Remove the common node's distance.#
                        ####################################
                        if ( $star_node eq $node )
                        {
                            if ( defined $distanceshash{"$star_node"} )
                            {
                                $node_distance -= $distanceshash{"$star_node"};
                                last;
                            }
                            else
                            {
                                $node_distance = $distanceshash{"$star_node"};
                                last;
                            }
                        }
                    }
                    #################################################
                    #For every sequence, except for the input one...#
                    #################################################
                    if ( $uni !~ /$starting_point/ )
                    {
                        ########################################
                        #Add the tip-to-closest-node distances.#
                        ########################################
                        $node_distance +=
                          $distanceshash2{"$uni"} +
                          $distanceshash2{"$starting_point"};
                       #########################################################
                       #If distance is 0, then the sequences are overlapping!?!#
                       #########################################################
                        if ( $node_distance == 0 )
                        {
                            $degree_of_confidence = $ginomenon / 0.1;
                            $candidate[$cancounter] =
                              $degree_of_confidence . q{ } . $uni;
                            $cancounter++;
                        }
                        #################################
                        #Calculate the confidence value.#
                        #################################
                        else
                        {
                            $degree_of_confidence = $ginomenon / $node_distance;
                            $candidate[$cancounter] =
                              $degree_of_confidence . q{ } . $uni;
                            $cancounter++;
                        }
                    }
                    $last_yes = 1;
                    last;
                }
                if ( $last_yes == 1 )
                {
                    last;
                }
            }
            if ( $last_yes == 1 )
            {
                last;
            }
        }
        if ( $last_yes == 0 )
        {
            $node_distance = 0;
            $ginomenon     = 1;
            ###################################################################
            #If the input sequence had only the root common with a sequence...#
            ###################################################################
            foreach my $node (@compare_seq)
            {
                if ( $node =~ /(\d+)\_?/ )
                {
                    ####################################################
                    #Add the distance from the sequence's node to root.#
                    ####################################################
                    $node_distance += $distanceshash{"$node"};
                    $ginomenon *= $1 / 1000.0;
                }
            }
            foreach my $star_node (@star_seq)
            {
                if ( $star_node =~ /(\d+)\_?/ )
                {
                    ##########################################################
                    #Add the distance from the input sequence's node to root.#
                    ##########################################################
                    $node_distance += $distanceshash{"$star_node"};
                    $ginomenon *= $1 / 1000.0;
                }
            }
            ########################################
            #Add the tip-to-closest-node distances.#
            ########################################
            $node_distance +=
              $distanceshash2{"$uni"} + $distanceshash2{"$starting_point"};
            $degree_of_confidence = $ginomenon / $node_distance;
            $candidate[$cancounter] = $degree_of_confidence . q{ } . $uni;
            $cancounter++;
        }
        $node_distance = 0;
        $ginomenon     = 1;
        @compare_seq   = ();
        $last_yes      = 0;
    }
    @candidate = map join( q{ }, @{$_} ),
      sort { $a->[0] <=> $b->[0] } map { [split] } @candidate;
    ############################################################################
    #Sort the sequences in descending order, according to the confidence value.#
    ############################################################################
    @candidate = reverse @candidate;
    $cand_sans[0] = $one;
    my $neocounter = 1;
    foreach my $can (@candidate)
    {
        if ( $can =~ /\s(\w+)/ )
        {
            $cand_sans[$neocounter] = $1;
            $neocounter++;
        }
    }
    return @candidate, @cand_sans, $starting_point;
}

############################################################
#Calculate the Z-values and the probabilities for each hit.#
############################################################
sub compute_probability
{
    my @numbers;
    my $numcounter = 0;
    my $re         = 0;
    ############################
    #Get the Confidence values.#
    ############################
    foreach my $can (@candidate)
    {
        if (   $can !~ /$starting_point/
            && $can =~ /(\d?\d?\d?.?\d+e?-?\d*) \w+/ )
        {
            $numbers[$numcounter] = $1;
            $numcounter++;
        }
    }
    ############################################
    #Calculate the mean and standard deviation.#
    ############################################
    my $mean               = mean(@numbers);
    my $standard_deviation = stddev(@numbers);
    foreach my $can (@candidate)
    {
        if (   $can !~ /$starting_point/
            && $can =~ /(\d?\d?\d?.?\d+e?-?\d*) (\w+)/ )
        {
            #####################
            #Calculate Z-values.#
            #####################
            my $z = ( $1 - $mean ) / $standard_deviation;
            ###########################################################
            #For positive Z-values, calculate the level of confidence.#
            ###########################################################
            if ( $z > 0.674 )
            {
                my $erf  = erfc( $z / sqrt(2) );
                my $prob = ( 1 - $erf ) * 100;
                $can = $z . q{ } . $prob . q{ } . $2;
                $realign[$re] = $2;
                $re++;
            }
            else
            {
                $can = $z . q{ } . $2;
            }
        }
    }
    return @candidate, @realign;
}

##########################################################################
#Find the GOs shared between the input sequence and every other sequence.#
##########################################################################
sub GOs_in_common
{
    my @input_gos;
    my $input_counter = 0;
    if ( $one ne "QUERY" )
    {
        foreach my $seq (@cand_sans)
        {
            chomp $seq;
            if ( $seq eq $one )
            {
                my $dbfetch = get("http://www.uniprot.org/uniprot/$one.txt");
                if ( !( defined $dbfetch ) )
                {
                    for ( 0 .. 3 )
                    {
                        $dbfetch =
                          get("http://www.uniprot.org/uniprot/$one.txt");
                        if ( defined $dbfetch )
                        {
                            last;
                        }
                    }
                }
                do
                {
                    ############################################################
                    #Parse the input sequence's UniProt flat file for GO terms.#
                    ############################################################
                    if ( $dbfetch =~ /GO; GO:(\d+);/ )
                    {
                        $input_gos[$input_counter] = $1;
                        $input_counter++;
                        $dbfetch = $';
                    }
                } while ( $dbfetch =~ /GO; GO:(\d+);/ );
                $one_gos = @input_gos;
                if ( $one_gos == 0 )
                {
                    last;
                }
            }
            else
            {
                my $common_prop  = q{};
                my $ncommon_prop = q{};
                my @go_list      = ();
                my $ncpi         = 0;
                my $eq_counter   = 0;
                my $res_counter  = 0;
                my $dbfetch = get("http://www.uniprot.org/uniprot/$seq.txt");
                if ( !( defined $dbfetch ) )
                {

                    for ( 0 .. 3 )
                    {
                        $dbfetch =
                          get("http://www.uniprot.org/uniprot/$seq.txt");
                        if ( defined $dbfetch )
                        {
                            last;
                        }
                    }
                }
                do
                {
        ########################################################################
        #Parse every other sequence's UniProt flat file for GO terms in common.#
        ########################################################################
                    if ( defined $dbfetch
                        && $dbfetch =~ /GO; GO:(\d+);\s(.+);/ )
                    {
                        $go_list[$ncpi] = $2;
                        $ncpi++;
                        foreach my $go (@input_gos)
                        {
                            if ( $go == $1 )
                            {
                                $eq_counter++;
                                $common_prop .= $2 . "\n";
                            }
                        }
                        $res_counter++;
                        $dbfetch = $';
                    }
                } while ( $dbfetch =~ /GO; GO:(\d+);\s(.+);/ );
                foreach my $term (@go_list)
                {
                    $term =~ s/\?/\\?/;
                    if ( $common_prop !~ /$term/ )
                    {
                        $ncommon_prop .= $term . "\n";
                    }
                }
                my $neq_counter = $res_counter - $eq_counter;
                #######################################
                #If no GO terms are found, put in "NA"#
                #######################################
                if ( $eq_counter == 0 && $res_counter == 0 )
                {
                    $eq_counter  = "NA";
                    $neq_counter = "NA";
                }
                $common{$seq} = $eq_counter . q{ } . $neq_counter;
                chomp $common_prop;
                chomp $ncommon_prop;
                $textcommon{$seq}  = $common_prop;
                $textncommon{$seq} = $ncommon_prop;
            }
        }
        return $one_gos, %common, %textcommon, %textncommon;
    }
    else
    {
        foreach my $seq (@cand_sans)
        {
            chomp $seq;
            if ( $seq ne $one )
            {
                my $ontologies  = q{};
                my @go_list     = ();
                my $ncpi        = 0;
                my $res_counter = 0;
                my $dbfetch = get("http://www.uniprot.org/uniprot/$seq.txt");
                if ( !( defined $dbfetch ) )
                {
                    for ( 0 .. 3 )
                    {
                        $dbfetch =
                          get("http://www.uniprot.org/uniprot/$seq.txt");
                        if ( defined $dbfetch )
                        {
                            last;
                        }
                    }
                }
                do
                {
        ########################################################################
        #Parse every other sequence's UniProt flat file for GO terms in common.#
        ########################################################################
                    if ( defined $dbfetch
                        && $dbfetch =~ /GO; GO:(\d+);\s(.+);/ )
                    {
                        $go_list[$ncpi] = $2;
                        $ncpi++;
                        $res_counter++;
                        $dbfetch = $';
                    }
                  } while ( defined $dbfetch
                    && $dbfetch =~ /GO; GO:(\d+);\s(.+);/ );
                foreach my $term (@go_list)
                {
                    $term =~ s/\?/\\?/;
                    $ontologies .= $term . "\n";
                }
                #######################################
                #If no GO terms are found, put in "NA"#
                #######################################
                if ( $res_counter == 0 )
                {
                    $res_counter = "NA";
                }
                $texts{$seq} = $res_counter;
                chomp $ontologies;
                $texts2{$seq} = $ontologies;
            }
        }
        return %texts, %texts2;
    }
}

###########################
#Email sending subroutine.#
###########################
sub send_email
{
    my $msg = MIME::Lite->new(
        Subject => "Pinda Job Result: $_[0]",
        From    => 'Pinda@orion.mbg.duth.gr',
        To      => $_[1],
        Type    => 'text/html',
        Data    => $email_data
    );
    $msg->send();
    return 0;
}
