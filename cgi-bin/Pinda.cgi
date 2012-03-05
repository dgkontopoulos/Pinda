#!/usr/bin/perl -w

#$VERSION = '0.01';

=head1 NAME

PINDA - Pipeline for INtraspecies Duplication Analysis

=head1 DESCRIPTION

A CGI program aiming to facilitate detection of specific gene
duplications in an organism species of choice.

=head1 AUTHOR

Dimitrios - Georgios Kontopoulos

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
use CGI qw(:standard);
use Data::Validate::Email qw(is_email);
use FreezeThaw qw(freeze thaw);
use List::MoreUtils qw(uniq);
use LWP::Simple qw(get);
use MIME::Lite;
use Sys::CPU;

use strict;
use warnings;

my ( $starting_point, @candidate );

open STDERR, '>', '/dev/null' or die $!;

my $query = CGI->new;    #Start the CGI environment.
print $query->header;
print <<"ENDHTML";
<!DOCTYPE html
	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	 "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" lang="en-US" xml:lang="en-US">
<head>
<title>Pinda - Pipeline for Intraspecies Duplication Analysis</title>
<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1" /> 
<link rel="stylesheet" href="../css/Pinda.css" type="text/css" media="screen">

<script src="https://ajax.googleapis.com/ajax/libs/jquery/1.6.2/jquery.min.js"
type="text/javascript" charset="utf-8"></script>

<script src="../js/jquery.dropkick.1-0.1.js" type="text/javascript" 
charset="utf-8"></script><script type="text/javascript" charset="utf-8">
\$(function () {
    \$('.default').dropkick();
});
</script>

</head>
<body background='../background.jpg'>
<LINK REL='SHORTCUT ICON' HREF='../pinda.ico'>
<center><br><a href='$0'><img src='../pindalogo.png'></a><br>
<p style='width: 500px; text-align:center;margin-bottom:1px;margin-top:1px'>
<br>
<hr/>
</p>
</center>
ENDHTML

my $email_data = <<'EMAIL_END';
<center><br>
<a href='http://localhost/cgi-bin/Pinda.cgi'>
<img src='http://localhost/pindalogo.png'></a>
<br>
EMAIL_END

if ( !$query->param )    #Index Page with input box, DB selection and email
{
    print <<"ENDHTML";
    <center>
    <p style='width: 470px; text-align:left;margin-bottom:1px;margin-top:1px'>
    <br><font size='3' face='Georgia' color='330033'>
    <i>Please enter a sequence.</i>
    </font><br>
    <form id = 'submit' method = 'post' action=$0>
    <textarea name='sequence' rows=7 cols=55 style='font-family:courier new'
    ></textarea></p><br>
    <p style='width: 470px; text-align:left;margin-bottom:1px;margin-top:1px'>
    <font size='3' face='Georgia' color='330033'><i>Please choose a database.
    </i></font>
    <div id="container">
        <ul class="menu">
            <li id="DNA" class="active">DNA</li>
            <li id="Protein">Protein</li>
        <ul>
        <span class="clear"></span>
        <div class="content DNA">
            <ul>
                <font side='2'>
                <input type="radio" name="db" value="nt" checked> <b>nt</b> 
                (<i>NCBI</i>)
                </font>
                
            <ul>
        </div>
        <div class="content Protein">
            <ul>
                <font size='2'>
                <input type="radio" name="db" value="Swiss-Prot" checked> <b>Swiss-Prot</b>
                <br>
                <input type="radio" name="db" value="UniProt"> <b>UniProt</b> (<i>Swiss-Prot
                + TrEMBL</i>)
                </font>
            <ul>
        </div>
    </div>
    <script type="text/javascript" src="../js/tabs.js"></script>
    </p><br>
    <p style='width: 470px; text-align:left;margin-bottom:1px;margin-top:1px'>
    <font size='3' face='Georgia' color='330033'><i>Please enter a valid email
    address so that you can be notified upon job completion.</i></font><br></p>
    <center><input type='text' name='email' size="40" maxlength="60"></input>
    <br><br><input id='submit' TYPE='submit' value=' Submit '></form></center>
    
    <script type="text/javascript">
    \$('#submit').submit(function(){
    \$('input[type=submit]', this).attr('disabled', 'disabled');
    }); 
    </script> 
    
    </body></html>
ENDHTML
}

#After the original input is given...#

elsif ( !$query->param('button') && !$query->param('dropdown') )
{
    my $SWISSPROT = '../databases/Swissprot/uniprot_sprot.fasta';
    my $UNIPROT   = '../databases/UniProt/UniProt.fasta';
    my $NT        = '../databases/nt/nt.fasta';
    my $list      = 0;
    my $list2     = 0;
    my $number    = 0;

    my (
        $accession, $alns_fh,   $dbfetch,    $des,       $hit,
        $hit_check, $input_hit, $match_line, $org,       $organism,
        $tmp_fh,    @input_org, @organism,   %organisms, %seq
    );
    print '<center>';
    my $string = $query->param('sequence');
    if ( $string =~ /^\s*$/ )
    {
        print <<"ENDHTML";
        <br><br<br><br><br>
        <font size='2' face='Georgia' color='330033'>
        <b>ERROR!</b> No sequence entered!<br><br>
        Please <a href='javascript:history.go(-1)'><u>go back</u></a> and enter
        a sequence.
        </font>
ENDHTML
        exit;
    }
    else
    {
        if ( $string =~ /(>)/ )
        {
            if ( $' =~ /(>)/ )
            {
                $string = '>' . $`;
            }
        }
        elsif ( $string !~ /^>/ )
        {
            $string = ">\n" . $string;
        }
        if (
            ( $string =~ /OS\=(\w+)\s+(\w+)/ )    #For proteins#
            || ( $string =~ /\[(\w+)\s+(\w+)/ )
          )
        {
            $organism = $1 . q{ } . $2;
            if (   ( $' =~ /^$/ )
                || ( $' =~ /^\s+$/ ) )
            {
                print <<"ENDHTML";
                <br><br<br><br><br><br>
                <font size='2' face='Georgia' color='330033'>
                <b>ERROR!</b> No sequence entered!<br><br>
                Please <a href='javascript:history.go(-1)'><u>go back</u></a>
                and enter a sequence.
                </font>
ENDHTML
                exit;
            }
            if ( $string =~ /[|](\w{6,})[|]/ )
            {
                $hit_check = $1;
                if ( $hit_check =~ /\d/ && $hit_check =~ /\D/ )
                {
                    $input_hit = $hit_check;
                }
            }
        }
        elsif ( $string =~ /ref[|](.+)[|]/ )    #For DNA from Refseq
        {
            $input_hit = $1;
            $dbfetch   = get(
                "http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=refseqn&id=$1");
            if ( $dbfetch =~ /ORGANISM\s+(\w+\s+\w+)/ )
            {
                $organism = $1;
            }
        }
        elsif (
            $string    =~ /gb[|](.+)[|]/        #For DNA from Genbank/EMBL/DBJ
            || $string =~ /dbj[|](.+)[|]/
            || $string =~ /emb[|](.+)[|]/
          )
        {
            $input_hit = $1;
            $dbfetch =
              get("http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=embl&id=$1");
            if ( $dbfetch =~ /OS\s+(\w+\s+\w+)/ )
            {
                $organism = $1;
            }
        }
    }

    # Database selection #
    my $db = $query->param('db');
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

    my $email = $query->param('email');
    if (   ( $email =~ /^$/ )
        || ( $email =~ /^\s+$/ ) )
    {
        print <<"ENDHTML";
        <br><br<br><br><br><br>
        <font size='2' face='Georgia' color='330033'>
        <b>ERROR!</b> No email address entered!<br><br>
        Please <a href='javascript:history.go(-1)'><u>go back</u></a> and enter 
        a valid email address.
        </font>
ENDHTML
        exit;
    }
    elsif ( !is_email($email) )
    {
        print <<"ENDHTML";
        <br><br<br><br><br><br>
        <font size='2' face='Georgia' color='330033'>
        <b>ERROR!</b> The email address does not appear to be valid!<br><br>
        Please <a href='javascript:history.go(-1)'><u>go back</u></a> and enter 
        a valid email address.
        </font>
ENDHTML
        exit;
    }

    print <<"ENDHTML";
    <div id="loading" class='unhidden'>
    <br><br<br><br><br><br>
    <center><img src="../loading.gif"></center><br>
    </div>
ENDHTML

    #Get process id and set BLAST parameters#
    my $prid = $$;
    my $tmp  = '../tmps/blast/' . $prid . '.tmp';
    my $out  = '../outs/blast/' . $prid . '.tmp';
    open $tmp_fh, '>', $tmp or die $!;
    print {$tmp_fh} "$string";
    close $tmp_fh or die $!;

    print "<!--\n";
    my $timeout = 60;    #Avoiding browser timeout.
    $SIG{ALRM} = sub { print ".\n"; alarm $timeout; };
    alarm $timeout;

    my $cpu_n = Sys::CPU::cpu_count();    # get the number of cpu cores

    if ( $db eq $NT )
    {

        #BLASTN for DNA.#
        system(
"blastn -query $tmp -db $db -evalue 0.00000001 -num_threads $cpu_n -out $out"
          ) == 0
          or die $?;
    }
    else
    {

        #BLASTP for Proteins.#
        system(
"blastp -query $tmp -db $db -evalue 0.00000001 -num_threads $cpu_n -out $out"
          ) == 0
          or die $?;
    }

    my $blast = Bio::SearchIO->new(
        -format => 'blast',
        -file   => "$out"
    );

    my $alns      = '../results/final_alns/' . $prid . '.tmp';
    my $alignment = Bio::AlignIO->new(
        -file   => ">>$alns",
        -format => 'clustalw'
    );

    my $hit_old = 0;

    #Start parsing BLAST output.#
    while ( my $result = $blast->next_result )
    {
        while ( my $hit = $result->next_hit )
        {
            while ( my $hsp = $hit->next_hsp )
            {
                if ( $hit->description =~ /(OS\=\w+)\s+(\w+)/ )
                {
                    $des = $1 . q{ } . $2 . $';
                    if ( $des =~ /OS\=(\w+)\s+(\w+)/ )
                    {
                        $org = $1 . q{ } . $2;

                        # Populate the organism dropdown list.#
                        $organism[$list] = $org;

                        $list++;

                        #Populate an array with results from defined organism.#
                        if (   defined $organism
                            && $org eq $organism
                            && defined $input_hit )
                        {
                            if ( $hit->accession =~ /tr[|](\w+)[|]/ )
                            {
                                $input_org[$list2] = $1;
                            }
                            else
                            {
                                $input_org[$list2] = $hit->accession;
                            }
                            $list2++;
                        }

                    #Populate a hash with the first result from every organism.#
                        if ( !( defined $organisms{$org} ) )
                        {
                            if ( $hit->accession =~ /tr[|](\w+)[|]/ )
                            {
                                $organisms{$org} = $1;
                            }
                            else
                            {
                                $organisms{$org} = $hit->accession;
                            }
                        }
                    }
                }
                else
                {
                    if ( $hit->name =~ /ref[|](.+)[|]/ )
                    {
                        if ( $1 =~ /[.]\d*/ )
                        {
                            $accession = $`;
                        }
                        $dbfetch = get(
"http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=refseqn&id=$accession"
                        );

                        # Populate the organism dropdown list.#
                        if ( $dbfetch =~ /ORGANISM\s+(\w+\s+\w+)/ )
                        {
                            $org = $1;
                            $organism[$list] = $org;
                            $list++;
                        }
                    }
                    elsif ($hit->name =~ /gb[|](.+)[|]/
                        || $hit->name =~ /dbj[|](.+)[|]/
                        || $hit->name =~ /emb[|](.+)[|]/ )
                    {
                        if ( $1 =~ /[.]\d*/ )
                        {
                            $accession = $`;
                        }
                        $dbfetch = get(
"http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=embl&id=$accession"
                        );

                        # Populate the organism dropdown list.#
                        if ( $dbfetch =~ /OS\s+(\w+\s+\w+)/ )
                        {
                            $org = $1;
                            $organism[$list] = $org;
                            $list++;
                        }
                    }

                    #Populate an array with results from defined organism.#
                    if ( $org eq $organism && defined $input_hit )
                    {
                        $input_org[$list2] = $accession;
                        $list2++;
                    }

                    #Populate a hash with the first result from every organism.#
                    if ( !( defined $organisms{$org} ) )
                    {
                        $organisms{$org} = $accession;
                    }
                    if ( $hit_old ne $accession )
                    {
                        my $aln = $hsp->get_aln;
                        $alignment->write_aln($aln);
                        open $alns_fh, '>>', $alns or die $!;
                        print {$alns_fh} " \n";
                        close $alns_fh or die $!;
                        $match_line = $hsp->hit_string();
                        $match_line =~ tr/- //d;
                        $seq{"$number"} =
                          ">$accession $org" . q{ } . $match_line;
                        $number++;
                    }
                    $hit_old = $accession;
                }
            }
        }
    }
    if ( defined $input_hit )
    {
        foreach my $hit (@input_org)
        {
            if ( $hit eq $input_hit )
            {
                $organisms{$organism} = $input_hit;
            }
        }
    }
    my @organism2 = uniq(@organism);      #Show each organism only once.
    my $mikos     = @organism2;
    my $orghash   = freeze %organisms;    #Pass the hash to the next CGI call.
    alarm 0;
    print "-->\n";
    if ( defined $organism )
    {
        print <<"ENDHTML";
        <font size='2' face='Georgia' color='330033'>
        <br><br<br><br><br><br>
        It seems that the source organism is <b>$organism</b>. Is this correct?
        </font>
ENDHTML
    }
    print <<"ENDHTML";
    <font size='2' face='Georgia' color='330033'>
    <br><br><br><br>
    Please select the correct source organism from the following list.
    </font><br><br>
    <form action='$0' method='POST'>
    <p style='width: 270px; text-align:center;margin-bottom:1px;margin-top:1px'>
    <div-align='center'><select name='organism' tabindex='1' class='default'>
ENDHTML

    #Set the dropdown list.#
    for ( 0 .. $mikos - 1 )
    {
        print "<option value='$organism2[$_]'>$organism2[$_]</option>\n";
    }
    print <<"ENDHTML";
    <input type=hidden name='prid' value='$prid'>
    <input type=hidden name='db' value='$db'>
    <input type=hidden name='email' value='$email'>
    <input type=hidden name='organisms' value='$orghash'>
ENDHTML
    if ( $db eq $NT )
    {
        my $seqhash = freeze %seq;
        print <<"ENDHTML";
    <input type=hidden name='seq' value='$seqhash'>
    <input type=hidden name='number' value='$number'>
ENDHTML
    }
    print <<"ENDHTML";
    </select></div><input type=submit name='dropdown' value='OK'>
    </form>
    </p>
        
    <script type="text/javascript">
    document.getElementById("loading").className = "hidden";
    </script>
    
    </body>
    </html>
ENDHTML
}

#On to the serious computation phase...#
elsif ($query->param('organism')
    && $query->param('prid')
    && $query->param('db')
    && $query->param('email') )
{
    my $seq2counter      = 0;
    my $hit_old          = 0;
    my $iteration_number = 0;
    my $resnum           = 0;
    my $resnum2          = 0;
    my $sequences_number = 0;
    my $tdcounter        = 0;

    my (
        $acnumber, $alignment, $alns_fh,    $dbfetch,    $des,
        $hit,      $i,         $line,       $match_line, $number,
        $org1,     $out,       $out_fh,     $sequence,   $sequences_fh,
        $tdbg,     @accession, @sequences2, @reslines,   @seq,
        @seq2,     %seq
    );
    my $start_timer = time;

    my $organism = $query->param('organism');
    my $prid     = $query->param('prid');
    my $db       = $query->param('db');
    my $email    = $query->param('email');

    if ( defined $query->param('seq') && defined $query->param('number') )
    {
        %seq    = thaw $query->param('seq');
        $number = $query->param('number');
        for ( 0 .. $number - 1 )
        {
            my $temp = $seq{$_} . "\n\n";
            if ( $temp =~ /(\s\w+\s\w+\s)/ )
            {
                $temp = $` . $1 . "\n" . $';
            }
            $seq[$_] = $temp;
        }
    }

    #Retrieve the hash from before.#
    my %organisms = thaw $query->param('organisms');

    my $one = $organisms{$organism};
    if ( $one =~ /(\w+)[.]\d*/ )
    {
        $one = $1;
    }
    my $tmp = '../tmps/blast/' . $prid . '.tmp';

    print "<!--\n";
    my $timeout = 60;    #Avoiding browser timeout.
    $SIG{ALRM} = sub { print ".\n"; alarm $timeout; };
    alarm $timeout;

    my $cpu_n = Sys::CPU::cpu_count();    #Get the number of cpu cores.
    my $e_th  = '0.00000001';             #E-value
    if ( $db !~ /nt[.]fasta/ )
    {
        $out = '../outs/psiblast/' . $prid . '.tmp';
        my $psib = 'blastpgp';            #Blast Algorithm.

    #Run the BLAST+ executable through the compatibility layer of the old BLAST
    #programs, so as to keep the maximum iterations option that, at first sight,
    #is missing from the BLAST+ ones.
        system(
"legacy_blast.pl $psib -i $tmp -b 7000 -j 50 -h $e_th -d $db -a $cpu_n -o $out"
          ) == 0
          or die $?;

        #Shorten the PSI-BLAST output file. We only need the last iteration.#
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

        my $alns = '../results/final_alns/' . $prid . '.tmp';

        #Start parsing PSI-BLAST output.#
        while ( my $result = $blast->next_result )
        {
            while ( my $it = $result->next_iteration )
            {
                $i = 1000.0;
                my $number = 0;
                while (( $hit = $it->next_hit_new )
                    || ( $hit = $it->next_hit_old ) )
                {
                    while ( my $hsp = $hit->next_hsp )
                    {
                        if ( ( $it->number ) > $iteration_number )
                        {
                            $iteration_number = $it->number;
                        }
                        if ( $hit->description =~ /(OS\=\w+)\s+(\w+)/ )
                        {
                            $des = $1 . q{ } . $2 . $';
                            if ( $des =~ /OS\=(\w+)\s+(\w+)/ )
                            {
                                $org1 = $1 . q{ } . $2;
                            }
                            if ( $org1 eq $organism )
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
                                if ( $hit_old ne $acnumber )
                                {
                                    $dbfetch = get(
"http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=uniprotkb&id=$acnumber&format=fasta&style=raw"
                                    );
                                    if ( $dbfetch =~ /\n/ )
                                    {
                                        $match_line = $';
                                    }
                                    $seq[$number] =
                                      ">$accession[$number] $org1\n"
                                      . $match_line . "\n";
                                    $number++;
                                }
                                $hit_old = $acnumber;
                            }
                        }
                    }
                }
            }
        }
    }

    #Append sequences to file.#
    @seq2 = uniq(@seq);
    my $sequences = '../seq/final_seq/' . $prid . '.tmp';
    open $sequences_fh, '>', $sequences or die $!;
    foreach my $sequence (@seq2)
    {
        if ( $sequence =~ /$organism/ )
        {
            if ( $sequence =~ /(\w+)[.]d*/ )
            {
                $sequence = $` . $1 . $';
            }
            print {$sequences_fh} $sequence;
        }
    }
    close $sequences_fh or die $!;
    open $sequences_fh, '<', $sequences or die $!;
    my $line2;
    {
        local $/ = "\n";
        while ( $line2 = <$sequences_fh> )
        {
            $line2 =~ s/$one/***$one***/;    #Add stars (***) to the input seq.
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

    my $fnaln     = '/web/results/final_alns/multalign/' . $prid . '.aln';
    my $ph        = '/web/results/trees/phs/' . $prid . '.ph';
    my $phb       = '/web/results/trees/phbs/' . $prid . '.phb';
    my $drawntree = '/web/results/trees/drawn/' . $prid . '.png';
    my $zip       = '/web/results/trees/zips/' . $prid . '.zip';

    open $sequences_fh, '<', $sequences or die $!;
    my $line3;
    {
        local $/ = "\n\n";
        while ( $line3 = <$sequences_fh> )
        {
            $sequences_number++;
        }
    }
    close $sequences_fh or die $!;

    my $email2 = $email;
    $email2 =~ s/([\@])/\\$1/;
    if ( $sequences_number > 3 )
    {

        #Alignment, NJ tree plotting, bootstrapping and parsing.#
        tcoffee( $sequences, $fnaln, $email2 );
        system("clustalw -INFILE=$fnaln -OUTFILE=$ph -tree") == 0 or die $?;
        system(
"clustalw -INFILE=$fnaln -OUTFILE=$phb -bootstrap=1000 -bootlabels=node"
          ) == 0
          or die $?;
        system("rm *.dnd") == 0 or die $?;
        rename "../results/final_alns/multalign/$prid.ph",
          "../results/trees/phs/$prid.ph"
          or die $!;
        rename "../results/final_alns/multalign/$prid.phb",
          "../results/trees/phbs/$prid.phb"
          or die $!;
        tree_manipulation1($phb);
        system("./Pinda.R -parser $phb > /web/parsing/$prid.tmp") == 0
          or die $?;    #Parse.
        system("./Pinda.R -lengths_1 $phb > /web/parsing/$prid\_1.tmp") == 0
          or die $?;    #Parse.
        system("./Pinda.R -lengths_2 $phb > /web/parsing/$prid\_2.tmp") == 0
          or die $?;    #Parse again.
        tree_manipulation2($phb);
        system("zip -j $zip $ph $phb") == 0 or die $?;
        system("./Pinda.R $drawntree $phb") == 0
          or die $?;    #Visually draw the tree.
        system("rm $ph $phb") == 0 or die $?;

        alarm 0;

        parser( "../parsing/$prid.tmp", "../parsing/$prid\_1.tmp",
            "../parsing/$prid\_2.tmp" );
        print <<"ENDHTML";
        -->
        <center><br><font size='3' face='Georgia' color='330033'>
        <a href=../results/final_alns/multalign/$prid.aln>T-Coffee Alignment</a>
        </font>
        <br><br>
        <table border='1'>
        <tr bgcolor=FFFF66><th><center>Possible duplications of
ENDHTML
        $email_data .= <<"EMAIL_END";
        <center><br><font size='3' face='Georgia' color='330033'>
        <a href=http://localhost/results/final_alns/multalign/$prid.aln>T-Coffee
         Alignment</a>
        </font>
        <br><br>
        <table border='1'>
        <tr bgcolor=FFFF66><th><center>Possible duplications of
EMAIL_END
        if ( $db !~ /nt[.]fasta/ )
        {
            print <<"ENDHTML";
            <a href=http://www.uniprot.org/uniprot/$one>
            $one</a>.</center></th>
            <th><center>Confidence Value</center></th></tr>
ENDHTML
            $email_data .= <<"EMAIL_END";
            <a href=http://www.uniprot.org/uniprot/$one>
            $one</a>.</center></th>
            <th><center>Confidence Value</center></th></tr>
EMAIL_END
        }
        else
        {
            if ( $one =~ /^\D{2}\_/ )
            {
                print <<"ENDHTML";
                <a href=http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=refseqn&id=$one&style=raw>
                $one</a>.</center></th>
                <th><center>Confidence Value</center></th></tr>
ENDHTML
                $email_data .= <<"EMAIL_END";
                <a href=http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=refseqn&id=$one&style=raw>
                $one</a>.</center></th>
                <th><center>Confidence Value</center></th></tr>
EMAIL_END
            }
            else
            {
                print <<"ENDHTML";
                <a href=http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=embl&id=$one&style=raw>
                $one</a>.</center></th>
                <th><center>Confidence Value</center></th></tr>
ENDHTML
                $email_data .= <<"EMAIL_END";
                <a href=http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=embl&id=$one&style=raw>
                $one</a>.</center></th>
                <th><center>Confidence Value</center></th></tr>
EMAIL_END
            }
        }

        #Table of Confidence Value.#
        my $top_can;
        if ( $candidate[0] =~ /(\d?\d?\d?.?\d+e?-?\d*) \w+/ )   #Input sequence.
        {
            $top_can = $1;
        }

        foreach my $can (@candidate)
        {
            if (
                $can !~ /$starting_point/    #Not the input sequence.
                && $can =~ /(\d?\d?\d?.?\d+e?-?\d*) (\w+)/
              )                              #All the rest.
            {
                my $conf_temp = $1;
                my $uni_temp  = $2;
                if ( $tdcounter == 1 )       #Color the html table alternately.
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
                    print
"<tr bgcolor=$tdbg><td><center><a href=http://www.uniprot.org/uniprot/$2>$2</a>
</center></td>";
                    $email_data .=
"<tr bgcolor=$tdbg><td><center><a href=http://www.uniprot.org/uniprot/$2>$2</a>
</center></td>";
                }
                else
                {
                    if ( $uni_temp =~ /^\D{2}\_/ )
                    {
                        print
"<tr bgcolor=$tdbg><td><center><a href=http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=refseqn&id=$uni_temp&style=raw>$uni_temp</a>
</center></td>";
                        $email_data .=
"<tr bgcolor=$tdbg><td><center><a href=http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=refseqn&id=$uni_temp&style=raw>$uni_temp</a>
</center></td>";
                    }
                    else
                    {
                        print
"<tr bgcolor=$tdbg><td><center><a href=http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=embl&id=$uni_temp&style=raw>$uni_temp</a>
</center></td>";
                        $email_data .=
"<tr bgcolor=$tdbg><td><center><a href=http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=embl&id=$uni_temp&style=raw>$uni_temp</a>
</center></td>";
                    }
                }
                printf '<td align=left><center>%5.1f%%</center></td></tr>',
                  $conf_temp / $top_can * 100;

                $email_data .=
                  sprintf '<td align=left><center>%5.1f%%</center></td></tr>',
                  $conf_temp / $top_can * 100;
            }
        }
        print '</table>';
        print <<"ENDHTML";
        <img src='../results/trees/drawn/$prid.png'><br>
        <a href=../results/trees/zips/$prid.zip>NJ Tree Produced</a>
ENDHTML
        $email_data .= <<"EMAIL_END";
        </table><img src='http://localhost/results/trees/drawn/$prid.png'><br>
        <a href=http://localhost/results/trees/zips/$prid.zip>NJ Tree Produced
        </a>
EMAIL_END

    }
    else
    {
        if ( $sequences_number == 3 )
        {
            alarm 0;
            print <<"ENDHTML";
            -->
            <font size='3' face='Georgia' color='330033'><br><br>
            <center>
            The number of sequences for this particular organism is three.
            <br><b>The dendrogram cannot be bootstrapped.</b></font>
            </center>
ENDHTML
            $email_data .= <<"EMAIL_END";
<font size='3' face='Georgia' color='330033'><br><br>
            <center>
            The number of sequences for this particular organism is three.
            <br><b>The dendrogram cannot be bootstrapped.</b></font>
            </center>
EMAIL_END

            #Alignment, NJ tree plotting and parsing.#
            tcoffee( $sequences, $fnaln, $email2 );
            system("clustalw -INFILE=$fnaln -OUTFILE=$ph -tree") == 0 or die $?;
            rename "../results/final_alns/multalign/$prid.ph",
              "../results/trees/phs/$prid.ph"
              or die $!;
            system("rm *.dnd") == 0 or die $?;
            tree_manipulation1($ph);
            system("./Pinda.R -parser $ph > /web/parsing/$prid.tmp") == 0
              or die $?;    #Parse.
            system("./Pinda.R -lengths_1 $ph > /web/parsing/$prid\_1.tmp") == 0
              or die $?;    #Parse.
            system("./Pinda.R -lengths_2 $ph > /web/parsing/$prid\_2.tmp") == 0
              or die $?;    #Parse again.

            tree_manipulation2($ph);

            system("zip -j $zip $ph") == 0 or die $?;
            system("./Pinda.R $drawntree $ph") == 0
              or die $?;    #Visually draw the tree.
            system("rm $ph") == 0 or die $?;

            parser( "../parsing/$prid.tmp", "../parsing/$prid\_1.tmp",
                "../parsing/$prid\_2.tmp" );
            print <<"ENDHTML";
            <center><br><font size='3' face='Georgia' color='330033'>
            <a href=../results/final_alns/multalign/$prid.aln>T-Coffee Alignment
            </a>
            </font>
            <br><br>
            <table border='1'>
            <tr bgcolor=FFFF66><th><center>Possible duplications of
ENDHTML
            $email_data .= <<"EMAIL_END";
            <center><br><font size='3' face='Georgia' color='330033'>
            <a href=http://localhost/results/final_alns/multalign/$prid.aln>T-Coffee
            Alignment</a>
            </font>
            <br><br>
            <table border='1'>
            <tr bgcolor=FFFF66><th><center>Possible duplications of
EMAIL_END
            if ( $db !~ /nt[.]fasta/ )
            {
                print <<"ENDHTML";
                <a href=http://www.uniprot.org/uniprot/$one>
                $one</a>.</center></th>
                <th><center>Confidence Value</center></th></tr>
ENDHTML
                $email_data .= <<"EMAIL_END";
                <a href=http://www.uniprot.org/uniprot/$one>
                $one</a>.</center></th>
                <th><center>Confidence Value</center></th></tr>
EMAIL_END
            }
            else
            {
                if ( $one =~ /^\D{2}\_/ )
                {
                    print <<"ENDHTML";
                    <a href=http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=refseqn&id=$one&style=raw>
                    $one</a>.</center></th>
                    <th><center>Confidence Value</center></th></tr>
ENDHTML
                    $email_data .= <<"EMAIL_END";
                    <a href=http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=refseqn&id=$one&style=raw>
                    $one</a>.</center></th>
                    <th><center>Confidence Value</center></th></tr>
EMAIL_END
                }
                else
                {
                    print <<"ENDHTML";
                    <a href=http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=embl&id=$one&style=raw>
                    $one</a>.</center></th>
                    <th><center>Confidence Value</center></th></tr>
ENDHTML
                    $email_data .= <<"EMAIL_END";
                    <a href=http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=embl&id=$one&style=raw>
                    $one</a>.</center></th>
                    <th><center>Confidence Value</center></th></tr>
EMAIL_END
                }
            }

            my $top_can;
            if ( $candidate[0] =~ /(\d?\d?\d?.?\d+e?-?\d*) \w+/ )    #Input seq.
            {
                $top_can = $1;
            }

            foreach my $can (@candidate)
            {
                if (
                    $can !~ /$starting_point/    #Not the input sequence.
                    && $can =~ /(\d?\d?\d?.?\d+e?-?\d*) (\w+)/
                  )                              #All the rest.
                {
                    if ( $tdcounter == 1 )    #Color the html table alternately.
                    {
                        $tdbg      = 'F8FBFE';
                        $tdcounter = 0;
                    }
                    else
                    {
                        $tdbg = 'EAF1FB';
                        $tdcounter++;
                    }
                    print
"<tr bgcolor=$tdbg><td><center><a href=http://www.uniprot.org/uniprot/$2>$2</a>
</center></td>";
                    printf '<td align=left><center>%5.1f%%</center></td></tr>',
                      $1 / $top_can * 100;
                    $email_data .=
"<tr bgcolor=$tdbg><td><center><a href=http://www.uniprot.org/uniprot/$2>$2</a>
</center></td>"
                      . sprintf
                      '<td align=left><center>%5.1f%%</center></td></tr>',
                      $1 / $top_can * 100;
                }
            }
            print '</table>';
            print <<"ENDHTML";
            <center><img src='../results/trees/drawn/$prid.png'><br>
            <a href=../results/trees/zips/$prid.zip>NJ Tree Produced</a>
            </center>
ENDHTML
            $email_data .= <<"EMAIL_END";
            </table><img src='http://localhost/results/trees/drawn/$prid.png'>
            <br>
            <a href=http://localhost/results/trees/zips/$prid.zip>NJ Tree
            Produced</a>
EMAIL_END

        }
        if ( $sequences_number == 2 )
        {
            alarm 0;
            print <<"ENDHTML";
            -->
            <center>
            <font size='3' face='Georgia' color='330033'><br><br>
            The number of sequences for this particular organism is less than
            three.<br><br>
            <b>A dendrogram cannot be plotted.</b></font>
            </center>
ENDHTML
            $email_data .= <<"EMAIL_END";
            <center>
            <font size='3' face='Georgia' color='330033'><br><br>
            The number of sequences for this particular organism is less than
            three.<br><br>
            <b>A dendrogram cannot be plotted.</b></font>
            </center>
EMAIL_END

            #Alignment only.#
            tcoffee( $sequences, $fnaln, $email );
            print <<"ENDHTML";
                <center><br><font size='3' face='Georgia' color='330033'><br>
                <a href=../results/final_alns/multalign/$prid.aln>
                T-Coffee Alignment</a></font></center>
ENDHTML
            $email_data .= <<"EMAIL_END";
            <center><br><font size='3' face='Georgia' color='330033'><br>
                <a href=../results/final_alns/multalign/$prid.aln>
                T-Coffee Alignment</a></font></center>
EMAIL_END
        }
        if ( $sequences_number == 1 )
        {
            alarm 0;

            #Nothing to do here.#
            print <<"ENDHTML";
            -->
            <center>
            <br><font size='3' face='Georgia' color='330033'><br><br>
            No possible duplications have been detected.</font>
            </center>
ENDHTML
            $email_data .= <<"EMAIL_END";
            <center>
            <br><font size='3' face='Georgia' color='330033'><br><br>
            No possible duplications have been detected.</font>
            </center>
EMAIL_END
        }
    }
    my $end_timer = time;
    my $run_time  = $end_timer - $start_timer;
    job_timer($run_time);    #Calculate how much time that took.

    send_email($email);      #Send email.
    print "</center>\n</body>\n</html>";
}

sub tcoffee                  #Pretty self-explanatory.
{
    system(
        "export HOME_4_TCOFFEE='/web/t-coffee/' ;
    export DIR_4_TCOFFEE='/web/t-coffee/dir_4_t-coffee/' ;
    export CACHE_4_TCOFFEE='/web/t-coffee/dir_4_t-coffee/cache/' ;
    export TMP='/web/t-coffee/dir_4_t-coffee/tmp/' ;
    export NO_ERROR_REPORT_4_TCOFFEE='1' ;
    export NO_WARNING_4_TCOFFEE='1' ;
t_coffee -infile $_[0] -output=clustal -outfile $_[1] -proxy -email=$_[2]"
      ) == 0
      or die $?;
    return 0;
}

#Removing the TRICHOTOMY word, setting negative values to zero.#
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
            $line =~ s/TRICHOTOMY//;    #Remove the 'TRICHOTOMY' word.
            $line =~ s/-\d[.]\d*/0/;

      #Fix NJ negative branch length artifact, by setting those numbers to zero.
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

#Dividing the bootstrap values by 10, to reach a maximum value of 100.#
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

sub parser    #Parsing.#
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

        #Get all the sequences.#
        while ( $line = <$tree_for_parsingfh> )
        {
            if ( $line =~ /\s"(\w\_\_\_\w+\_\_\_)"/ )    #Input sequence.
            {
                $starting_point = $1;
                $uni_ids[$unicounter] = $1;
                $unicounter++;
            }
            if ( $line =~ /"(\w{6,})"/ )                 #Other sequences.
            {
                if ( $1 !~ /\D\d{1,4}\_/ )
                {
                    $uni_ids[$unicounter] = $1;
                    $unicounter++;
                }
            }
            if ( $' =~ /"(\w{6,})"/ && $1 !~ /\D\d{1,4}\_/ )   #Other sequences.
            {
                $uni_ids[$unicounter] = $1;
                $unicounter++;
            }
            $parsing_lines[$linocounter] = $line;
            $linocounter++;
        }
    }

    open $tree_node_distancesfh, '<', $_[1]
      or die $!;    #Distance parser between nodes.
    my $neoline;
    {
        local $/ = "\n";
        while ( $neoline = <$tree_node_distancesfh> )
        {
            if ( $neoline !~ /\$/ && $neoline =~ /\w/ )
            {
                do
                {
                    if ( $neoline !~ /[.]/ ) #This line contains sequence names.
                    {
                        if ( $neoline =~ /\s?(\w+)/ || $neoline =~ /(Root)/ )
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
                        $neoline =~ /(\d+[.]\d+)/ )  #This line contains values.
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
    open $tree_tip_distancesfh, '<', $_[2]
      or die $!;    #Distance parser from tip to node.
    my $neoline2;
    {
        local $/ = "\n";
        while ( $neoline2 = <$tree_tip_distancesfh> )
        {
            if ( $neoline2 !~ /\$/ && $neoline2 =~ /\w/ )
            {
                do
                {
                    if (
                        $neoline2 !~ /[.]/ ) #This line contains sequence names.
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
                    elsif (
                        $neoline2 =~ /(\d+[.]\d+)/ ) #This line contains values.
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
            if ( $line =~ /parts\$(\w)(\w+)/ )   #Get the node of the input seq.
            {
                $star_seq[$sscounter] = $1 . $2;
                $search_id = $1 . $2;
                $sscounter++;
            }
            for ( 0 .. $p_l_s - 1 )
            {
                if ( $parsing_lines[$_] =~ /$search_id"/ )
                {
                    if ( $parsing_lines[$_] =~
                        /parts\$(\w)(\w+)/ )     #Get the node's node...
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
                if ( $line =~ /parts\$(\w)(\w+)/ )    #Get the node of a seq.
                {
                    $compare_seq[$sscounter] = $1 . $2;
                    $search_id = $1 . $2;
                    $sscounter++;
                }
                for ( 0 .. $p_l_s - 1 )
                {
                    if ( $parsing_lines[$_] =~ /$search_id"/ )
                    {

                        #Get the node's node...#
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
            if ( $node =~ /(\d+)\_?/ )    #Add the distance for this node.
            {
                $node_distance += $distanceshash{"$node"};
                $ginomenon *= $1 / 1000.0;
            }
            foreach my $star_node (@star_seq)
            {

                #Find the common node between input seq and every other seq.#
                if ( $node eq $star_node )
                {
                    foreach my $star_node (@star_seq)
                    {

                        #Add distances until the common node.#
                        if ( $star_node ne $node && $star_node =~ /(\d+)\_?/ )
                        {
                            $node_distance += $distanceshash{"$star_node"};
                            $ginomenon *= $1 / 1000.0;
                        }

                        #Remove the common node's distance.#
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

                    #For every sequence, except for the input one...
                    if ( $uni !~ /$starting_point/ )
                    {

                        #Add the tip-to-closest-node distances.
                        $node_distance +=
                          $distanceshash2{"$uni"} +
                          $distanceshash2{"$starting_point"};

                        #If distance is 0, then the sequences are overlapping!?!
                        if ( $node_distance == 0 )
                        {
                            $candidate[$cancounter] = $ginomenon . q{ } . $uni;
                            $cancounter++;
                        }

                        #Calculate the confidence value.
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

            #If the input sequence had only the root common with a sequence...#
            foreach my $node (@compare_seq)
            {
                if ( $node =~ /(\d+)\_?/ )
                {

                    #Add the distance from the sequence's node to root.#
                    $node_distance += $distanceshash{"$node"};
                    $ginomenon *= $1 / 1000.0;
                }
            }
            foreach my $star_node (@star_seq)
            {
                if ( $star_node =~ /(\d+)\_?/ )
                {

                    #Add the distance from the input sequence's node to root.#
                    $node_distance += $distanceshash{"$star_node"};
                    $ginomenon *= $1 / 1000.0;
                }
            }

            #Add the tip-to-closest-node distances.#
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
    @candidate = map join( q{ }, @{$_} ), sort { $a->[0] <=> $b->[0] } map {[split]}
      @candidate;
    @candidate = reverse @candidate;    #Reverse the array.#
    return @candidate, $starting_point;
}

sub job_timer    #Calculates the execution time from PSI-BLAST to the end.#
{
    my ( $hours, $minutes, $seconds );
    if ( $_[0] > 3600 )
    {
        $hours = $_[0] / 3600;
        $_[0] %= 3600;
    }
    if ( $_[0] > 60 )
    {
        $minutes = $_[0] / 60;
        $_[0] %= 60;
    }
    if ( $_[0] > 0 )
    {
        $seconds = $_[0];
    }
    print '<br><br><font size="2" face="Courier New"><center>This job took ';
    if ( defined $hours && defined $minutes && defined $seconds )
    {
        if ( $hours >= 2 )
        {
            printf '<b>%.0f hours</b>,', $hours;
        }
        else
        {
            printf '<b>1 hour</b>,';
        }
        if ( $minutes >= 2 )
        {
            printf ' <b>%.0f minutes</b>', $minutes;
        }
        else
        {
            printf ' <b>1 minute</b>';
        }
        if ( $seconds >= 2 )
        {
            printf ' and <b>%.0f seconds</b>.', $seconds;
        }
        else
        {
            printf ' and <b>1 second</b>.';
        }
    }
    elsif ( defined $hours && defined $minutes )
    {
        if ( $hours >= 2 )
        {
            printf '<b>%.0f hours</b>', $hours;
        }
        else
        {
            printf '<b>1 hour</b>';
        }
        if ( $minutes >= 2 )
        {
            printf ' and <b>%.0f minutes</b>.', $minutes;
        }
        else
        {
            printf ' and <b>1 minute</b>.';
        }
    }
    elsif ( defined $hours && defined $seconds )
    {
        if ( $hours >= 2 )
        {
            printf '<b>%.0f hours</b>', $hours;
        }
        else
        {
            printf '<b>1 hour</b>';
        }
        if ( $seconds >= 2 )
        {
            printf ' and <b>%.0f seconds</b>.', $seconds;
        }
        else
        {
            printf ' and <b>1 second</b>.';
        }
    }
    elsif ( defined $minutes && defined $seconds )
    {
        if ( $minutes >= 2 )
        {
            printf '<b>%.0f minutes</b>', $minutes;
        }
        else
        {
            printf '<b>1 minute</b>';
        }
        if ( $seconds >= 2 )
        {
            printf ' and <b>%.0f seconds</b>.', $seconds;
        }
        else
        {
            printf ' and <b>1 second</b>.';
        }
    }
    elsif ( defined $hours )
    {
        if ( $hours >= 2 )
        {
            printf '<b>%.0f hours</b>.', $hours;
        }
        else
        {
            printf '<b>1 hour</b>.';
        }
    }
    elsif ( defined $minutes )
    {
        if ( $minutes >= 2 )
        {
            printf '<b>%.0f minutes</b>.', $minutes;
        }
        else
        {
            printf '<b>1 minute</b>.';
        }
    }
    elsif ( defined $seconds )
    {
        if ( $seconds >= 2 )
        {
            printf '<b>%.0f seconds</b>.', $seconds;
        }
        else
        {
            printf '<b>1 second</b>.';
        }
    }
    print '</font></center>';
    return 0;
}

sub send_email
{
    my $msg = MIME::Lite->new(
        Subject => 'Pinda Job Result',
        From    => 'pinda@mbg.duth.gr',
        To      => $_[0],
        Type    => 'text/html',
        Data    => $email_data
    );
    $msg->send();
    return 0;
}
