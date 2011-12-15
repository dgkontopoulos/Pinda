#!/usr/bin/perl -w

###################
#Pinda source code#
###################

use autodie;
use Bio::Perl;
use Bio::Tools::Run::StandAloneBlastPlus;
use Bio::Tools::Run::StandAloneBlastPlus::BlastMethods;
use Bio::SeqIO;
use Bio::SearchIO;
use Bio::Search::Result::BlastResult;
use Bio::Root::IO;
use Bio::Search::Iteration::GenericIteration;
use Bio::AlignIO;
use CGI qw(:standard);
use Data::Validate::Email qw(is_email);
use List::MoreUtils qw(uniq);
use Sys::CPU;
use Bio::TreeIO;
use Bio::Tree::TreeI;
use Bio::Tree::TreeFunctionsI;

#use strict;
#use warnings;

my $Swissprot        = '../databases/Swissprot/uniprot_sprot.fasta';
my $UniProt          = '../databases/TrEMBL/uniprot.fasta';
my $alncounter       = 0;
my $bscounter        = 1;
my $cancounter       = 0;
my $gcounter         = 0;
my $ginomenon        = 1;
my $hit_old          = 0;
my $it_counting      = 0;
my $linocounter      = 0;
my $sequences_number = 0;
my $iteration_number = 0;
my $it_exit          = 0;
my $jac              = 0;
my $list             = 0;
my $p                = 0;
my $p2               = 0;
my $pathcounter      = 0;
my $parse_counter    = 0;
my $pc_id            = 0;
my $resnum           = 0;
my $resnum2          = 0;
my $sscounter        = 0;
my $unicounter       = 0;
my $yacounter        = 0;
my (
    $acnumber,       $alignment,  $alns,      $alns_fh,
    $bvdiv,          $des,        $drawntree, $fnaln,
    $ge_url,         $gene,       $genepo,    $hit,
    $i,              $line,       $line2,     $match_line,
    $matchreg,       $mikos,      $needed,    $nod,
    $number_of_cpus, $one,        $original,  $org,
    $org1,           $organism,   $out,       $out_fh,
    $ph,             $phb,        $pos,       $prid,
    $results,        $results_fh, $sequence,  $sequences,
    $sequences_fh,   $tmp,        $tmp_fh,    $tree_fh,
    @accession,      @alignments, @nodes,     @nodes2,
    @numero,         @organism,   @organism2, @possible,
    @possible2,      @reslines,   @score,     @seq,
    @seq2,           @evalue,     @tree
);

open STDERR, '>', '/dev/null';

my $query = CGI->new;
print $query->header;
print <<"ENDHTML";
<!DOCTYPE html
	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	 "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" lang="en-US" xml:lang="en-US">
<head>
<title>Pinda - Pipeline for Intraspecies Duplications Analysis</title>
<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1" /> 
<link rel="stylesheet" href="../css/Pinda.css" type="text/css">

ENDHTML
print <<"ENDHTML";
<script src="https://ajax.googleapis.com/ajax/libs/jquery/1.6.2/jquery.min.js" type="text/javascript" charset="utf-8"></script>

<script src="../js/jquery.dropkick.1-0.1.js" type="text/javascript" charset="utf-8"></script>
<script type="text/javascript" charset="utf-8">
\$(function () {
    \$('.default').dropkick();
});
</script>
</head>
<body background='../background.jpg'>
<LINK REL='SHORTCUT ICON' HREF='../pinda.ico'>
<center>
<br>
<a href='$0'><img src='../pindalogo.png'></a>
<br>
<p style='width: 500px; text-align:center;margin-bottom:1px;margin-top:1px'>
<br>
<hr/>
</p>
</center>
ENDHTML

if ( !$query->param )    #Homepage
{
    print <<"ENDHTML";
    <center>
    <p style='width: 470px; text-align:left;margin-bottom:1px;margin-top:1px'>
    <br><font size='3' face='Georgia' color='330033'>
    <i>Please enter a sequence.</i>
    </font><br>
    <form id = 'submit' method = 'post' action=$0>
    <textarea name='sequence' rows=7 cols=55 style='font-family:courier new'></textarea>
    </p><br>
    <p style='width: 470px; text-align:left;margin-bottom:1px;margin-top:1px'>
    <font size='3' face='Georgia' color='330033'><i>Please choose a database.</i></font><br>
    <center>
    <font size='2'>
    <input type="radio" name="db" value="Swiss-Prot" checked> <b>Swiss-Prot</b>
    <input type="radio" name="db" value="UniProt"> <b>UniProt</b> (<i>Swiss-Prot + TrEMBL</i>)
    </font>
    </center>
    </p><br>
    <p style='width: 470px; text-align:left;margin-bottom:1px;margin-top:1px'>
    <font size='3' face='Georgia' color='330033'><i>Please enter a valid email address so that you can be notified upon job completion.</i></font><br>
    </p>
    <center>
    <input type='text' name='email' size="40" maxlength="60"></input>
    <br><br>
    <input id='submit' TYPE='submit' value=' Submit '></form>
    </center>
    <script type="text/javascript">
    \$('#submit').submit(function(){
    \$('input[type=submit]', this).attr('disabled', 'disabled');
    }); 
    </script> 
    </body>
    </html>
ENDHTML
}

#After the original input is given...#

elsif ( !$query->param('button') && !$query->param('dropdown') )
{
    print '<center>';
    my $string = $query->param('sequence');
    if ( $string =~ /^$/ )
    {
        print <<"ENDHTML";
        <br><br<br><br><br>
        <font size='2' face='Georgia' color='330033'>
        <b>ERROR!</b> No sequence entered!<br><br>
        Please <a href='javascript:history.go(-1)'><u>go back</u></a> and enter a sequence.
        </font>
ENDHTML
        exit;
    }
    if ( $string !~ /[DEFHIKLMNPQRSVWY]/i )
    {
        print <<"ENDHTML";
        <br><br<br><br><br><br>
        <font size='2' face='Georgia' color='330033'>
        THAT... IS... DNA!
        </font>
ENDHTML
        exit;
    }
    else
    {
        if ( $string !~ /^>/ )
        {
            $string = ">\n" . $string;
        }
        if (   ( $string =~ /OS\=(\w+)\s+(\w+)/ )
            || ( $string =~ /\[(\w+)\s+(\w+)\]/ ) )
        {
            $organism = $1 . ' ' . $2;
            if (   ( $' =~ /^$/ )
                || ( $' =~ /^\s+$/ ) )
            {
                print <<"ENDHTML";
                <br><br<br><br><br><br>
                <font size='2' face='Georgia' color='330033'>
                <b>ERROR!</b> No sequence entered!<br><br>
                Please <a href='javascript:history.go(-1)'><u>go back</u></a> and enter a sequence.
                </font>
ENDHTML
                exit;
            }
            elsif ( $' !~ /[DEFHIKLMNPQRSVWY]/i )
            {
                print <<"ENDHTML";
                <br><br<br><br><br><br>
                <font size='2' face='Georgia' color='330033'>
                THAT... IS... DNA!
                </font>
ENDHTML
                exit;
            }
        }

    }
    my $db = $query->param('db');
    if ( $db eq 'Swiss-Prot' )
    {
        $db = $Swissprot;
    }
    elsif ( $db eq 'UniProt' )
    {
        $db = $UniProt;
    }

    my $email = $query->param('email');
    if (   ( $email =~ /^$/ )
        || ( $email =~ /^\s+$/ ) )
    {
        print <<"ENDHTML";
        <br><br<br><br><br><br>
        <font size='2' face='Georgia' color='330033'>
        <b>ERROR!</b> No email address entered!<br><br>
        Please <a href='javascript:history.go(-1)'><u>go back</u></a> and enter a valid email address.
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
        Please <a href='javascript:history.go(-1)'><u>go back</u></a> and enter a valid email address.
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
    $prid = $$;
    $tmp  = '../tmps/blastp/' . $prid . '.tmp';
    $out  = '../outs/blastp/' . $prid . '.tmp';
    open $tmp_fh, '>', $tmp;
    print {$tmp_fh} "$string";
    close $tmp_fh;

    my $timeout = 60;
    print "<!--\n";
    $SIG{ALRM} = sub { print ".\n"; alarm $timeout; };
    alarm $timeout;

    $number_of_cpus = Sys::CPU::cpu_count();    # get the number of cpu cores
`blastp -query $tmp -db $db -evalue 0.00000001 -num_threads $number_of_cpus -out $out`;

    my $blast = Bio::SearchIO->new(
        -format => 'blast',
        -file   => "$out"
    );

    #Start parsing BLAST output.#

    while ( my $result = $blast->next_result )
    {
        while ( my $hit = $result->next_hit )
        {
            while ( my $hsp = $hit->next_hsp )
            {
                if ( $hit->description =~ /(OS\=\w+)\s+(\w+)/ )
                {
                    $des = $1 . ' ' . $2 . $';
                    if ( $des =~ /OS\=(\w+)\s+(\w+)/ )
                    {
                        $org = $1 . ' ' . $2;
                        $organism[$list] =
                          $org;    # Populate the organism dropdown list.
                        $list++;
                    }
                }
            }
        }
    }
    @organism2 = uniq(@organism);
    $mikos     = @organism2;
    alarm 0;
    print "-->\n";
    if ( defined $organism )
    {
        print <<"ENDHTML";
        <font size='2' face='Georgia' color='330033'>
        <br><br<br><br><br><br>
        It seems that the organism source is <b>$organism</b>. Is this correct?
        </font><br>
ENDHTML
    }
    print <<"ENDHTML";
    <font size='2' face='Georgia' color='330033'>
    <br><br>
    Please select the correct organism source from the following list.
    </font><br><br>
    <form action='$0' method='POST'>
    <p style='width: 250px; text-align:center;margin-bottom:1px;margin-top:1px'>
    <div-align='center'><select name='organism' tabindex='1' class='default'>
ENDHTML
    for ( 0 .. $mikos - 1 )
    {
        print "<option value='$organism2[$_]'>$organism2[$_]</option>\n";
    }
    print <<"ENDHTML";
    <input type=hidden name='prid' value='$prid'>
    <input type=hidden name='db' value='$db'>
    <input type=hidden name='email' value='$email'>
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

elsif ($query->param('organism')
    && $query->param('prid')
    && $query->param('db')
    && $query->param('email') )
{
    my $organism = $query->param('organism');
    my $prid     = $query->param('prid');
    my $db       = $query->param('db');
    my $email    = $query->param('email');

    $tmp = '../tmps/blastp/' . $prid . '.tmp';
    $out = '../outs/psiblast/' . $prid . '.tmp';
    my $seqfblast = Bio::SeqIO->newFh( -file => $tmp, -format => 'fasta' );
    my $seqfblast2 = <$seqfblast>;

    my $timeout = 60;
    print "<!--\n";
    $SIG{ALRM} = sub { print ".\n"; alarm $timeout; };
    alarm $timeout;

    $number_of_cpus = Sys::CPU::cpu_count();    # get the number of cpu cores
`legacy_blast.pl blastpgp -i $tmp -b 7000 -j 50 -h 0.00000001 -d $db -a $number_of_cpus -o $out`;

    open $out_fh, '<', $out;
    while ( $line = <$out_fh> )
    {
        if ( $line =~ /Results from round (\d+)/ )
        {
            $resnum = $1;
        }
    }
    close $out_fh;
    open $out_fh, '<', $out;
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
    close $out_fh;
    open $out_fh, '>', $out;
    foreach my $line (@reslines)
    {
        print {$out_fh} $line;
    }
    close $out_fh;

    my $blast = Bio::SearchIO->new(
        -format => 'blast',
        -file   => "$out"
    );

    my $alns      = '../results/final_alns/' . $prid . '.tmp';
    my $results   = '../results/final_results/' . $prid . '.tmp';
    my $sequences = '../seq/final_seq/' . $prid . '.tmp';

    $alignment = Bio::AlignIO->new(
        -file   => ">>$alns",
        -format => 'clustalw'
    );

    #Start parsing BLAST output.#

    while ( my $result = $blast->next_result )
    {
        while ( my $it = $result->next_iteration )
        {
            $i = 1000.0;
            my $number = 1;
            while (( $hit = $it->next_hit_new )
                || ( $hit = $it->next_hit_old ) )
            {
                while ( my $hsp = $hit->next_hsp )
                {
                    if ( ( $it->number ) > $iteration_number )
                    {
                        open $results_fh, '>', $results;
                        print {$results_fh} "\n";
                        open $alns_fh, '>', $alns;
                        print {$alns_fh} " \n";
                        close $alns_fh;
                        $iteration_number = $it->number;
                    }
                    if ( $hit->description =~ /(OS\=\w+)\s+(\w+)/ )
                    {
                        $des = $1 . ' ' . $2 . $';
                        if ( $des =~ /OS\=(\w+)\s+(\w+)/ )
                        {
                            $org1 = $1 . ' ' . $2;
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
                                if ( $hsp->percent_identity > $pc_id )
                                {
                                    $one   = $acnumber;
                                    $pc_id = $hsp->percent_identity;
                                }
                                $score[$number]  = $hit->score;
                                $evalue[$number] = $hit->significance;
                                $numero[$number] = $it->number;
                                my $aln = $hsp->get_aln;
                                $alignment->write_aln($aln);
                                open $alns_fh, '>>', $alns;
                                print {$alns_fh} "\n";
                                close $alns_fh;
                                $match_line = $hsp->hit_string() . "\n";
                                $match_line =~ tr/- //d;
                                $seq[$number] = ">$accession[$number] $org1\n"
                                  . $match_line . "\n";
                                print {$results_fh}
"<tr><td><p style='text-align:right;margin-bottom:1px;margin-top:1px'><center>$accession[$number]</center></td> <td><center>$score[$number]</center></td> <td><center>$evalue[$number]</center></td></tr>\n\n";
                                $number++;
                            }
                            $hit_old = $acnumber;
                        }
                    }
                }
            }
        }
    }
    close $results_fh;
    alarm 0;
    print "-->\n";

    #Append sequences to file.#
    @seq2 = uniq(@seq);
    open $sequences_fh, '>', $sequences;
    foreach my $sequence (@seq2)
    {
        if ( $sequence =~ /$organism/ )
        {
            print {$sequences_fh} $sequence;
        }
    }
    close $sequences_fh;

    open $sequences_fh, '<', $sequences;
    $/ = "\n";
    while ( $line = <$sequences_fh> )
    {
        $line =~ s/$one/***$one***/;
        $alignments[$alncounter] = $line;
        $alncounter++;
    }
    close $sequences_fh;

    open $sequences_fh, '>', $sequences;
    foreach my $alncounter (@alignments)
    {
        print {$sequences_fh} $alncounter;
    }
    close $sequences_fh;

    print <<"ENDHTML";
    <form method ='post' action=$0>
    <center>
    <input type=hidden name='results' value='$results'>
    <input type=hidden name='alns' value='$alns'>
    <input type=hidden name='sequences' value='$sequences'>
    <input type=hidden name='prid' value='$prid'>
    <input type=hidden name='email' value='$email'>
    <input type=submit name='button' value=' Click here to view the results. '>
    </FORM>
    </center>
    </body>
    </html>
ENDHTML
}

#Results page for defined organism.#
elsif ( $query->param('button') eq ' Click here to view the results. ' )
{
    my $results   = $query->param('results');
    my $alns      = $query->param('alns');
    my $sequences = $query->param('sequences');
    my $prid      = $query->param('prid');
    my $email     = $query->param('email');
    my $fnaln     = '/web/results/final_alns/multalign/' . $prid . '.aln';
    my $ph        = '/web/results/trees/phs/' . $prid . '.ph';
    my $phb       = '/web/results/trees/phbs/' . $prid . '.phb';
    my $drawntree = '/web/results/trees/drawn/' . $prid . '.png';

    $/ = "\n\n";
    open $sequences_fh, '<', $sequences;
    while ( $line = <$sequences_fh> )
    {
        $sequences_number++;
    }
    close $sequences_fh;

    my $timeout = 60;

    print "<!--\n";
    $SIG{ALRM} = sub { print ".\n"; alarm $timeout; };
    alarm $timeout;

    $email =~ s/([\@])/\\$1/;
    if ( $sequences_number > 3 )
    {

        #Alignment, NJ tree plotting and bootstrapping.#
        tcoffee( $sequences, $fnaln, $email );
        `clustalw -INFILE=$fnaln -OUTFILE=$ph -tree`;
`clustalw -INFILE=$fnaln -OUTFILE=$phb -bootstrap=1000 -bootlabels=node`;
        `rm *.dnd`;
        `mv /web/results/final_alns/multalign/$prid.ph /web/results/trees/phs/`;
`mv /web/results/final_alns/multalign/$prid.phb /web/results/trees/phbs/`;

        open $tree_fh, '<', $phb;
        $/ = "\n";
        while ( $line = <$tree_fh> )
        {
            $line =~ s/TRICHOTOMY//;    #Remove the 'TRICHOTOMY' word.
            $line =~ s/-\d[.]\d*/0/
              ; #Fix NJ negative branch length artifact, by setting those numbers to zero.
            if ( $line =~ /^(\d+):/ )
            {
                $bvdiv = $1 / 10.0;
                $line  = $` . $bvdiv . ":" . $';
            }
            $tree[$yacounter] = $line;
            $yacounter++;
        }
        close $tree_fh;
        open $tree_fh, '>', $phb;
        foreach my $yacounter (@tree)
        {
            print {$tree_fh} $yacounter;
        }
        close $tree_fh;
        `./Pinda.R $drawntree $phb > /web/parsing/$prid.tmp`
          ;    #Visually draw the tree.

        alarm 0;
        open $tree_for_parsingfh, '<', "../parsing/$prid.tmp";
        parser();
        close $tree_for_parsingfh;
        @candidate = sort @candidate;
        @candidate = reverse @candidate;
        print <<"ENDHTML";
        -->
        <center><br><font size='3' face='Georgia' color='330033'>
        <a href=../results/final_alns/multalign/$prid.aln>T-Coffee Alignment</a> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
        <a href=../results/trees/phs/$prid.ph>NJ Tree</a> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
        <a href=../results/trees/phbs/$prid.phb>Bootstrapped NJ Tree</a>
        </font>
        <br><br>
        <table border='1'>
        <tr bgcolor=FFFF66><th><center>Possible duplications of <a href=http://www.uniprot.org/uniprot/$starting_point>$starting_point</a>.</center></th>
        <th><center>P</center></th></tr>
ENDHTML
        foreach $can (@candidate)
        {
            if ($can !~ /$starting_point/ && $can =~ /(\d?.?\d+) (\w+)/)
            {
                print "<tr><td><center><a href=http://www.uniprot.org/uniprot/$2>$2</a></center></td><td align=left>$1</td></tr>";
            }
        }
        print "</table>";
        print "<img src='../results/trees/drawn/$prid.png'>";
    }
    else
    {
        if ( $sequences_number == 3 )
        {
            print <<"ENDHTML";
            <font size='3' face='Georgia' color='330033'><br><br>
            The number of sequences for this particular organism is three.
            <br><b>The phylogenetic tree cannot be bootstrapped.</b></font>
ENDHTML
            tcoffee( $sequences, $fnaln, $email );
            `clustalw -INFILE=$fnaln -OUTFILE=$ph -tree`;
`mv /web/results/final_alns/multalign/$prid.ph /web/results/trees/phs/`;
            `rm *.dnd`;
            `./Pinda.R $drawntree $ph`;
            print <<"ENDHTML";
            <center><br><font size='3' face='Georgia' color='330033'><br>
            <a href=../results/final_alns/multalign/$prid.aln>T-Coffee Alignment</a> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
            <a href=../results/trees/phs/$prid.ph>NJ Tree</a></font>
            <br><br><center><img src='../results/trees/drawn/$prid.png'></center>
ENDHTML
        }
        if ( $sequences_number < 3 )
        {
            alarm 0;
            print <<"ENDHTML";
            <font size='3' face='Georgia' color='330033'><br><br>
            The number of sequences for this particular organism is less than three.<br><b>
            A phylogenetic tree cannot be plotted.</b></font>
ENDHTML
            if ( $sequences_number > 1 )
            {
                tcoffee( $sequences, $fnaln, $email );
                print <<"ENDHTML";
                <br><font size='3' face='Georgia' color='330033'><br>
                <a href=../results/final_alns/multalign/$prid.aln>T-Coffee Alignment</a></font>
ENDHTML
            }
        }
    }
    print "</center>\n</body>\n</html>";

}

sub tcoffee    #Pretty self-explanatory.
{
    `export HOME_4_TCOFFEE="/web/t-coffee/" ;
    export DIR_4_TCOFFEE="/web/t-coffee/dir_4_t-coffee/" ;
    export CACHE_4_TCOFFEE="/web/t-coffee/dir_4_t-coffee/cache/" ;
    export TMP="/web/t-coffee/dir_4_t-coffee/tmp/" ;
    export NO_ERROR_REPORT_4_TCOFFEE='1' ;
    export NO_WARNING_4_TCOFFEE='1' ;
    t_coffee -infile $_[0] -output=clustal,fasta_aln -outfile $_[1] -proxy -email=$_[2]`;
    return;
}

sub parser
{
    $/ = "\n\n";
    while ( $line = <$tree_for_parsingfh> )
    {
        if ( $line =~ /\_\_\_(\w+)\_\_\_/ )
        {
            $starting_point = $1;
            $uni_ids[$unicounter] = $1;
            $unicounter++;
        }
        elsif ( $line =~ /"(\w{6})"/ )
        {
            if ( $1 !~ /_/ )
            {
                $uni_ids[$unicounter] = $1;
                $unicounter++;
            }
        }
        if ( $' =~ /"(\w{6})"/ && $1 !~ /_/ )
        {

            $uni_ids[$unicounter] = $1;
            $unicounter++;
        }
        $parsing_lines[$linocounter] = $line;
        $linocounter++;
    }
    foreach $line (@parsing_lines)
    {
        if ( $line =~ /$starting_point/ )
        {
            if ( $line =~ /parts\$X(\w+)/ )
            {
                $star_seq[$sscounter] = $1;
                $search_id = $star_seq[$sscounter];
                $sscounter++;
            }
            foreach $line (@parsing_lines)
            {
              LOOP0:
                if ( $line =~ /$search_id"/ )
                {
                    if ( $line =~ /parts\$X(\w+)/ )
                    {
                        $star_seq[$sscounter] = $1;
                        $search_id = $star_seq[$sscounter];
                        $sscounter++;
                        goto LOOP0;
                    }
                }
            }
        }
    }
    foreach $uni (@uni_ids)
    {
        $sscounter = 0;
        foreach $line (@parsing_lines)
        {
            if ( $line =~ /$uni/ )
            {
                if ( $line =~ /parts\$X(\w+)/ )
                {
                    $compare_seq[$sscounter] = $1;
                    $search_id = $compare_seq[$sscounter];
                    $sscounter++;
                }
                foreach $line (@parsing_lines)
                {
                  LOOP:
                    if ( $line =~ /$search_id"/ )
                    {
                        if ( $line =~ /parts\$X(\w+)/ )
                        {
                            $compare_seq[$sscounter] = $1;
                            $search_id = $compare_seq[$sscounter];
                            $sscounter++;
                            goto LOOP;
                        }
                    }
                }
            }
        }
        foreach $node (@compare_seq)
        {
            if ( $node =~ /(\d+)\_?/ )
            {
                $ginomenon *= $1 / 100;
            }
            foreach $star_node (@star_seq)
            {
                if ( $node eq $star_node )
                {
                    foreach $star_node (@star_seq)
                    {
                        if ( $star_node ne $node && $star_node =~ /(\d+)\_?/ )
                        {
                            $ginomenon *= $1 / 100;
                        }
                        if ( $star_node eq $node )
                        {
                            goto EXIT0;
                        }
                    }
                  EXIT0:
                    if ($uni !~ /$starting_point/)
                    {
                        $candidate[$cancounter] = $ginomenon . " " . $uni;
                        $cancounter++;
                    }
                    goto EXIT;
                }
            }
        }
        $ginomenon = 1;
        foreach $node (@compare_seq)
        {
            if ( $node =~ /(\d+)\_?/ )
            {
                $ginomenon *= $1 / 100;
            }
        }
        foreach $star_node (@star_seq)
        {
            if ( $star_node =~ /(\d+)\_?/ )
            {
                $ginomenon *= $1 / 100;
            }
        }
        $candidate[$cancounter] = $ginomenon . " " . $uni;
        $cancounter++;
      EXIT:
        $ginomenon   = 1;
        @compare_seq = ();
    }
    return @candidate, $starting_point;
}