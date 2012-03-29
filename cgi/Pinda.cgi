#!/usr/bin/perl -w

####################
#$VERSION = '0.01';#
####################

=head1 NAME

PINDA - Pipeline for INtraspecies Duplication Analysis

=head1 DESCRIPTION

A Web service aiming to facilitate detection of specific gene 
duplications in an organism species of choice.

=head1 AVAILABILITY

Pinda is located at http://orion.mbg.duth.gr/Pinda/cgi/Pinda.cgi, 
whereas the Source Code can be obtained at 
https://github.com/dgkontopoulos/Pinda.

=head1 USAGE

At Pinda's index page, you may enter a sequence, select a database 
for the analysis and submit your email address, so that you can be 
notified when the job is complete. After clicking the "Submit" button, 
Pinda will take you to the next page, where you may choose the organism 
within which the duplication analysis will take place. Past this point, 
your task has entered the queue. You will be notified after completion 
via email, including a possible duplications table, a multiple sequence 
alignment and a dendrogram.

=head1 AUTHOR

Dimitrios - Georgios Kontopoulos

=head1 LICENSE

This program is free software: you can redistribute it and/or modify 
it under the terms of the GNU Affero General Public License as 
published by the Free Software Foundation, either version 3 of the 
License, or (at your option) any later version.

For more information, see http://www.gnu.org/licenses/.

=cut

use Bio::DB::Taxonomy;
use Bio::Search::HSP::GenericHSP;
use Bio::Tools::Run::StandAloneBlastPlus::BlastMethods;
use CGI qw(:standard);
use Data::Validate::Email qw(is_email);
use File::stat;
use FreezeThaw qw(freeze thaw);
use List::MoreUtils qw(uniq);
use LWP::Simple qw(get);
use Sys::CPU;
use Time::localtime;

use strict;
use warnings;

open STDERR, '>', '/tmp/err' or die $!;

###################################
#Initializing the CGI environment.#
###################################
my $query = CGI->new;
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
<center><a href='http://orion.mbg.duth.gr/Pinda/cgi/Pinda.cgi'>
<img src='http://orion.mbg.duth.gr/Pinda/pindalogo.png'></a><br><br>
<font size=3 face='Georgia' color='330033'>
<i>Pinda is a Web service aiming to facilitate detection of<br>
specific gene duplications in an organism species of choice.</i>
</font><font size=2 face='Georgia' color='330033'>
<p style='width: 500px; text-align:center;margin-bottom:1px;margin-top:1px'>
<hr/>
</p>
</center>
<div style="position:fixed;bottom:0;width:100%">
<div style="width:300;margin:0px auto;">
<hr /><center>
<font size='3'>
<a href='http://orion.mbg.duth.gr/Pinda/documentation.html'>Documentation</a>
| 
<a href="https://github.com/dgkontopoulos/Pinda/">Source Code</a>
| Built with <a href="http://www.perl.org/">Perl</a></center>
</div></div>
ENDHTML

############################
#Index Page with input box.#
############################
if ( !$query->param )
{
    ###########################################
    #Get the timestamp for the database files.#
    ###########################################
    print <<"ENDHTML";
    <center>
    <br><br>
    <form id = 'submit' method = 'post' action=http://orion.mbg.duth.gr/Pinda/cgi/Pinda.cgi>
    <fieldset class="fieldset-auto-width"><legend><font size='2'>Sequence Input</legend>
    <font size='2' face='Georgia' color='330033'>
    <p style='width: 570px; text-align:left;margin-bottom:1px;margin-top:1px'>
    <b>Please enter a sequence.</b><br>
    </font><center>
    <textarea name='sequence' rows=7 cols=70 style='font-family:courier new'
    ></textarea></center></p></fieldset><br><br>
    <input id='submit' TYPE='submit' value=' Continue '></form></center>
    
    <script type="text/javascript">
    \$('#submit').submit(function(){
    \$('input[type=submit]', this).attr('disabled', 'disabled');
    }); 
    </script>
    </body></html>
ENDHTML
}

######################################
#After the original input is given...#
######################################
elsif ( !$query->param('button') && !$query->param('dropdown') )
{
    my $SWISSPROT = '/usr/local/databases/Swissprot/uniprot_sprot.fasta';
    my $UNIPROT   = '/usr/local/databases/UniProt/UniProt.fasta';
    my $NT        = '/usr/local/databases/nt/nt.fasta';
    my $list      = 0;
    my $list2     = 0;
    my $number    = 0;

    my (
        $accession, $db,         $dbfetch,  $hit,       $hit_check,
        $input_hit, $match_line, $org,      $organism,  $orghash,
        $tmp_fh,    @input_org,  @organism, %organisms, %seq
    );
    print '<center>';
    my $string = $query->param('sequence');

    #####################
    #Database selection.#
    #####################
    my $string2;
    if ( $string =~ /^>.*\n/ )
    {
        $string2 = $';
    }
    else
    {
        $string2 = $string;
    }
    chomp $string2;
    do
    {
        $string2 =~ s/\n//;
        $string2 =~ s/\s//;
    } while ( $string2 =~ /\n/ || $string2 =~ /\s/ );
    if ( $string2 =~ /^[A|C|G|T]+$/ )
    {
        $db = $NT;
    }
    elsif ( $string2 =~ /^[R|H|K|D|E|S|T|N|Q|C|U|G|P||A|V|I|L|M|F|Y|W]+$/ )
    {
        $db = $SWISSPROT;
    }
    else
    {
        print <<"ENDHTML";
        <br><br<br><br><br>
        <font size='3' face='Georgia' color='330033'>
        <b>ERROR!</b><br>This sequence looks neither like an amino acid one nor like a nucleotide one.<br><br>
        Please <a href='javascript:history.go(-1)'><u>go back</u></a> and enter
        a Protein or DNA sequence.
        </font>
ENDHTML
        exit;
    }

    ##########################
    #Handling empty sequences#
    ##########################
    if ( $string =~ /^\s*$/ || $string !~ /\n?\w+/ )
    {
        print <<"ENDHTML";
        <br><br<br><br><br>
        <font size='3' face='Georgia' color='330033'>
        <b>ERROR!</b> No sequence entered!<br><br>
        Please <a href='javascript:history.go(-1)'><u>go back</u></a> and enter
        a sequence.
        </font>
ENDHTML
        exit;
    }
    #####################################################
    #Making sure that sequences will be in FASTA format.#
    #####################################################
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
        ############################################
        #Handling organism declaration in proteins.#
        ############################################
        if (   ( $string =~ /OS\=(\w+)\s+(\w+)/ )
            || ( $string =~ /\[(\w+)\s+(\w+)/ ) )
        {
            $organism = $1 . q{ } . $2;
            if (   ( $' =~ /^$/ )
                || ( $' =~ /^\s+$/ ) )
            {
                print <<"ENDHTML";
                <br><br<br><br><br><br>
                <font size='3' face='Georgia' color='330033'>
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
    }
    print <<"ENDHTML";
    <div id="loading" class='unhidden'>
    <br><br<br><br><br><br>
    <center><img src="../loading.gif"></center><br>
    </div>
ENDHTML

    #########################################
    #Get process id and set BLAST parameters#
    #########################################
    my $prid = time . $$;
    my $tmp  = '../tmps/blast/' . $prid . '.tmp';
    my $out  = '../outs/blast/' . $prid . '.tmp';
    open $tmp_fh, '>', $tmp or die $!;
    print {$tmp_fh} "$string";
    close $tmp_fh or die $!;
    my $query_line;

    if ( $string =~ />.*\n/ )
    {
        $query_line = $';
        do
        {
            $query_line =~ s/\n//;
            $query_line =~ s/\s//;
        } while ( $query_line =~ /\n/ || $query_line =~ /\s/ );
    }

    ###########################
    #Avoiding browser timeout.#
    ###########################
    my $timeout = 60;
    $SIG{ALRM} = sub { print ".\n"; alarm $timeout; };
    alarm $timeout;

    ##############################
    #Get the number of cpu cores.#
    ##############################
    my $cpu_n = Sys::CPU::cpu_count();

    ######################
    #BLASTP for Proteins.#
    ######################
    if ( $db !~ /nt[.]fasta/ )
    {
        print "<!--\n";
        system(
"blastp -query $tmp -db $db -evalue 0.1 -num_threads $cpu_n -out $out -seg yes"
          ) == 0
          or die $?;

        my $blast = Bio::SearchIO->new(
            -format => 'blast',
            -file   => "$out"
        );

        my $hit_old = 0;

        #############################
        #Start parsing BLAST output.#
        #############################
        while ( my $result = $blast->next_result )
        {
            while ( my $hit = $result->next_hit )
            {
                if ( $hit->description =~ /OS\=\w+\s+\w+/ )
                {
                    if ( $hit->accession =~ /tr[|](\w+)[|]/ )
                    {
                        my $ac = $1;
                        local $/ = undef;
                        open my $out_fh, '<', $out or die $!;
                        while ( my $readline = <$out_fh> )
                        {
                            if ( $readline =~ />tr[|]$ac[|]/ )
                            {
                                $readline = $';
                                if ( $readline =~ /OS\=(\w+\s+\w+)\s+/ )
                                {
                                    $org = $1;
                                    if ( $org =~ /\n/ )
                                    {
                                        $org = $` . $';
                                    }
                                }
                            }
                        }
                        close $out_fh or die $!;
                    }
                    else
                    {
                        my $ac = $hit->accession;
                        local $/ = undef;
                        open my $out_fh, '<', $out or die $!;
                        while ( my $readline = <$out_fh> )
                        {
                            if ( $readline =~ />sp[|]$ac[|]/ )
                            {
                                $readline = $';
                                if ( $readline =~ /OS\=(\w+\s+\w+)\s+/ )
                                {
                                    $org = $1;
                                    if ( $org =~ /\n/ )
                                    {
                                        $org = $` . $';
                                    }
                                }
                            }
                        }
                    }
                    ######################################
                    #Populate the organism dropdown list.#
                    ######################################
                    $organism[$list] = $org;
                    $list++;

                    ############################################################
                    #Populate a hash with the first result from every organism.#
                    ############################################################
                    while ( my $hsp = $hit->next_hsp )
                    {
                        my $match_line = $hsp->hit_string;
                        do
                        {
                            $match_line =~ s/\n//;
                            $match_line =~ s/\s//;
                          } while ( $match_line =~ /\n/
                            || $match_line =~ /\s/ );
                        if ( $match_line eq $query_line )
                        {
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
                        else
                        {
                            $organisms{$org} = q{};
                        }
                    }
                }
            }
        }
        ###############################
        #Show each organism only once.#
        ###############################
        my @organism2 = uniq(@organism);
        my $mikos     = @organism2;
        #####################################
        #Pass the hash to the next CGI call.#
        #####################################
        $orghash = freeze %organisms;
        alarm 0;
        print "-->\n";
        if ( defined $organism )
        {
            print <<"ENDHTML";
        <font size='3' face='Georgia' color='330033'>
        <i>It seems that the source organism is <b>$organism</b>.
        <br>Is this correct?</i><br><br></font>
ENDHTML
        }
        ################################
        #Creation of the dropdown list.#
        ################################
        print <<"ENDHTML";
    <form id = 'submit' method = 'post' action=http://orion.mbg.duth.gr/Pinda/cgi/Pinda.cgi>
    <fieldset class="fieldset-auto-width"><legend><font size='2'>
    Organism Selection </legend>
ENDHTML

        if ( $mikos == 0 )
        {
            print <<"ENDHTML";
		<font size='2' face='Georgia' color='330033'>
		<p style='width: 270px; text-align:left;margin-bottom:1px;margin-top:1px'>
		<center>
		No organism could be identified.<br><br>
		<b>Please enter an organism name or a Taxonomy ID here.</b><br>
		<center><input type='text' name='organism2' size="30" maxlength="60">
		</fieldset>
ENDHTML
        }
        else
        {
            print <<"ENDHTML";
	<p style='width: 270px; text-align:left;margin-bottom:1px;margin-top:1px'>
    <font size='2' face='Georgia' color='330033'>
    <b>Either select the correct source organism from the following list
    </b></font><br><br></p>
    <p style='width: 270px; text-align:right;margin-bottom:1px;margin-top:1px'>
    <center><select name='organism' tabindex='1' class='default'>
    <option value=''>---Select an organism---</option>\n
ENDHTML

            for ( 0 .. $mikos - 1 )
            {
                print
                  "<option value='$organism2[$_]'>$organism2[$_]</option>\n";
            }
            print <<"ENDHTML";
		</select></center></p>
		<br><br><br><font size='2' face='Georgia' color='330033'>
		<p style='width: 270px; text-align:left;margin-bottom:1px;margin-top:1px'>
		<b>or enter an organism name or a Taxonomy ID here.</b><br>
		<center><input type='text' name='organism2' size="30" maxlength="60">
		</input></center></fieldset>
ENDHTML
        }
        my $sw_timestamp =
          ctime( stat('/usr/local/databases/Swissprot/')->mtime );
        my $uni_timestamp =
          ctime( stat('/usr/local/databases/UniProt/')->mtime );
        print <<"ENDHTML";
        <br><br><fieldset class="fieldset-auto-width"><legend>
         Parameters </legend>
		<p style='width: 270px; text-align:left;margin-bottom:1px;margin-top:1px'>
		<b>Please select a database.</b><br>
		<font size='2'>
        <input type="radio" name="db" value="Swiss-Prot" checked><b>Swiss-Prot</b>
        </font><font size='1'>(Updated: $sw_timestamp)</font><br>
        <font size='2'>
        <input type="radio" name="db" value="UniProt"> <b>UniProt</b> 
        </font><font size='1'>(Updated: $uni_timestamp)</font>
        <br><br>
        <input type=checkbox name='lcr_filtering' value='1'>
		Disable low complexity region filtering
        </fieldset><br>
        
ENDHTML
    }
    else
    {
        print <<"ENDHTML";
    <form id = 'submit' method = 'post' action=http://orion.mbg.duth.gr/Pinda/cgi/Pinda.cgi>
    <fieldset class="fieldset-auto-width"><legend><font size='2'>
    Organism Selection </legend>
	<font size='2' face='Georgia' color='330033'>
	<p style='width: 270px; text-align:left;margin-bottom:1px;margin-top:1px'>
	<center>
	<b>Please enter an organism name or a Taxonomy ID here.</b><br>
	<center><input type='text' name='organism2' size="30" maxlength="60">
	</fieldset>
ENDHTML
        my $nt_timestamp = ctime( stat('/usr/local/databases/nt/')->mtime );
        print <<"ENDHTML";
        <br><br><fieldset class="fieldset-auto-width"><legend>
         Parameters </legend>
		<p style='width: 270px; text-align:left;margin-bottom:1px;margin-top:1px'>
		<b>Database</b><br>
		<font size='2'>
        <input type="radio" name="db" value="nt" checked><b>nt</b>
        </font><font size='1'>(Updated: $nt_timestamp)</font><br>
        <br>
        <input type=checkbox name='lcr_filtering' value='1'>
		Disable low complexity region filtering
        </fieldset><br>
        
ENDHTML
    }

    print <<"ENDHTML";
    <br><fieldset class="fieldset-auto-width"><legend> Email Address </legend>
    <p style='width: 270px; text-align:left;margin-bottom:1px;margin-top:1px'>
    <font size='2' face='Georgia' color='330033'><b>Please enter a valid email
    address so that you can be notified upon job completion.</b></font><br></p>
    <center><input type='text' name='email' size="30" maxlength="60"></input>
    <br><font size='1'>Your email address will <b><u>NOT</u></b> be stored.
    </fieldset></font>
    <br><br><input type=submit name='dropdown' value='Submit'></form>
    <br><br><br><br></p>
ENDHTML

    if ( $db !~ /nt[.]fasta/ )
    {
        print <<"ENDHTML";
    <input type=hidden name='prid' value='$prid'>
    <input type=hidden name='db' value='$db'>
    <input type=hidden name='organisms' value='$orghash'>
ENDHTML
    }
    else
    {
        print <<"ENDHTML";
		<input type=hidden name='prid' value='$prid'>
		<input type=hidden name='db' value='$db'>
ENDHTML
    }

    print <<"ENDHTML";
    <script type="text/javascript">
    document.getElementById("loading").className = "hidden";
    </script>
    
    </body>
    </html>
ENDHTML
}
else
{
    my $email = $query->param('email');
    if ( !is_email($email) )
    {
        print <<"ENDHTML";
        <br><br<br><br><br><br>
        <center><font size='3' face='Georgia' color='330033'>
        <b>ERROR!</b> The email address does not appear to be valid!<br><br>
        Please <a href='javascript:history.go(-1)'><u>go back</u></a> and enter 
        a valid email address.</center></font>
ENDHTML
        exit;
    }

    my ( $organism, $note );
    if ( defined $query->param('organism') && $query->param('organism') ne q{} )
    {
        $organism = $query->param('organism');
    }
    elsif ( $query->param('organism2') ne q{} )
    {
        my $org_temp = $query->param('organism2');
        if ( $org_temp =~ /\d+/ )
        {
            my $dbfetch = get(
"http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=taxonomy&id=$org_temp&style=raw"
            );
            if ( $dbfetch =~ /SCIENTIFIC NAME\s+[:] (\w.+)\n/ )
            {
                $organism = $1;
                if ( $dbfetch !~ /RANK\s+[:] species/ )
                {
                    $note = 1;
                }
            }
            else
            {
                print <<"ENDHTML";
				<br><br<br><br><br>
				<center><font size='3' face='Georgia' color='330033'>
				<b>ERROR!</b> This Taxonomy ID is not valid!<br><br>
				Please <a href='javascript:history.go(-1)'><u>go back</u></a> and
				select an organism or enter a valid Taxonomy ID.</center></font>
ENDHTML
                exit;
            }
        }
        else
        {
            my $db      = Bio::DB::Taxonomy->new( -source => 'entrez' );
            my $taxonid = $db->get_taxonid($org_temp);
            my $dbfetch = get(
"http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=taxonomy&id=$taxonid&style=raw"
            );
            if ( $dbfetch =~ /SCIENTIFIC NAME\s+[:] (\w.+\n)/ )
            {
                if ( $dbfetch !~ /RANK\s+[:] species/ )
                {
                    $note = 2;
                }
            }
            elsif ( !defined $taxonid )
            {
                $note = 2;
            }
            $org_temp =~ s/\s+$//;
            $organism = $org_temp;
        }
    }
    else
    {
        print <<"ENDHTML";
        <br><br<br><br><br>
        <center><font size='3' face='Georgia' color='330033'>
        <b>ERROR!</b> No organism specified!<br><br>
        Please <a href='javascript:history.go(-1)'><u>go back</u></a> and
        specify the source organism.</center></font>
ENDHTML
        exit;
    }
    my $database = $query->param('db');
    my $one;
    if ( $database !~ /nt/ )
    {
        my %organisms = thaw $query->param('organisms');
        if ( defined $organisms{$organism} && $organisms{$organism} ne q{} )
        {
            $one = $organisms{$organism};
        }
        else
        {
            $one = 'QUERY';
        }
    }
    else
    {
        $one = 'QUERY';
    }
    my $prid = $query->param('prid');
    my $lcr_filtering;
    if ( defined $query->param('lcr_filtering')
        && $query->param('lcr_filtering') eq 1 )
    {
        $lcr_filtering = '0';
    }
    else
    {
        $lcr_filtering = '1';
    }

    my $job =
        "../Pinda_exec.pl $email " . '"'
      . $organism . '"'
      . " $prid $database $lcr_filtering $one";
    my $job_temp = "/tmp/$prid";
    open my $job_fh, '>', $job_temp;
    print {$job_fh} $job;
    close $job_fh;
    system("at -q Z -f $job_temp now");
    system("rm $job_temp");

    print <<"ENDHTML";
    <p style='width: 470px; text-align:center;margin-bottom:1px;margin-top:1px'>
    <font size='3' face='Georgia' color='330033'>
    <center><br><br><br>
    <b>Your job has entered the queue.
    <br>Results will be sent to you via email.</b></p></center>
    </p></font>
ENDHTML

    if ( defined $note && $note == 1 )
    {
        print <<"ENDHTML";
		<br><center><font size='2' face='Georgia' color='330033'>
		<b>Note:</b> This Taxonomy ID does not belong to a species.
		<br>We hope you know what you're doing...</font></center>
ENDHTML
    }
    elsif ( defined $note && $note == 2 )
    {
        print <<"ENDHTML";
		<br><center><font size='2' face='Georgia' color='330033'>
		<b>Note:</b> The name that you entered does not correspond to one exact species.
		<br>We hope you know what you're doing...</font></center>
ENDHTML
    }
}
