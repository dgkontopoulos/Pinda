#!/usr/bin/perl

####################
#$VERSION = '0.02';#
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

At Pinda's starting page, you may enter a Protein or DNA sequence. After
clicking the "Continue" button, Pinda will take you to the next page,
where you may choose the organism and the database within which the
duplication analysis will take place. You also have to enter your email
address, so that you can be notified after the job is complete. Past
this point, your task has entered the queue. You will receive its
results via email, including a possible duplications table, a multiple
sequence alignment and a dendrogram.

=head1 FLOWCHART

http://orion.mbg.duth.gr/Pinda/flowchart.png

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

use feature qw(say);

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
<img src='http://orion.mbg.duth.gr/Pinda/pindalogo.png' width=354
height=74 alt='Pipeline for INtraspecies Duplication Analysis'></a><br><br>
<font size=3 face='Georgia' color='330033'>
<i>Pinda is a Web service aiming to facilitate detection of<br>
specific gene duplications in an organism species of choice.</i>
</font>
<p style='width: 500px; text-align:center;margin-bottom:1px;margin-top:1px'>
<hr/>
</p>
</center>
<div style="position:fixed;bottom:0;width:100%">
<div style="width:300;margin:0px auto;">
<hr /><center>
<font size='3' face='Georgia' color='330033'>
<a href='http://orion.mbg.duth.gr/Pinda/documentation.html'>Documentation</a>
| 
<a href="https://github.com/dgkontopoulos/Pinda/">Source Code</a>
| Built with <a href="http://www.perl.org/">Perl</a></font></center>
</div></div>
ENDHTML

###############################
#Starting page with input box.#
###############################
if ( !$query->param )
{
    ###########################################
    #Get the timestamp for the database files.#
    ###########################################
    print <<"ENDHTML";
    <center>
    <br><br>
    <form id = 'submit' method = 'post' action=http://orion.mbg.duth.gr/Pinda/cgi/Pinda.cgi>
    <fieldset class="fieldset-auto-width"><legend><font size='2'
    face='Georgia' color='330033'>Sequence Input</font></legend>
    <font size='2' face='Georgia' color='330033'>
    <p style='width: 570px; text-align:left;margin-bottom:1px;margin-top:1px'>
    <b>Please enter a sequence.</b><br>
    </font><center>
    <textarea name='sequence' rows=7 cols=70 style='font-family:courier new'
    ></textarea></center></fieldset><br><br>
    <input id='submit' TYPE='submit' value=' Continue '></form></center>
    
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
    if ( $string2 =~ /^[A|C|G|T]+$/i )
    {
        $db = $NT;
    }
    elsif ( $string2 =~ /^[R|H|K|D|E|S|T|N|Q|C|U|G|P||A|V|I|L|M|F|Y|W]+$/i )
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

    my $timezone = `date`;
    if ( $timezone =~ /EEST/ )
    {
        $timezone = 'EEST';
    }
    else
    {
        $timezone = 'EET';
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
    $SIG{ALRM} = sub { say "."; alarm $timeout; };
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
        say "<!--";
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
        say "-->";
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
    <fieldset class="fieldset-auto-width"><legend><font size='2' face='Georgia' color='330033'>
    Organism Selection </font></legend>
ENDHTML

        if ( $mikos == 0 )
        {
            print <<"ENDHTML";
		<font size='2' face='Georgia' color='330033'>
		<p style='width: 270px; text-align:left;margin-bottom:1px;margin-top:1px'>
		<center>
		No organism could be identified.<br><br>
		<b>Please enter an organism name or a
		<a href=http://www.ncbi.nlm.nih.gov/Taxonomy/>Taxonomy ID</a> here.
		</b><br>
		<center><input type='text' name='organism2' size="30" maxlength="60">
		</fieldset>
ENDHTML
        }
        else
        {
            say <<"ENDHTML";
	<p style='width: 270px; text-align:left;margin-bottom:1px;margin-top:1px'>
    <font size='2' face='Georgia' color='330033'>
    <b>Either select the correct source organism from the following list
    </b></font><br><br></p>
    <p style='width: 270px; text-align:right;margin-bottom:1px;margin-top:1px'>
    <center><select name='organism' tabindex='1' class='default'>
    <option value=''>---Select an organism---</option>
ENDHTML

            for ( 0 .. $mikos - 1 )
            {
                say "<option value='$organism2[$_]'>$organism2[$_]</option>";
            }
            print <<"ENDHTML";
		</select></center></p>
		<br><br><br><font size='2' face='Georgia' color='330033'>
		<p style='width: 270px; text-align:left;margin-bottom:1px;margin-top:1px'>
		<b>or enter an organism name or a
		<a href=http://www.ncbi.nlm.nih.gov/Taxonomy/>Taxonomy ID</a>
		here.</b><br>
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
        </font><font size='1'>(Updated: $sw_timestamp $timezone)</font><br>
        <font size='2'>
        <input type="radio" name="db" value="UniProt"> <b>UniProt</b> 
        </font><font size='1'>(Updated: $uni_timestamp $timezone)</font>
        <br><br><input type=checkbox name='lcr_filtering' value='1'>
		Disable low complexity region filtering
        <br><input type=checkbox name='masking' value='1'>
		Disable alignment masking
		</fieldset><br>
ENDHTML
    }
    else
    {
        print <<"ENDHTML";
    <form id = 'submit' method = 'post' action=http://orion.mbg.duth.gr/Pinda/cgi/Pinda.cgi>
    <fieldset class="fieldset-auto-width"><legend><font size='2' face='Georgia' color='330033'>
    Organism Selection </font></legend>
	<font size='2' face='Georgia' color='330033'>
	<p style='width: 270px; text-align:left;margin-bottom:1px;margin-top:1px'>
	<center>
	<b>Please enter an organism name or a
	<a href=http://www.ncbi.nlm.nih.gov/Taxonomy/>Taxonomy ID</a> 
	here.</b><br>
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
        </font><font size='1'>(Updated: $nt_timestamp $timezone)</font><br>
        <br><input type=checkbox name='lcr_filtering' value='1'>
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
    if ( `ls /var/www/Pinda/job_files` =~ /$prid/ )
    {
        print <<"ENDHTML";
        <br><br<br><br><br>
        <center><font size='3' face='Georgia' color='330033'>
        <b>ERROR!</b> You have already submitted this sequence to Pinda.<br><br>
        If you want to resubmit it with different parameters,<br>
        please submit it again to the 
        <a href='http://orion.mbg.duth.gr/Pinda'>starting page</a>.</center></font>
ENDHTML
        exit;
    }
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
    my $masking;
    if ( defined $query->param('masking')
        && $query->param('masking') eq 1 )
    {
        $masking = '0';
    }
    else
    {
        $masking = '1';
    }

    my $job =
        "nice -n +19 ../Pinda_exec.pl $email " . '"'
      . $organism . '"'
      . " $prid $database $lcr_filtering $one $masking";
    my $job_temp = "/var/www/Pinda/job_files/$prid";
    open my $job_fh, '>', $job_temp;
    print {$job_fh} <<"END";
#!/bin/sh
cd /var/www/Pinda/slurm_errors/	
END
    print {$job_fh} $job;
    close $job_fh;
    my $job_temp_pl = $job_temp . '.pl';
    open my $job_pl_fh, '>', $job_temp_pl;
    print {$job_pl_fh} <<"END";
#!/usr/bin/perl -w
use MIME::Lite;
system("chmod +x $job_temp");
open my \$job_fh, '<', "$job_temp";
my \$email;
my \$db;
my \$organism;
my \$prid;
while ( my \$line = <\$job_fh> )
{
	if ( \$line =~ /Pinda_exec.pl (.+) ["]/ )
	{
		\$email = \$1;
	}
	if ( \$line =~ /nt[.]fasta/ )
	{
		\$db = "nt.fasta";
	}
	elsif ( \$line =~ /Swiss/ )
	{
		\$db = "Swiss-prot";
	}
	else
	{
		\$db = "UniProt";
	}
	if ( \$line =~ /["](.+)["]\\s(\\d+)\\s/ )
	{
		\$organism = \$1;
		\$prid = \$2;
	}
}
close \$job_fh;
system("$job_temp");
if ( \$? != 0 )
{
	my \$error_data = "Pinda has run into an error while processing your job.";
	\$error_data .= "\\nWe have been notified and are looking into it.";
	send_email(\$email, \$error_data);
	restore_job_count(\$db);
	open my \$input_fh, '<', "/var/www/Pinda/tmps/blast/\$prid.tmp";
	\$/ = undef;
	my \$whole = <\$input_fh>;
	close \$input_fh;
	my \$program_output = `cd /var/www/Pinda/slurm_errors && ls -1t|tail -1`;
	my \$output = `cd /var/www/Pinda/slurm_errors && cat < \$program_output`;
	my \$error = "/tmp/error_Pinda";
	open my \$error_fh, '>', \$error;
	print {\$error_fh} "<b>Pinda has run into an error.</b>\\n\\n";
	print {\$error_fh} "<b>Database:</b> \$db\\n\\n";
	print {\$error_fh} "<b>Organism:</b> \$organism\\n\\n";
	print {\$error_fh} "<b>Input:</b> \$whole\\n\\n";
	print {\$error_fh} "<b>Pinda's Output:</b> \$output\\n\\n";
	close \$error_fh;
	system("mail Pinda -s 'Pinda Error' < /tmp/error_Pinda");
	system("rm /tmp/error_Pinda");
	system("rm /var/www/Pinda/slurm_errors/\$program_output");
}
system("rm $job_temp");
system("rm $job_temp_pl");


sub send_email
{
    my \$msg = MIME::Lite->new(
        Subject => "Pinda ERROR",
        From    => "Pinda\\\@orion.mbg.duth.gr",
        To      => "\$_[0]",
        Type    => 'text/html',
        Data    => \$_[1]
    );
    \$msg->send();
    return 0;
}

sub restore_job_count
{
	my \$job_counting = "/var/www/Pinda/running_jobs";
	my \$protein_jobs;
	my \$dna_jobs;
	open my \$job_counting_fh, '<', \$job_counting;
	local \$/ = "\\n";
	while ( my \$line = <\$job_counting_fh> )
	{
		if ( \$line =~ /Protein[:]\\s(\\d+)/ )
		{
			\$protein_jobs = \$1;
		}
		elsif ( \$line =~ /DNA[:]\\s(\\d+)/ )
		{
			\$dna_jobs = \$1;
		}
	}
	close \$job_counting_fh;
	if ( \$_[0] =~ /nt[.]fasta/ )
	{
		\$dna_jobs--;
	}
	else
	{
		\$protein_jobs--;
	}
	open \$job_counting_fh, '>', \$job_counting;
	say {\$job_counting_fh} "Protein: \$protein_jobs";
	say {\$job_counting_fh} "DNA: \$dna_jobs";
	close \$job_counting_fh;
	return 0;
}

END
    close $job_pl_fh;
    system(
"cd /var/www/Pinda/slurm_errors/ && sbatch --mail-type=FAIL $job_temp_pl > /dev/null"
    );

    my $job_counting = "../running_jobs";
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
    if ( $database =~ /nt/ )
    {
        $dna_jobs++;
    }
    else
    {
        $protein_jobs++;
    }
    open $job_counting_fh, '>', $job_counting;
    print {$job_counting_fh} "Protein: $protein_jobs\n";
    print {$job_counting_fh} "DNA: $dna_jobs\n";
    close $job_counting_fh;
    my $rank = $protein_jobs + $dna_jobs;

    my ( $pr_time, $dn_time );
    my $job_average = "/var/www/Pinda/job_times";
    open my $job_average_fh, '<', $job_average;
    local $/ = "\n";
    while ( my $line = <$job_average_fh> )
    {
        if ( $line =~ /Protein Jobs[:] \d+ Average Time[:] (\d+)/ )
        {
            $pr_time = $1;
        }
        elsif ( $line =~ /DNA Jobs[:] \d+ Average Time[:] (\d+)/ )
        {
            $dn_time = $1;
        }
    }
    close $job_average_fh;
    my $estimated_time;
    $estimated_time =
      ( ( $protein_jobs * $pr_time ) + ( $dna_jobs * $dn_time ) ) * 2;

    print <<"ENDHTML";
    <p style='width: 470px; text-align:center;margin-bottom:1px;margin-top:1px'>
    <font size='3' face='Georgia' color='330033'>
    <center><br><br><br>
    <b>Your job has entered the queue with a rank of
    <font size=5>#$rank</font>.
    <br>Results will be sent to you via email.</b></p></center>
    </p></font><font color='black'>
ENDHTML

    if ( defined $note && $note == 1 )
    {
        print <<"ENDHTML";
		<br><center><font size='2' face='Georgia' color='330033'>
		<b>Note:</b> This Taxonomy ID does not belong to a species.
		<br>We hope you know what you're doing...</font></center>
		<font color='black'>
ENDHTML
    }
    elsif ( defined $note && $note == 2 )
    {
        print <<"ENDHTML";
		<br><center><font size='2' face='Georgia' color='330033'>
		<b>Note:</b> The name that you entered does not correspond to one exact species.
		<br>We hope you know what you're doing...</font></center>
		<font color='black'>
ENDHTML
    }

    job_timer($estimated_time);
    print <<"ENDHTML";
    <br><br><br><center><font size='2'>
    [<a href="http://orion.mbg.duth.gr/Pinda">Return to submission form</a>]
    <br><br><br>
    <fieldset class="fieldset-auto-width"; style="background-color: #FFFFCC; width:300px"><legend>
    Please, do <b>NOT</b> use your browser's back button to resubmit this sequence with other parameters.
    <br><br>Rather, submit it again via the <a href="http://orion.mbg.duth.gr/Pinda">starting page</a>.</font>
    </legend></fieldset>
    <br><br><br></center>
ENDHTML
}

sub job_timer
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
    print '<br><br><font size="2" face="Courier New"><center>
    Estimated time until job completion: ';
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
            printf ' <b>and %.0f minutes</b>.', $minutes;
        }
        else
        {
            printf ' <b>and 1 minute</b>.';
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
            printf '<b>%.0f hours</b>.', $hours;
        }
        else
        {
            printf '<b>1 hour</b>.';
        }
    }
    elsif ( defined $minutes && defined $seconds )
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
