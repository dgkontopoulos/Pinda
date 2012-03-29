#!/bin/bash

###############################################
#Bash script to keep the databases up to date.#
###############################################

##########################################
#Author: Dimitrios - Georgios Kontopoulos#
##########################################

###################################################################### 
#This program is free software: you can redistribute it and/or modify#
#it under the terms of the GNU General Public License as             #
#published by the Free Software Foundation, either version 3 of the  #
#License, or (at your option) any later version.                     #
#                                                                    #
#For more information, see http://www.gnu.org/licenses/.             #
######################################################################

##########################################################################
#Get the Swiss-Prot db, format it and also copy it to the UniProt folder.#
##########################################################################
mkdir /usr/local/databases/Swissprot.1/
cd /usr/local/databases/Swissprot.1/
nice -n +19 wget ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.fasta.gz

################################################################
#If wget does not finish properly, mail someone to take a look.#
################################################################
if [ "$?" -ne "0" ]; then
   echo 'The Swiss-prot database failed to download properly.' > /tmp/tmp
   mail dimikont2@mbg.duth.gr -s 'DB fetching error' < /tmp/tmp
   rm /tmp/tmp
else
   nice -n +19 gunzip uniprot_sprot.fasta.gz
   mkdir /usr/local/databases/UniProt.1/
   cp uniprot_sprot.fasta ../UniProt.1/
   nice -n +19 makeblastdb -in uniprot_sprot.fasta -dbtype prot -parse_seqids
   rm uniprot_sprot.fasta
   
   ##########################################################################
   #Get the TrEMBL db, concatenate it with the Swiss-Prot one and format it.#
   ##########################################################################
   cd ../UniProt.1/
   nice -n +19 wget ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.fasta.gz
   
   ################################################################
   #If wget does not finish properly, mail someone to take a look.#
   ################################################################
   if [ "$?" -ne "0" ]; then
      echo 'The TrEMBL database failed to download properly.' > /tmp/tmp
      mail dimikont2@mbg.duth.gr -s 'DB fetching error' < /tmp/tmp
      rm /tmp/tmp
   else
      nice -n +19 gunzip uniprot_trembl.fasta.gz
      cat uniprot_sprot.fasta uniprot_trembl.fasta > UniProt.fasta
      rm -rf uniprot_sprot.fasta uniprot_trembl.fasta
      nice -n +19 makeblastdb -in UniProt.fasta -dbtype prot -parse_seqids
      rm UniProt.fasta
      
      ################################################
      #Check for running processes (blastp/psiblast).#
      ################################################
      while true
      do
         if (ps ax | grep -v grep | grep blastp) > /dev/null && ( ps ax | grep -v grep | grep psiblast) > /dev/null;
         then
            sleep 60
         else
            rm -rf /usr/local/databases/Swissprot/ && mv /usr/local/databases/Swissprot.1/ /usr/local/databases/Swissprot/
            rm -rf /usr/local/databases/UniProt/ && mv /usr/local/databases/UniProt.1/ /usr/local/databases/UniProt/
            break
         fi
      done
   fi
fi

####################################
#Get the nt database and format it.#
####################################
mkdir /usr/local/databases/nt.1/
cd /usr/local/databases/nt.1/
nice -n +19 wget ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nt.gz

################################################################
#If wget does not finish properly, mail someone to take a look.#
################################################################
if [ "$?" -ne "0" ]; then
   echo 'The nt database failed to download properly.' > /tmp/tmp
   mail dimikont2@mbg.duth.gr -s 'DB fetching error' < /tmp/tmp
   rm /tmp/tmp
else
   nice -n +19 gunzip nt.gz
   mv nt nt.fasta
   nice -n +19 makeblastdb -in nt.fasta -dbtype nucl -parse_seqids
   rm nt.fasta
   
   while true
   do
      if (ps ax | grep -v grep | grep blastn) > /dev/null;
      then
         sleep 60
      else
         rm -rf /usr/local/databases/nt/ && mv /usr/local/databases/nt.1/ /usr/local/databases/nt/
         break
      fi
   done
fi
exit 0
