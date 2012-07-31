#!/bin/bash

########################################################
#Bash script to check for remaining db temporary files.#
########################################################

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

#####################################################
#Check for remaining database temporary directories.#
#####################################################

set -e

touch /tmp/whut
if [ -d "/usr/local/databases/Swissprot.1" ]
then
   echo "Two days after the last update, the Swissprot.1 directory is still there.
   Did a power failure take place? Please, take a look." >> /tmp/whut
fi

if [ -d "/usr/local/databases/UniProt.1" ]
then
   echo "Two days after the last update, the UniProt.1 directory is still there.
   Did a power failure take place? Please, take a look." >> /tmp/whut
fi

if [ -d "/usr/local/databases/nt.1" ]
then
   echo "Two days after the last update, the nt.1 directory is still there.
   Did a power failure take place? Please, take a look." >> /tmp/whut
fi

if [ -s "/tmp/whut" ]
then
   mail Pinda -s 'DB temp files are still there!' < /tmp/whut
fi

rm /tmp/whut
exit 0
