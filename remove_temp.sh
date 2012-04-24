#!/bin/bash

########################################################
#Bash script for temporary files removal after 10 days.#
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

remove_temp_files ()
{
#########################################################################
#Find files older than 10 days in the current directory and remove them.#
#########################################################################
   find . -mtime +11 -name "*" -exec rm {} \;
}

cd /var/www/Pinda/outs/blast/
remove_temp_files

cd /var/www/Pinda/outs/psiblast/
remove_temp_files

cd /var/www/Pinda/parsing/
remove_temp_files

cd /var/www/Pinda/results/final_alns/multalign/
remove_temp_files

cd /var/www/Pinda/results/final_alns/multalign/conf/
remove_temp_files

cd /var/www/Pinda/results/trees/drawn/
remove_temp_files

cd /var/www/Pinda/results/trees/phs/
remove_temp_files

cd /var/www/Pinda/results/trees/phbs/
remove_temp_files

cd /var/www/Pinda/results/trees/zips/
remove_temp_files

cd /var/www/Pinda/seq/final_seq/
remove_temp_files

cd /var/www/Pinda/tmps/blast/
remove_temp_files

exit 0
