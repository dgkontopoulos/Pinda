#!/usr/bin/Rscript

##########################################
#Author: Dimitrios - Georgios Kontopoulos#
##########################################

######################################
#Script in the R programming language#
# for drawing and parsing NJ-trees,  #
#       bootstrapped or not          #
######################################

###################################################################### 
#This program is free software: you can redistribute it and/or modify#
#it under the terms of the GNU General Public License as             #
#published by the Free Software Foundation, either version 3 of the  #
#License, or (at your option) any later version.                     #
#                                                                    #
#For more information, see http://www.gnu.org/licenses/.             #
######################################################################

args <- commandArgs(TRUE)
#r <- getOption('repos')
#r['CRAN'] <- 'http://cran.cc.uoc.gr/'
#options(repos=r)
#install.packages(c('ape','ade4') lib='../R')
if (args[1] == "-parser") {
    library(ade4, lib.loc = "../R/")  #load ADE4 package
    tree <- newick2phylog(scan(args[2], what = ""), add.tools = FALSE)
    treestructure <- tree[4]
    print(treestructure)
} else if (args[1] == "-lengths_1") {
    library(ade4, lib.loc = "../R/")
    tree <- newick2phylog(scan(args[2], what = ""), add.tools = FALSE)
    treestructure <- tree[3]
    print(treestructure)
} else if (args[1] == "-lengths_2") {
    library(ade4, lib.loc = "../R/")
    tree <- newick2phylog(scan(args[2], what = ""), add.tools = FALSE)
    treestructure <- tree[2]
    print(treestructure)
} else {
    library(ape, lib.loc = "../R/")  # load APE package
    png(filename = args[1], width = 750, height = 550, units = "px", 
        bg = "transparent")  # set output file
    tree <- read.tree(args[2])  # read tree
    plot(tree, use.edge.length = TRUE, show.node.label = TRUE, 
        font = 2, no.margin = TRUE, node.pos = 2)  # plot tree
    add.scale.bar(length = 0.05, col = "red", lcol = "red")  # add scale bar
    dev.off()  #quit
}
