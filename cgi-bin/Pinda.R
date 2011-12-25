#!/usr/bin/Rscript

######################################
#Script in the R programming language#
# for drawing NJ-trees, bootstrapped #
#               or not               #
######################################
args <- commandArgs(TRUE)
#r <- getOption("repos")
#r["CRAN"] <- "http://cran.stat.ucla.edu"
#options(repos=r)
#install.packages(c('ape','ade4') lib="../R")
if(args[1] == '-parser')
{
	library(ade4, lib.loc='../R/') #load ADE4 package
	tree2 <- newick2phylog(scan(args[2], what = ""))
	treestructure <- tree2[4]
	print(treestructure)
}else
{
	library(ape, lib.loc='../R/') # load APE package
	png(filename = args[1], width = 750, height = 550, units ="px", bg = "transparent") # set output file
	tree <- read.tree(args[2]) # read tree
	plot(tree, use.edge.length = TRUE, show.node.label = TRUE, font = 2, no.margin = TRUE, node.pos = 2) # plot tree
	add.scale.bar(length = 0.05, col = "red", lcol = "red") # add scale bar
	dev.off() #quit
}