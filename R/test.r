#!/data/binaries/R-3.1.0/bin/Rscript
#
# This file is part of peakAnalysis,
# http://github.com/alexjgriffith/alpha-score/, 
# and is Copyright (C) University of Ottawa, 2015. It is Licensed under 
# the three-clause BSD License; see LICENSE.txt.
# Author : Alexander Griffith
# Contact: griffitaj@gmail.com
#
#
######################################################################
######################################################################

loadProject<-function(projectDirectory,files){
    tempDir<-getwd()
    args <- commandArgs(trailingOnly = F)
    scriptDir<-dirname(sub("--file=","",args[grep("--file",args)]))
    setwd(scriptDir)
    for(file in files){
        source(paste(projectDirectory,file,sep=""),chdir=TRUE)}
    setwd(tempDir)}

loadProject("../src/r/", c("pcaAnalysis.r","html.r"))
