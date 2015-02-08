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


######################################################################
######################################################################


library('getopt')
#loadProject("../src/r/","pcaAnalysis.r")
spec = matrix(c(
    'fileLocation','i', 1,"character",
    'catagories','c', 1,"character",
    'cols','b', 1,"character",
    'print','p',0,"logic"
    ),byrow=TRUE,ncol=4)
args=getopt(spec)
if(is.null(args$print)){prints=FALSE}else{prints=TRUE}
rootFile<-args$fileLocation
catFile<-args$catagories
column<-as.numeric(strsplit(args$cols,",")[[1]])


values<-loadData(rootFile)
cats<-t(read.table(catFile))
data<-values$data
normData<-standPCAPrep(data,"colQn")
pcs<-prcomp(t(normData))
if(prints==TRUE){
    X11()
    print(args$print)
    plotPCs(pcs,column,normData,cats)
    locator(1)}
if(prints==FALSE){
    print("FALSE")
    write.table(pcs$rotation[,column],stdout(),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")}



