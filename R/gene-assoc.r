#!/usr/bin/env R
#
# This file is part of peakAnalysis,
# http://github.com/alexjgriffith/alpha-score/, 
# and is Copyright (C) University of Ottawa, 2015. It is Licensed under 
# the three-clause BSD License; see LICENSE.txt.
# Author : Alexander Griffith
# Contact: griffitaj@gmail.com
#
#######################################################################
######################################################################
#
# Usage:
#
# fileLocation="~/Dropbox/UTX-Alex/jan/"
# bedData<-read.delim(paste(fileLocation,"combined_sorted.bed",sep=""),header=0)
# geneList<-read.delim(paste(fileLocation,"hg19.RefSeqGenes.csv",sep=""))
# genes<-geneAssociation(bedData,geneList,  c(50000,0,0,0))
# write.table(cbind(bedData,unlist(t(lapply(genes, function(x) {if(identical(x,character(0))){"None"} else{x}})))) ,"combined_tagged_genes.bed" ,col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")

library(parallel)

member<-function(posibList,mem,test="default"){
    fun<-switch(test,default=function(x,y){all(y %in% x)})
    sapply(posibList,fun,mem)}

memberWrapper<-function(bedData,lis ,n=4,sep="-"){
    member(strsplit(as.character(bedData[,n]),sep),lis)}

geneAssoc<-function(point,geneList,bounds){
    # Currently Returns the enhancer locations
    peak<-(as.numeric(point[[2]])+as.numeric(point[[3]]))/2
    tss<-geneList$txStart
    ess<-geneList$txEnd
    chrom<-geneList$chrom
    a<-which(chrom==as.character(point[[1]]))
    b<-a[which(tss[a]-bounds[1]<peak)]
    r<-b[which(ess[b]+bounds[3]>peak)]
    e<-r[c(which(tss[r]-bounds[2]>peak),which(ess[r]+bounds[4]<peak))]
    e}

leneAssoc<-function(point,geneList,bounds){
    peak<-(as.numeric(point[[2]])+as.numeric(point[[3]]))/2
    tss<-geneList$txStart
    ess<-geneList$txEnd
    chrom<-geneList$chrom
    a<-which(chrom==as.character(point[[1]]))
    b<-a[which(tss[a]-bounds[1]<peak)]
    r<-b[which(ess[b]+bounds[3]>peak)]
    locations<-r[c(which(tss[r]-bounds[2]>peak),which(ess[r]+bounds[4]<peak))]
    t<-geneList$txStart[unlist(locations)]
    gene<-as.character(geneList$name[locations[order(abs(t-peak))]])
    if(identical(gene,character(0))){"None"}else(gene)}

minGene<-function(x,y){
    x[order(abs(x-y))]}

geneAssociation<-function(bedData,geneList,bounds,n=FALSE){
    if(is.logical(n)){
        peak<-(as.numeric(bedData[,2])+as.numeric(bedData[,3]))/2    
        locations<-apply(bedData,1,geneAssoc,geneList,bounds)
        t<-lapply(n,locations,function(locations) geneList$txStart[unlist(locations)])
        
        tssGenes<-Map(minGene,t,peak)
        a<-lapply( tssGenes, function(x){as.character(geneList$name2[x])})}
    else{
        cs<-makeForkCluster(n,renice=0)
        a<-parApply(cs,bedData,1,leneAssoc,geneList=geneList,bounds=bounds)
    stopCluster(cs)}
    return (a)}
