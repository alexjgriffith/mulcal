#!/usr/bin/env R
#
# This file is part of peakAnalysis,
# http://github.com/alexjgriffith/alpha-score/, 
# and is Copyright (C) University of Ottawa, 2015. It is Licensed under 
# the three-clause BSD License; see LICENSE.txt.
# Author : Alexander Griffith
# Contact: griffitaj@gmail.com

#' @import parallel
library(parallel)

#' loadBedFile
#' 
#' Calls file_length and read_bed from hg19Height.c and returns a 3 column bed file
#' @param file the location of the bed file of interest
#' @return data.frame $chro $start $end
#' @export
#' @template authorTemplate
#' @examples
#'
#'  data<-loadBedFile(file="filelocation")
#' 
loadBedFile<-function(file){
    file<-normalizePath(file)
    if(!file.exists(file)){
        sprintf("Can't find file %s.",file)}
    else{
        fileLength<-as.integer(.C("file_length",file,stringLength=integer(1))[[2]])
        results<-.C("read_bed",file,chro=character(fileLength),start=integer(fileLength),end=integer(fileLength))
        data.frame(chro=results$chro,start=results$start,end=results$end)
    }}

#' convertChroms
#' 
#' Calls valueChromosome from hg19Height.c which converts a chromosme rank to a string
#' @seealso For the reverse see \code{\link{rankChroms}}
#' @param lis a list of chromosmes to
#' @template authorTemplate
#' @export
convertChroms<-function(lis){
    n<-length(lis)
    .C("valueChromosomes",character(n),as.integer(n), lis)[[1]]}

#' rankChroms
#' 
#' Calls rankChromosome from hg19Height.c which converts a chromosme string to its equivelent rank
#' @seealso For the reverse see \code{\link{convertChroms}}
#' @template authorTemplate
#' @export
rankChroms<-function(lis){
    n<-length(lis)
    .C("rankChromosomes",as.character(lis),as.integer(n), integer(n))[[3]]}

#' getPileUp
#'
#' A minimal wrapper for the pileup function from hg19Height.c
#' returns a integer vector of computed read pileups
#'  @export
#' @template authorTemplate
#' @param file The raw data file of interest
#' @param bed The preloaded bed infromation including bed$start and bed$end
#' @param chroms a list of chromosomes whos string values have been repalced with ranks
#' @param peakLength The length of the bed data provided
getPileUp<-function(file,bed,chroms,peakLength){
    results<-.C("pileup",file,chroms,bed$start,bed$end,as.integer(peakLength),score=integer(peakLength))
    results$score}

#' PileUp
#'
#' Generates a pile up matrix from a unified set of peaks and a list of raw data sets
#' @export
#' @template authorTemplate
#' @param data A preloaded bed data.frame which includes slots $chro $start $end
#' @param rawdata a list of raw data files
#' @param n the number of nodes to use. If 0 then the parrallel package is not used
pileUp<-function(data,rawdata,n=0){
    for(file in rawdata){
        if(!file.exists(file)){
            sprintf("Can't find file %s.",file)
            return}}
    peakLength<-length(data$chro)
    chroms<-rankChroms(data$chro)
    if(n>0){
        cs<-makeForkCluster(n,renice=0)
        ret<-matrix(unlist(parLapply(cs,rawdata,getPileUp,data,chroms,peakLength)),nrow=peakLength)
        stopCluster(cs)}
    else{
    ret<-matrix(unlist(lapply(rawdata,getPileUp,data,chroms,peakLength)),nrow=peakLength)}
    ret}

#' hg19Sort
#'
#' Reorders the chro factor in data to that which is outputed by bwa
#' @template authorTemplate
#' @export
hg19Sort<-function(data){
    neworder<-levels(data$chro)[c(1,12,16,17,18,19,20,24,21,22,2,3,4,5,6,7,8,9,10,11,13,23,14,15,25)]
    data$chro<-factor(data$chro,neworder)
    data<-data[with(data,order(chro,start)),]
    data}

#' getData
#'
#' loads data from file returns a sorted set with a height matrix
#'
#' @export
#' @template authorTemplate
#' @param file peak file
#' @param rawdata a list of raw data files
#' @param n the number of nodes used default=0
getData<-function(file,rawdata,n=0){
    data<-hg19Sort(loadBedFile(file))
    score<-pileUp(data,rawdata,n=n)
    c(data=data,score=score)}

#' mainProf
#'
#' A simple test function for speed
#'
#' @export
#' @template authorTemplate
mainProf<-function(file="/home/agriffith/Dropbox/UTX-Alex/jan/combined_sorted.bed",
                   rawdata=rep(file,22),
                   n=0){
    Rprof("temp.prof")
    data<-hg19Sort(loadBedFile(file))
    score<-pileUp(data,rawdata,n=n)
    print(str(score))
    Rprof(NULL)
    summaryRprof("temp.prof")}

testHeights<-function(){
    file<-"/home/griffita/Dropbox/UTX-Alex/jan/combined_sorted.bed"

    cats<-read.table("/home/griffita/Dropbox/UTX-Alex/jan/catagories")
    prefix<-"/mnt/brand01-00/mbrand_analysis/data_sets/"
    suffix<-"_sorted.bed"
    rawdata<-apply(cats,1,function(x){paste(prefix,x,"/",x,suffix,sep="")})

    data<-hg19Sort(loadBedFile(file))
    score<-pileUp(data,rawdata,n=22)
    temp<-cor(score)
    rownames(temp)<-t(cats)
    colnames(temp)<-t(cats)
    pdf("test.pdf")
    plot(hclust(dist(temp)),hang=-1)
    dev.off()}
