#!/usr/bin/env R
#
# This file is part of peakAnalysis,
# http://github.com/alexjgriffith/alpha-score/, 
# and is Copyright (C) University of Ottawa, 2015. It is Licensed under 
# the three-clause BSD License; see LICENSE.txt.
# Author : Alexander Griffith
# Contact: griffitaj@gmail.com

testZero<-function(value,width){
    shift<-value-width
    if(shift>=0)
        return (c(value-width,value+width))
    else{
        return (c(value-width-shift,value+width-shift))}}

outputFile<-function(peaks,chrdata,width=150){
    output<-mapply(function(i1,i2,width=300){
        summit<-as.character(testZero(floor(mean(as.double(chrdata[i1:i2,5]))),width) )
        filevalues<-lapply(chrdata[i1:i2,10],function(x){chrcats[which(sub("_[[:digit:]]*$","",x,perl=TRUE)==fileOrder)]})
        files<-do.call(paste,c(as.list(sort(unique(as.character(unlist(filevalues))))),sep="-"))
        chro<-chrdata[i1,1]
        c(chro,summit,files)
    },peaks[1,],peaks[2,],MoreArgs=list(width=width))
    data.frame(t(output))}


tags<-function(fileList,cats,overlapWidth,outputWidth){
    chrcats<-as.character(unlist(cats))
    temp<-as.character(do.call(rbind,lapply(fileList,read.table,header=TRUE,nrow=1))$name)
    fileOrder<-unlist(lapply(temp ,function(x){substr(x,1,nchar(x)-2)}))
    rm(temp)
    dataset<-do.call(rbind,lapply(fileList,read.table,header=TRUE))
    dataset<-dataset[with(dataset,order(chr,abs_summit)),]
    chrdata<-as.matrix(dataset)
    l<-nrow(dataset)
    t1<-dataset[1:l-1,]
    t2<-dataset[2:l,]
    b<-mapply(all,chr=t1$chr==t2$chr , summit=abs(t1$abs_summit-t2$abs_summit)<overlapWidth)
    prePeaks<-unlist(lapply(which(b==FALSE),function(x){c(x,x+1)}))
    peaks<-matrix(c(1,prePeaks[1:length(prePeaks)-1]),nrow=2)
    outputFile(peaks,chrdata,outputWidth)}
    
testTags<-function(file="test.bed"){
    cats<-read.table("/home/griffita/Dropbox/UTX-Alex/jan/catagories")
    prefix<-"/mnt/brand01-00/mbrand_analysis//peaks/october/"
    suffix<-"/combined_mock_peaks.xls"
    fileList<-apply(cats,1,function(x){paste(prefix,x,suffix,sep="")})
    data<-tags(fileList,cats,600,150)
    write.table(data,file,quote=FALSE,row.names=FALSE,col.names = FALSE)}
