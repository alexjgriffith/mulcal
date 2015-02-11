#!/usr/bin/env R
#
# This file is part of peakAnalysis,
# http://github.com/alexjgriffith/alpha-score/, 
# and is Copyright (C) University of Ottawa, 2015. It is Licensed under 
# the three-clause BSD License; see LICENSE.txt.
# Author : Alexander Griffith
# Contact: griffitaj@gmail.com
#

#' Stacked Contribution
#'
#' @examples
#' tags<-read.table("~/Dropbox/UTX-Alex/jan/combined_sorted.bed")[,4]
#' catagories<-unique(unlist(lapply(levels(tags),strsplit,"-")))
#' values<-loadHeightFile("~/Dropbox/UTX-Alex/jan/combined_heights.bed")
#' data<-values$data
#' info<-values$info
#' pcs<-as.matrix(mulcal::ascore(data,c(1:10)))
#' 
#' 
#' swCat<-c("normal","normal","abnormal","abnormal","abnormal","abnormal","abnormal","abnormal","meka","ecfc","stem","stem")
#' names(swCat)<-c("eryt","k562","jurk","cem","rpmi","tall_p1","tall_p2","tall_p3","meka","ecfc","cd133","cd34")
#' tSwCat<-as.character(unlist(lapply(catagories,function(x) swCat[x])))
#' stackedContrib(pcs[,1],catagories)
#'
#' @export
#' @template authorTemplate
stackedContrib<-function(data,catagories,swCat=catagories,steps=40,n=10,f=c(function(x){mean(x)+sqrt(var(x))*n},function(x){mean(x)-sqrt(var(x))*n})){
    boundLimit<-function(min,vector,width){
        f<-cbind(function(x){x>=min},
                 function(x){x<min+width})
        a<-lapply(f,function(f)f(vector))
        unlist(a[[1]]==a[[2]])}
    names(swCat)<-catagories
    catagories2<-unique(swCat)
    bounds<-lapply(f,function(x) x(data))
    names(bounds)<-c("max","min")
    range<-bounds$max-bounds$min
    iter<-range/steps
    x<-seq(from=bounds$min,to=bounds$max,by=iter)
    a<-lapply(x,boundLimit,data,iter)
    b<-lapply(a,function(loc,tag){tag[loc]},tags)
    ca<-lapply(catagories, function(y) {unlist(lapply(b,function(x,y) length(grep(y,x)),y))} )
    c1<-matrix(unlist(ca),steps+1)
    colnames(c1)<-catagories
    c2<-t(apply(c1,1,function(x) x/sum(x)))
    c4<-unlist(lapply(catagories2,function(cats) {unlist(apply(cbind(c2[,unlist(lapply(catagories, function(x) swCat[x]))==cats],rep(0,steps+1)),1,sum))} ))
    c3<-unlist(lapply(catagories2,rep,steps+1))
    dt<-data.frame(x=x,data=c4,cats=c3)
    ggplot(data=dt,aes(x=x,y=data,fill=cats))+geom_area(position="stack")}
