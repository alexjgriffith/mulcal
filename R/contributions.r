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
#' swCat<-c("Erythroid","K562","Jurkat","CEM","RPMI","Prima5(Brand)","Prima2","Prima5","Megakaryocite","ECFC","CD133","CD34")
#' names(swCat)<-c("eryt","k562","jurk","cem","rpmi","tall_p1","tall_p2","tall_p3","meka","ecfc","cd133","cd34")
#' tSwCat<-as.character(unlist(lapply(catagories,function(x) swCat[x])))
#' stackedContrib(data,pcs[,1],catagories,tSwCat)
#' @export
#' @template authorTemplate
stackedContrib<-function(data,
                         tags,
                         catagories,
                         swCat=catagories,
                         colors=FALSE,
                         inord=seq(length(swCat)),
                         steps=40,
                         n=10,
                         f=c(function(x){mean(x)+sqrt(var(x))*n},
                             function(x){mean(x)-sqrt(var(x))*n})
                        ,sum=FALSE){
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
    if(sum==FALSE) {c1<-t(apply(c1,1,function(x) x/sum(x)))}
    c4<-unlist(lapply(catagories2,function(cats) {unlist(apply(cbind(c1[,unlist(lapply(catagories, function(x) swCat[x]))==cats],rep(0,steps+1)),1,sum))} ))
    b2<-matrix(c4,steps+1)
    colnames(b2)<-catagories2
    if(length(inord)==1)
        rord<-order(apply(b2,2,sum))
    else
        rord<-inord
    c4<-as.numeric(b2[,rord])
    catagories2<-catagories2[rord]
    c3<-unlist(lapply(catagories2,rep,steps+1))
    dt<-data.frame(x=x,data=c4,cats=factor(c3,levels=rev(catagories2)),colors=colors[rord])
    ggplot(data=dt,aes(x=x,y=data,fill=cats,order=-as.numeric(cats)))+
        geom_area(position="stack")+
        theme(panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.border=element_blank(),axis.title.y=element_blank()) +
            scale_y_continuous(breaks=NULL)+
                scale_x_continuous(name="Principle Component Value") +
                    scale_fill_discrete(name="Cell Conditions")}

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
#' swCat<-c("Erythroid","K562","Jurkat","CEM","RPMI","Prima5(Brand)","Prima2","Prima5","Megakaryocite","ECFC","CD133","CD34")
#' names(swCat)<-c("eryt","k562","jurk","cem","rpmi","tall_p1","tall_p2","tall_p3","meka","ecfc","cd133","cd34")
#' tSwCat<-as.character(unlist(lapply(catagories,function(x) swCat[x])))
#' stackedContrib(data,pcs[,1],catagories,tSwCat)
#' @export
#' @template authorTemplate
stackedContrib2<-function(data,
                         tags,
                         catagories,
                         swCat=catagories,
                          repls=FALSE,
                         colors=FALSE,
                         inord=seq(length(unique(swCat))),
                         steps=40,
                         n=10,
                         f=c(function(x){mean(x)+sqrt(var(x))*n},
                             function(x){mean(x)-sqrt(var(x))*n})
                        ,sum=FALSE){
    boundLimit<-function(min,vector,width){
        f<-cbind(function(x){x>=min},
                 function(x){x<min+width})
        a<-lapply(f,function(f)f(vector))
        unlist(a[[1]]==a[[2]])}
    names(swCat)<-catagories
    bounds<-lapply(f,function(x) x(data))
    names(bounds)<-c("max","min")
    range<-bounds$max-bounds$min
    iter<-range/steps
    x<-seq(from=bounds$min,to=bounds$max,by=iter)
    a<-lapply(x,boundLimit,data,iter)

    #b<-lapply(a,function(loc,tag){tag[loc]},tags)
    #print(str(b))
    #bd<-data.frame(matrix(unlist(lapply(catagories,function(x)as.numeric(grep(x,tags))))),ncol=length(catagories) )
    #findKeys<-function(x,i) length(as.numeric(grep(paste(unlist(x),collapse="|"),tags[unlist(i)])))
        
    findKeys<-function(x,i)
        sum(as.numeric(
            unlist(lapply(x,function(x) length(grep(x,tags[unlist(i)]))))))
    
    
    bt<-data.frame(unlist(lapply(swCat,findKeys,a[1])))
    for(i in a[2:length(a)]){
        l<-unlist(lapply(swCat,findKeys,i))
                    bt<-cbind(bt,l)
    }
    bt<<-bt
    if(sum)
        bt<-apply(bt,2,function(x) x/sum(x))
    
    rord<-order(apply(bt,1,sum))

    c4<-as.numeric(t(bt[rord,]))
    catagories2<-catagories[rord]
    c3<-unlist(lapply(catagories2,rep,steps+1))
    dt<-data.frame(x=x,data=c4,cats=factor(c3,levels=rev(catagories2)))
    #dt

    #test.dt<<-dt
    ggplot(data=dt,aes(x=x,y=data,fill=cats,order=-as.numeric(cats)))+
        geom_area(position="stack")+
        theme(panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.border=element_blank(),axis.title.y=element_blank()) +
            scale_y_continuous(breaks=NULL)+
                scale_x_continuous(name="Principle Component Value") +
                    scale_fill_discrete(name="Cell Conditions")}

