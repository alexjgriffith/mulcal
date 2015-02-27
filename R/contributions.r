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
#' stackedContrib(pcs[,1],catagories)
#'
#' @export
#' @template authorTemplate
stackedContrib<-function(data,catagories,swCat=catagories,steps=40,n=10,f=c(function(x){mean(x)+sqrt(var(x))*n},function(x){mean(x)-sqrt(var(x))*n}),sum=FALSE){
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
    rord<-order(apply(b2,2,sum))
    c4<-as.numeric(b2[,rord])
    catagories2<-catagories2[rord]
    c3<-unlist(lapply(catagories2,rep,steps+1))
    dt<-data.frame(x=x,data=c4,cats=factor(c3,levels=rev(catagories2)))
    ggplot(data=dt,aes(x=x,y=data,fill=cats,order=-as.numeric(cats)))+geom_area(position="stack")+
        theme(panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.border=element_blank(),axis.title.y=element_blank()) +
            scale_y_continuous(breaks=NULL)+
            scale_x_continuous(name="Principle Component Value") +
             scale_fill_discrete(name="Cell Conditions")}


alternatPCAPLot<-function(){
catagories<-t(read.table("~/Dropbox/UTX-Alex/jan/catagories"))


swCat<-c("Abnormal","Abnormal","Abnormal","Abnormal","Abnormal","Abnormal","Abnormal","Abnormal","Abnormal","Abnormal","Abnormal","Abnormal","ECFC","MEKA","Stem","Stem","Stem","Normal","Normal","Normal","Normal","Normal")
names(swCat)<-catagories

cats<-as.character(unlist(lapply(catagories,function(x) swCat[x])))

values<-loadHeightFile("~/Dropbox/UTX-Alex/jan/combined_heights.bed")
data<-values$data
info<-values$info

width<-dim(data)[2]
#normData<-standPCAPrep(data)
pcs<-as.matrix(mulcal::ascore(data,c(1,3)))
#plotPCs(pcs,c(1,2),normData,catagories)

pc1<-pcs[,1]
pc2<-pcs[,2]

x<-t(as.matrix(cbind(pc1,pc2))) %*% as.matrix(normData)
dt<-data.frame(x=x[1,],y=x[2,],catagories=as.character(t(cats)),vals=as.character(t(catagories)))

colour=factor(vals)


vals<-c(Abnormal=21,Normal=22,MEKA=23,Stem=24,ECFC=25)
shapelist<-unlist(lapply(cats[order(catagories)],function(x) vals[x]))

ggplot(dt,aes(x,y))+
    geom_point(aes(shape=catagories,colour=vals,fill=vals),size=2)+
    geom_dl(aes(label=catagories),method=list(smart.grid))+
    scale_shape_manual(values=c(Abnormal=21,Normal=22,MEKA=23,Stem=24,ECFC=25),guide="none")+
     theme(panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.border=element_rect(fill=NA))+
guides(shape=FALSE,colour=guide_legend(override.aes=list(shape=shapelist,size=6)))

ggsave("test.pdf",width=1024,height=780)

}

