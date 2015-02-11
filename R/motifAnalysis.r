#!/usr/bin/env R
#
# This file is part of peakAnalysis,
# http://github.com/alexjgriffith/alpha-score/, 
# and is Copyright (C) University of Ottawa, 2015. It is Licensed under 
# the three-clause BSD License; see LICENSE.txt.
# Author : Alexander Griffith
# Contact: griffitaj@gmail.com
#

#library(Biostrings)

#' Turn IUPAC lables into basic neucleotides
#'
#' @param char A string of IUPAC characters
#' c("A","G","C","T","R","Y","S","W","K","M","B","D","H","V","N")
#' @param rl choice to return a value
#' \describe{
#' \item{\strong{TRUE}}{Each posible string is genearated}
#' \item{\strong{FALSE}}{Direct translation to brakets R=[AG]}}
#' @export
#' @template authorTemplate
IUPACtoBase<-function(char,rl=FALSE){
    IUPAC<-strsplit(char,"")[[1]]
    IUPACCharacters<-list("A","G","C","T",c("A","G"),c("C","T"),c("G","C"),c("A","T"),c("G","T"),c("A","C"),c("C","G","T"),c("A","G","T"),c("A","C","T"),c("A","C","G"),c("A","C","G","T"))
    names(IUPACCharacters)<-c("A","G","C","T","R","Y","S","W","K","M","B","D","H","V","N")

    if(any( ! IUPAC %in% names(IUPACCharacters)))
        {
            stop("Input string contains non IUPAC characters.")
        }
    vals<-sapply(IUPAC,function(x){(IUPACCharacters[x])})
    if(rl){
        Base<-c("")
        for(i in vals){
            kit<-c()
            for(j in Base)
                kit<-c(kit,paste(j,i,sep=""))
            Base<-kit}}
    else{
        Base<-c("")
        for(i in vals)
            if(length(i)==1)
                Base<-paste(Base,i,sep="")
            else
                Base<-paste(Base,"[",do.call(paste, as.list(c(i,sep=""))),"]",sep="")
    }
    Base}

#' reverse compliment
#'
#' Takes a string which can include base nucleaotides (AGCT) and brakets ("[","]") and returns the reverse compliment including brakets
#' 
#' @export
#' @template authorTemplate
compliment<-function(string){
    chars<-c("A","G","C","T","[","]")
    names(chars)<-c("T","C","G","A","]","[")
    paste(rev(sapply(strsplit(string,"")[[1]], function(x){(chars[x])})),collapse="")
}

#' Motif Testing
#' @export
motifTest<-function(fastaFile="~/Dropbox/UTX-Alex/jan/combined.fasta",
                    heightFile="~/Dropbox/UTX-Alex/jan/combined_heights.bed",
                    addmotifs=c("CGNNGC")){
    values<-loadHeightFile(heightFile)
    data<-values$data
    reg4<-ascore(data,1,"top",3)
    reg4<-rep(TRUE,length(data[,1]))
    test<-readDNAStringSet(fastaFile,use.names = TRUE)
    motifs<-as.matrix(read.table("data/normal_not_abnormal_motifs"))
    mList<-unlist(lapply(c(motifs,addmotifs),IUPACtoBase))
    cList<-unlist(lapply(lapply(c(motifs,addmotifs),IUPACtoBase),compliment))
    locationsM<-lapply(mList,grep,test)
    locationsC<-lapply(cList,grep,test)
    l<-length(cList)
    h<-motifHist(test,mList,cList,locationsM,locationsC,24,25,reg4)
    h1<-unlist(lapply(seq(from=-300,to=300),function(x){length(which(x==h))}))
               hist(h1[! h1==0],breaks=max(h1))
       }
    
#' Motif histogram
#'
#' @examples
#' 
#' values<-loadHeightFile(heightFile)
#' data<-values$data
#' reg<-ascore(data,1,"top",3)
#' test<-readDNAStringSet(fastaFile,use.names = TRUE)
#' motifs<-as.matrix(read.table("data/normal_not_abnormal_motifs"))
#' mList<-unlist(lapply(c(motifs,addmotifs),IUPACtoBase))
#' cList<-unlist(lapply(lapply(c(motifs,"CGNNGC"),IUPACtoBase),compliment))
#' locationsM<-lapply(mList,grep,test)
#' locationsC<-lapply(cList,grep,test)
#' l<-length(cList)
#' motifHist(mList,cList,locationsM,locationsC,4,l,reg)
#' 
#' @export
 motifHist<-function(data,mList,cList,locationsM,locationsC,n1,n2,reg){
    lM<-intersect(intersect(locationsM[[n1]],locationsM[[n2]]),which(reg))
    lC<-intersect(intersect(locationsC[[n1]],locationsC[[n2]]),which(reg))
    bM<-c()
    bC<-c()
    if(length(lM)>0)
        bM<-lapply(c(mList[n1],mList[n2]),function(x) lapply(gregexpr(x, data[lM]),as.numeric))
    if(length(lC)>0)        
        bC<-lapply(c(cList[n1],cList[n2]),function(x) lapply(gregexpr(x, data[lC]),as.numeric))
    h<-c(getDistance(bM[[1]],bM[[2]]),getDistance(bC[[1]],bC[[2]]))
    hist(h,breaks=600,xlim=range(-30,30),xlab=paste(mList[n1],mList[n2],sep=" - "))
    h}




#' @export
getDistance<-function(x,y){
    as.numeric(
        mapply(function(x,y){
        temp<-outer(x, y,"-")
        temp<-temp[upper.tri(x=temp,diag=TRUE)]
        temp[which.min(abs(temp))]
        #temp
    },x,y))}



