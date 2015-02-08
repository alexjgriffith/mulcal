#!/usr/bin/env R
#
# This file is part of peakAnalysis,
# http://github.com/alexjgriffith/alpha-score/, 
# and is Copyright (C) University of Ottawa, 2015. It is Licensed under 
# the three-clause BSD License; see LICENSE.txt.
# Author :Alexander Griffith
# Contact: griffitaj@gmail.com
#
#
######################################################################
######################################################################

strList<-"(intersection jurk cem rpmi (not (union k562 eryt)))"

subList<-list("union","jurk","cem","rpmi",list("not",list("union","k562","eryt")))

for(i in subList[2:length(subList)]){if(is.list(i)){print i}}

website<-htmlDoc(htmlTags("body",collapse(htmlTags("H1","Welcome"),
         htmlTags("p","This is a quick sample pagargaph"),
         htmlTable(matrix(c(1,2,3,4),2,2)))))
write(website,"test.html")

nb<-buildAnotations("cellpadding", "0")
cat(htmlTable(matrix(c(htmlTable(matrix(c(1,2,3),3,1),nb),htmlTable(matrix(c(1,2,3),3,1),nb),htmlTable(matrix(c(1,2,3),3,1),nb),htmlTable(matrix(c(1,2,3),3,1),nb)),2,2)))

matrix(list(list(1,2,3),list(1,2),list(1),list(1,2,3,4)),ncol=2)[1,1]

#allFasta<-readDNAStringSet(paste(fileLocation,"combined.fasta",sep=""),use.names=TRUE)
#foreg<-allFasta[member(strsplit(as.character(types[,4]),"-"),c("k562","eryt"))]
#backg<-allFasta[member(strsplit(as.character(types[,4]),"-"),c("jurk","cem","rpmi"))]



interlace<-function(...){
    ma<-max(sapply(list(...),length))
    m<-mapply( function(x,y){c(y,rep("",x))},ma-sapply(list(...),length),list(...))
    unlist(apply(m,1,function(x){list(x)}))}


wraps<-function(c,val="'"){
    apply(as.matrix(as.character(c)),1,function(x){paste(val,x,val,sep="")})}

unwrap<-function(string,n=1){
    a<-strsplit(string,"")
    lcollapse(a[[1]][(n+1):(length(a[[1]])-n)])}

makeOutDict<-function(values,keys,len=length(values)){
    if(len>1){
        lcollapse(c("{",interlace(wraps(keys)[1:len],rep(":",len),wraps(values)[1:len],rep(",",len-1)),"}"))}
    else{
        lcollapse(c("{","'",keys[1],"'",":","'",values[1],"'","}"))}}

readInDict<-function(string){
    dict<-strsplit( unwrap(string),",")[[1]]
    r<-c()
    for (i in dict){
        t<-strsplit(i,":")[[1]]
        v<-as.numeric(unwrap(t[2]))
        names(v)<-t[1]
        r<-c(r,v)}
    r}

vals=c(1,3)
keys<-attributes(pcs$rotation)$dimnames[[2]][vals]
lemp<-apply(pcs$rotation[,vals],1,makeOutDict,keys,2)
as.data.frame(t(apply(matrix(lemp),1,readInDict)))

