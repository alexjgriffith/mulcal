#!/usr/bin/env R
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
#
# Usage:
#
# root_file<-"~/Dropbox/UTX-Alex/october/all-data/"
# mock<-cbind("combined")

# for( i in c(1)){
# values<-loadData( paste(root_file,mock[i],"_heights.bed",sep="" ))
# cats<-t(read.table(paste(root_file,"catagories",sep="" )))
# data<-values$data
# normData<-standPCAPrep(data,"colQn")
# pcs<-prcomp(t(normData))
# plotPCs(pcs,cbind(1,2),normData,cats,c( "Abnormal <-> Normal","Stem <-> Differentiated",paste(mock[i]," PC1 Vs PC3", sep="")))}

# library(ape)
# temp<-cor(data)
# rownames(temp)<-cats
# colnames(temp)<-cats
## writes corelations to a table
# write.table(temp,"test.table",quote=FALSE,sep="\t")
# Generates and plots a dendogram
# plot(as.phylo(hclust(dist(temp),method="average")),type="fan")
# title(main="Combinded Mock Dendogram")

# n=0
# for(i in sequence(length(cats))){
#    plotBox(pcs,i,normData,cats)
#    if(n>5){break}
#    n=n+1}

# x<-(t(as.matrix(pcs$rotation)) %*% as.matrix(normData))[c(1,3),]
# cluster<-kmeans(t(x),3,algorithm="Lloyd",iter.max=10)
# lapply(seq(3), function(x) cats[cluster$cluster==x])

ascore<-function(pca,lis){
interlace<-function(...){t(cbind(...))}
lines<-apply(pcs$rotation,2,function(x){c(x[c(1,2,3)],colnames(pcs$rotation)[c(1,2,3)])})
}

loadData<-function(file="~/masters/normal-abnormal/single_heights.bed"){  
  cdata<-read.table(file)
  stats<-cdata[c(1,2,3)]
  l=length(cdata)
  list(stats=as.matrix(stats),
            data=as.matrix(apply(cdata[4:l],2, function(x) as.vector((unlist(x/(stats[3]-stats[2])))))))}

qn <-function(data){
shape<-dim(data)
sequence<-apply(data,2,order)
reverseSequence<-unlist(apply(sequence,2,order))
ranks<-apply(matrix(unlist(lapply(seq(shape[2]),function(i,x,y) x[y[,i],i],data,sequence)),ncol=shape[2]),1,sum)/shape[2]
apply(reverseSequence,2,function(x) ranks[x])
}

standPCAPrep <-function(data,v){
  switch(v,
         rowSumOne=t(apply(data,1, function(x) x/sum(x))),
         colSumOne=apply(data,2, function(x) x/sum(x)),
         row=t(apply(data,1, function(x) (x-mean(x))/var(x))),
         col=apply(data,2, function(x) (x-mean(x))/var(x)),
         rows1=t(apply(data,1, function(x) x/var(x))),
         cols1=apply(data,2, function(x) x/var(x)),
         colQn=qn(data),
         non=data)}

altPCA<-function(data){
  prcomp(1-cor(data))}

plotScatter<-function(pca,r,cats){
  x<-as.vector(pca$rotation[,r[1]])
  y<-as.vector(pca$rotation[,r[2]])
  plot(x,y)
  text(x,y, labels=cats,cex=0.7,pos=3)}

plotHist<-function(pcs,pos){
  x<-pcs$rotation[,pos] * pcs$sdev[pos]
  hist(x,1000)}

plotBox<-function(pcs,pos,data,cats){
  x<-t(as.matrix(pcs$rotation)) %*% as.matrix(data)
  d<-data.frame(x[pos,],row.names=cats)
  boxplot(t(d),las=2)
  text(seq(length(cats)),x[pos,],labels=cats,cex=0.7,pos=3)}

plotPCs<-function(pcs,pos,data,cats,lab=c("xlabel","ylable","Title")){
  x<-t(as.matrix(pcs$rotation)) %*% as.matrix(data)
  d1<-data.frame(x[pos[1],])
  d2<-data.frame(x[pos[2],])
  plot(t(d1),t(d2),xlab=lab[1],ylab=lab[2])
  title(main=lab[3])
  text(x[pos[1],],x[pos[2],],labels=cats,cex=0.7,pos=3)}


interClass<-function(pca,data,pos,C=3,r=13){
  sm<-data %*%  pca$rotation[,pos]
  m<-mean(sm)
  v<-var(sm)
  S<-which(unlist(lapply(sm,function(x) x< m-v*C)))
  G<-which(unlist(lapply(sm,function(x) x> m+v*C)))
  data.frame(data=unlist(lapply(seq(r),function(x) mean(data[G,x])/mean(data[S,x]))),row.names=cats)}

pcaAnalysisTest<- function()
{
root_file<-"~/Dropbox/UTX-Alex/jan/"
mock<-cbind("combined")
i=1

values<-loadData( paste(root_file,mock[i],"_heights.bed",sep="" ))
cats<-t(read.table(paste(root_file,"catagories",sep="" )))
data<-values$data
normData<-standPCAPrep(data,"colQn")
pcs<-prcomp(t(normData))
pc1=3
pc2=6    
plotPCs(pcs,cbind(pc1,pc2),normData,cats,c( paste("PC:",as.character(pc1),sep="" ),paste("PC:",as.character(pc2),sep="" ),mock[i]))

library(ape)
temp<-cor(data)
rownames(temp)<-cats
colnames(temp)<-cats
png("~/Dropbox/dend_corr_average_combiend.png")
dend=hclust(dist(temp),method="average")
plot(dend,hang=-1,main="Combinded Mock Dendogram",xlab="Cell Type")
dev.off()
#title(main="Combinded Mock Dendogram")

png("~/Dropbox/heatmap.png")
heatmap(1-log(cor(values$data)) ,hclustfun=function(x)hclust(x,method="complete") ,distfun=function(x)as.dist(x),scale="column")
dev.off()

}

