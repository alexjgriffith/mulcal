library(mulcal)
library(grid)
library(Biostrings)
library(ggplot2)
library(data.table)


loadFlatData<-function(filename){
    colnames<-as.numeric(strsplit(readLines(filename,n=1)," ")[[1]])
    t0f<-fread(filename,header=FALSE,skip=1)
    f<-as.data.frame(t0f)
    rownames<-f[,1]
    f<-f[,2:dim(f)[2]]
    colnames(f)<-colnames
    rownames(f)<-rownames
    f    
}

floorWidth<-function(x,y,a=6,b=-a,out=0){
    mapply(function(x,y){
        if(y>a & y<b) {0}
        else {x}},x,y) }


singleScore<-function(m1,m2,c,fun){
        v<-quickSelect(m1,m2,c)
        score<-do.call(fun,list(v))
        score<-floorWidth(score,v$x,-nchar(mots[m1]),nchar(mots[m2]))
        score
    }

distributionT<-function(x){
    r<-max(abs(range(x$x)))
    l<-length(x$x)
    d<-seq(r)/sum(seq(r))*sum(x$y)
    c(d,rev(d)[1],rev(d))    
}
binomTri<-function(v){
    width=300
    d<-distributionT(data.frame(x=seq(-width,width),y=1/(width*2+1)))
    binomDist(v,d)
    }
binomDist<-function(v,d){
    n<-sum(v$y)
    unlist(lapply(seq(length(v$x)),function(x)dbinom(v$y[x],n,d[x])))}

quickSelect<-function(x,y,loc=1)
{
    m1<-mots[x]
    m2<-mots[y]
    name=paste(m1,m2,loc,sep="-")
    tip=t0f[name,]
    x<-as.numeric(colnames(t0f))
    data.frame(y=as.numeric(tip),x=x)
}


t0f<-loadFlatData("motifMatrix.txt")
rownames<-rownames(t0f)
colnames<-colnames(t0f)
mots<-unique(unlist(lapply(rownames,function(x) strsplit(x,"-")[[1]][1])))

opt<-t(combn(length(mots),2))
write.table("motifA-motifB-Condition-Score","binomTriScore.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)
for (i in seq(dim(opt)[1])){
    m1<-opt[i,1]
    m2<-opt[i,2]
    score<-singleScore(m1,m2,1,function(x) -log10(binomTri(x)))
    value<-paste(mots[m1],mots[m2],1,max(score),sep="-")
    print(value)
    write.table(value,"binomTriScore.txt",append=TRUE,quote=FALSE,col.names=FALSE,row.names=FALSE)
           }
