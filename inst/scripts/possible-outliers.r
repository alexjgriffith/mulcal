library(mulcal)
library(grid)
library(Biostrings)
library(ggplot2)
library(data.table)


t0f<-loadFlatData("motifMatrix.txt")
rownames<-rownames(t0f)
colnames<-colnames(t0f)
mots<-unique(unlist(lapply(rownames,function(x) strsplit(x,"-")[[1]][1])))

opt<-t(combn(n,2))
results<-c()
for (i in seq(dim(opt)[1]))
    results<-c(results,
               singleScore(opt[1,1],opt[1,2],function(x) -log10(binomTri(x))))

#y<-t(matrix(unlist(lapply(unlist(Map(as.character,read.table("test.csv"))),
#       function(x) {
#           y<-strsplit(x,"-")[[1]]
#           c(as.numeric(y[1]),as.numeric(y[2]))
#       })),nrow=2))

#lines<-as.numeric(read.table("handScores.txt")[,4])

testScore<-function(y,lines,fun){
    n=length(lines)
    scores<-unlist(lapply( seq(n),function(i){
        v<-quickSelect(y[i,1],y[i,2])
        v$y<-floorWidth(v$y,v$x,-nchar(mots[y[i,2]]),nchar(mots[y[i,1]]))
        max(floorWidth(do.call(fun,list(v)),v$x,-nchar(mots[y[i,2]]),nchar(mots[y[i,1]])))
    }))
    data.frame(actual=lines ,observed=scores)}

singleScore<-function(m1,m2,c,fun){
        v<-quickSelect(m1,m2,c)
        score<-do.call(fun,list(v))
        score<-floorWidth(score,v$x,-nchar(mots[m1]),nchar(mots[m2]))
        score
    }

foldChange<-function(x){
    m<-mean(x$y)
    if(m>0)
        max(x$y)/m
    else
        1
}

normalizedFoldChange<-function(x){
    m<-mean(x$y)
    m2<-mean(x$d)
    if (sum(x$y)<length(x$y))
        0
    else if(m>0 & m2>0)
        m/m2
    else
        0    
}


distributionTVar<-function(x,norm=var(x$y))
{
    r<-seq(max(abs(range(x$x))))
    d<-r/length(r)*norm
    c(d,rev(d)[1],rev(d))    
}

    
distributionT<-function(x){
    r<-max(abs(range(x$x)))
    l<-length(x$x)
    d<-seq(r)/sum(seq(r))*sum(x$y)
    c(d,rev(d)[1],rev(d))    
}

windowNormalizedFoldChange<-function(x){
    d<-distributionT(x)
    y<-cbind(x,d)
    max(windowing(y,normalizedFoldChange,10))
    
}

windowFoldChange<-function(x) max(windowing(x,foldChange,100))


binomDist<-function(v,d){
    n<-sum(v$y)
    unlist(lapply(seq(length(v$x)),function(x)dbinom(v$y[x],n,d[x])))}

binomTri<-function(v){
    width=300
    d<-distributionT(data.frame(x=seq(-width,width),y=1/(width*2+1)))
    binomDist(v,d)
    }

meanNorm<-function(window=50){
    aux<- function(x) windowing(x,function(x) mean(unlist(x$y)),window,TRUE)    
    function(x){
        aux(x)
    }    
}
   
visualizeWeightedMeanNorm<-function(v){
    compDist<-data.frame(y=distributionT(v),x=seq(-300,300))
    a<-meanNorm(1)(v)
    b<-meanNorm(1)(compDist)
    plot(a,ylim=c(0,max(c(a,b))))
    lines(b)
    c<-a/b
    which.max(c)
}


quickplot<-function(n,...){
    v<-do.call(quickSelect,as.list(y[n,]))
    motifs<-mots[y[n,]]
    v$y2<-floorWidth(v$y,v$x,-nchar(motifs[2]),nchar(motifs[1]))
    n<-paste("max=",v$x[which.max(v$y2)],sep="")
    plot(v$x,v$y2,xlab=paste(c(motifs,n),collapse="-"),...)
}

