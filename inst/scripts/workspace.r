library(mulcal)
library(grid)
library(Biostrings)
library(ggplot2)
library(data.table)


######################
date<-"May 12 2014"
author<-"Alexander Griffith"

plus<-function(...){paste(...,collapse = "",sep="")}

FastaFile<-"~/Dropbox/UTX-Alex/jan/combined.fasta"
v<-loadHeightFile("~/Dropbox/UTX-Alex/jan/combined_heights.bed")
data<-v$data
stats<-v$stats
formatOut<-apply(sub("   ","",stats[,1:3]), 1,paste,collapse="\t")

reg<-mapply(function(pc,loc)buildRegions(data,pc,loc)[,1],
            list(1,1,3,c(3,5),c(3,7),c(3,7),7),
            list("top","bottom","top",c("top","top"),c("top","bottom"),c("top","bottom"),"top"))


Sequences <- readDNAStringSet(FastaFile, "fasta")

folder="~/Dropbox/UTX-Alex/analysis/workspace/"
motifFile<-sapply(c("tall_union_hspc_vs_not_tall_union_hspc.pwm",
             "tall_inter_hspc_vs_not_tall_union_hspc.pwm",
             "tall_vs_not_tall_union_hspc.pwm",
             "hspc_vs_not_tall_union_hspc.pwm"), function(x) paste(folder,x ,sep=""))


regJ<-cbind((reg[,2] | reg[,5]),reg[,2] & reg[,5] , reg[,2],reg[,5])
#regN<-apply(regJ, 1 ,"!")
regN<-(! regJ[,1])

for (i in seq(4))
    motifs<-homerWrapper(Sequences,regJ[,i],regN,"~/Masters/mulcal/inst/lib/homer-4.7/cpp/homer2",motifFile[i],opts=" -S 25 -len 6")

# 1+   Eryt/K562  1-
# 1-   T-All      1+
# 3+   ECFC       3-
# 3+5+ Other      3-5-, 3+5- , 3-5+
# 3+7- HSPC       3-7-, 3+,7+, 3-,7+
# 3+7+ MEKA       3-,7-, 3-,7+ 3+,7-
# 7+   Eryt/Meka  7-

#a<-motifs2View("AGNCAGAC","CANNTG",reg[,2],Sequences,nearHeights=TRUE)
quickplot(a)
length(a)
#b<-motifs2View("AGNCAGAC","CANNTG",reg[,2],Sequences,nearHeights=FALSE)
quickplot(b)
length(b)
#560
#1301

#e<-motifs2View("AGNCAGAC","CANNTG",reg[,5],Sequences,nearHeights=TRUE)
quickplot(e)
length(e)
#f<-motifs2View("AGNCAGAC","CANNTG",reg[,5],Sequences,nearHeights=FALSE)
quickplot(f)
length(f)
# 163
# 305
# -6 AGNCAGAC-CANNTG AGNCAGACANNTG
#c<-motifs2View("AGNCAGAC","CANNTG",reg[,1],Sequences,nearHeights=TRUE)
quickplot(c)
length(c)

#d<-motifs2View("AGNCAGAC","CANNTG",reg[,1],Sequences,nearHeights=FALSE)
quickplot(d)
length(d)
# 748
# 612


runA<-motifs2View("ACCACA","CANNTG",reg[,2],Sequences,nearHeights=FALSE)
quickplot(runA)
length(runA)

runB<-motifs2View("ACCACA","CANNTG",reg[,2],Sequences,nearHeights=TRUE)
quickplot(runB)
length(runB)


maxLoc<-function(x) {which.max(getHeights(x))+min(x)}
quickplot<-function(a,...){plot(do.call(seq,as.list(range(a)+1)),getHeights(a),...); maxLoc(a)}


m12<-c("ACCACA", "CANNTG","AGNCAGAC")
mList<-unlist(lapply(m12,IUPACtoBase))
cList<-lapply(unlist(lapply(m12,IUPACtoBase)),compliment)
locationsM<-lapply(mList,grep,Sequences)
locationsC<-lapply(cList,grep,Sequences)

states<-function(reg,...){
    do.call(rbind,
            lapply(reg,
                   function(x,y)
                       {do.call(rbind,lapply(y, append,x))}
                  ,list(...))
            )
}

st<-states(c(1,2),c(1,2),c(3,2))

h<-apply(st,1,function(x)motifHist(Sequences,mList,cList,locationsM,locationsC,x[1],x[2],reg[,x[3]]))

#histVisualize(h[[2]],m12[1],m12[2])

histVisualize(h,m1,m2)
    h<-unlist(h)
    
    else{
        h<-nearSummit(Sequences,mList,cList,locationsM,locationsC,1,reg)}
 
h
######################



smad<-rep(FALSE,length(Sequences))
smad[grep(IUPACtoBase("AGNCAGAC"),Sequences)]=TRUE
smad[compliment(grep(IUPACtoBase("AGNCAGAC"))),Sequences)]=TRUE

write.table(formatOut[smad&reg[,2]],"tall-smad_peaks.bed",quote=FALSE,row.names=FALSE,col.names=FALSE)


######################

findLocations<-function(ml,Sequences,subset=seq(length(Sequences))){    
        subset[grep(ml,Sequences[subset])]
}

combineCombinations<-function(m1,m2,Sequences,reg){
    l1<-findLocations(m1,Sequences,which(reg))
    findLocations(m2,Sequences,l1)}

singleCombinations<-function(m1,Sequences,reg)
    findLocations(m1,Sequences,which(reg))


complimentLocations<-function(m1,m2,Sequences,reg){
    ml<-unlist(lapply(c(m1,m2),IUPACtoBase))
    cl<-unlist(lapply(ml,compliment))
    l2<-combineCombinations(ml[1],ml[2],Sequences,reg)
    l4<-combineCombinations(cl[1],cl[2],Sequences,reg)
    sort(union(l2,l4))
}

getMembers<-function(m1,Sequences){
        mL1<-IUPACtoBase(m1)
        cL1<-paste(mL1,compliment(mL1),sep="|")
        gregexpr(cL1,Sequences)}

getDistanceWrap<-function(m1,m2,Sequences){
    getDistance(getMembers(m1,Sequences),getMembers(m2,Sequences))
}

slowDist<-function(m1,m2,Sequences,reg){
    l5<-complimentLocations(m1,m2,Sequences,reg)
    getDistanceWrap(m1,m2,Sequences[l5])}

slowIndv<-function(m1,Sequences,reg){
    l2<-sort(union(singleCombinations(IUPACtoBase(m1),Sequences,reg),singleCombinations(compliment(IUPACtoBase(m1)),Sequences,reg)))
    mL1<-IUPACtoBase(m1)
    cL1<-paste(mL1,compliment(mL1),sep="|")
    unlist(lapply(gregexpr(cL1,Sequences[l2]),as.numeric))}


jointDist<-function(m1,m2,Sequences,reg){
    lM<-complimentLocations(m1,m2,Sequences,reg)
    tbM<-list(getMembers(m1,Sequences[lM]),getMembers(m2,Sequences[lM]))
    val1<-do.call(rbind,Map(function(x,y){
        x<-unlist(x)
        y<-unlist(y)
        c<-length(x)*length(y)
        a<-c/length(x)
        b<-c/length(y)
        cbind(unlist(lapply(x,function(x)rep(x,a))),rep(y,b))},
                            tbM[[1]],tbM[[2]]))
}

rotation<-function(val1){
    val2<-val1%*%matrix(c(-1,1,0,0),ncol=2)
    val3<-val1%*%matrix(c(1,1,0,0),ncol=2)   
    cbind(val2[,1],val3[,1])}

projectDiff<-function(val1P){
    getHeights(val1P[,1],c(-300,300))-getHeights(val1P[,2],c(0,600))
}

windowProjectDiff<-function(y,m1,m2,fun,window=10){
    #a<-floorWidth(getHeights(y[,1],c(-300,300)),seq(-300,300),-nchar(m1),nchar(m2))
    a<-getHeights(y[,1],c(-300,300))
    a1<-a[(301+nchar(m2)):600]
    a2<-a[1:(300-nchar(m1))]
    b<-getHeights(y[,2],c(0,600))
    b1<-b[(301+nchar(m2)):600]
    b2<-b[1:(300-nchar(m1))]
    aP<-c(rev(windowing(rev(a2),fun,window,sliding=TRUE)),
          windowing(a1,fun,window,sliding=TRUE))
    bP<-c(rev(windowing(rev(b2),fun,window,sliding=TRUE)),
          windowing(b1,fun,window,sliding=TRUE))    
    cbind(bP,aP)
}

floorFunction<-function(fun,m1,m2,Sequences,reg,pass=FALSE){
    y<-fun(m1,m2,Sequences,reg)
    if(! pass)
        floorWidth( y,seq(-300,300),-nchar(m1),nchar(m2))
    else
        y
}



lapply(mots[1:10], function(x){cor(windowProjectDiff(rotation(jointDist(x,m2,Sequences,reg[,1])),m1,m2))[[2]]^2})

m2<-mots[301]

for(x in mots[1:20]){
    y<-rotation(jointDist(x,m2,Sequences,reg[,1]))
    v<-windowProjectDiff(y,m1,m2,fun="var",window=5)
    fit<-lm(y~x,data.frame(y=v[,2],x=v[,1]))
    c<-coefficients(fit)
    x=seq(length(v[,1]))/length(v[,1])*max(v[,1])
    l<-data.frame(x=x,y=x*c[2]+c[1])
    par(mfrow=c(2,2))
    plot(v)
    lines(l)    
    #plot(v[,2]/l$y)
    #plot((v[,2]-l$y[order(v[,1])])/l$y[order(v[,1])])
    #plot((v[,2][order(v[,1])]/l$y))
    plot(unlist(lapply(v[,2],function(x)max(0.01,x)))/unlist(lapply(v[,1],function(x)max(0.01,x))))
    plot(seq(-300,300),floorWidth(getHeights(y[,1],c(-300,300)),seq(-300,300),-nchar(m1),nchar(m2)) )
    plot(seq(0,600),getHeights(y[,2],c(0,600)))
    readKey()
}

fit<-lm(y~x,data.frame(y=v[,2],x=v[,1]))


plot(getHeights(jointDist(m1,m2,Sequences,reg[,1])%*%matrix(c(0,1,0,1),ncol=2)[,1] ) )

plot(floorFunction( function(...)projectDiff(rotation(jointDist(...))),m1,m2,Sequences,reg[,1]))


y<-slowDist(m1,m2,Sequences,reg[,1])
a<-floorWidth( getHeights(y,c(-300,300)),seq(-300,300),-nchar(m1),nchar(m2))
plot(a)


rev(order(getHeights(a))+min(a))[1:10]

m1<-"GATAAG"
#m1<-"KNWVNAN"
m2<-"CANNTG"

m1<-mots[runif(1,1,302)]
m2<-mots[runif(1,1,302)]

y<-projectDiff(m2,m1,Sequences,reg[,1])
plot(y[,1],y[,2])
a<-floorWidth( getHeights(y[,1],c(-300,300))-getHeights(y[,2],c(1,601)),seq(-300,300),-nchar(m1),nchar(m2))
var<-distributionT(data.frame(x=seq(1:300),y=rep(1/601,300)))*var(a)*601
#plot(-log10(dnorm(a,mean(a),sqrt(var))))



plot(windowing(a,function(x) var(x),50))
lines(var[25:(length(var)-25)])

plot(windowing(a,function(x) var(x),50)-var[25:(length(var)-25)])

m<-mean(a)



plot(mapply(function(x,v)-log10(dnorm(a,mean(a),var(a))),var,a))

mean(a)


plot(windowing(floorFunction(function(...) getHeights(slowDist(...),c(-300,300)),"AGTANBA","AAHATN",Sequences,reg[,1]),function(x) -log10(dnorm(x[25],mean(x),var(x))),50 )[150:400],ylab="")

plot(getHeights(slowDist("AGTANBA","AAHATN",Sequences,reg[,1])))


#mots<-c(mots[1:(n-length(addmotifs))],addmotifs)
mL<-unlist(lapply(m1,IUPACtoBase))
cL<-unlist(lapply(lapply(m1,IUPACtoBase),compliment))
locationsM<-lapply(mL,grep,Sequences)
locationsC<-lapply(cL,grep,Sequences)

tm<-IUPACtoBase(m2)
loc<-which(reg[,1])
sub<-intersect(loc,grep(tm,Sequences))
a<-unlist(lapply(gregexpr(tm, Sequences[sub]),as.numeric))
plot(getHeights(a))


t0f<-loadFlatData("motifMatrix.txt")
rownames<-rownames(t0f)
colnames<-colnames(t0f)
mots<-unique(unlist(lapply(rownames,function(x) strsplit(x,"-")[[1]][1])))





m<-motifs2View("CANNTG","TGACCT",reg[,3],Sequences)

m1<- "RGACAT"
m2<- "AAHATN"
c=1

m1<- "NTNSGVSYNG"
m2<- "GMTGG"
c=1

c=2
m1<-"MAGDS"
m2<-"KGVNNGNAA"

m1<-"KNWVNAN"
m2<-"HTYKNABNK"

m1<-"TAATCC"
m2<-"YACYWKGN"    
c<-7

m1<-"CANNTG"
m2<-"MAGDS"
c<-1
unlist(Map(function(x) as.character(jasparVals[which(as.character(jasparVals[,1])==x),2]),c(m1,m2)))    
    lables<-paste(m1,m2,c,sep="-")
    y<-as.numeric(t0f[which(rownames==lables),])
    x<-as.numeric(colnames)
    v<-data.frame(x=x,y=y)
v[rev(order(floorWidth(v$y,v$x,-nchar(m1),nchar(m2))))[1:10],]
plot(v$x,floorWidth(v$y,v$x,-nchar(m1),nchar(m2)))

m<-singleScore(m1,m2,5,function(x) -log10(binomTri(x)))
plot(m)





singleScore<-function(m1,m2,c,fun){
    lables<-paste(m1,m2,c,sep="-")
    y<-as.numeric(t0f[which(rownames==lables),])
    x<-as.numeric(colnames)
    v<-data.frame(x=x,y=y)
    score<-do.call(fun,list(v))
    score<-floorWidth(score,v$x,-nchar(m1),nchar(m2))
    score
    }


#plot(x$x[x$x>-50 & x$x<50],x$y[x$x>-50 & x$x<50])

############


results<-testScore(y[lines>4,],lines[lines>4],windowNormalizedFoldChange)

results<-testScore(y[1:100,],lines[1:100],function(x) -log10(binomTri(x)))

cut<-3
TP<-length(which(vals>cut & lines>1))
FP<-length(which(vals>cut & lines<2))
TN<-length(which(vals<cut & lines<2))
FN<-length(which(vals<cut & lines>1))



b<-do.call(plotPort,as.list(y[86,] ))


# 117
# 94
# 97
v<-do.call(quickSelect,as.list(y[97,]))


r<-singleScore(y,lines,122,function(x) -log10(binomTri(x)))

l<-c()
n<-301
m<-1000
y<-t(combn((1:n),2))[sample(seq(dim(combn((1:n),2))[2]))[1:m],]
#y<-cbind(runif(n,1,301),runif(n,1,301))
for (i in seq(m)){
    r<-singleScore(y,c(1),i,function(x) -log10(binomTri(x)))
    #plot(do.call(quickSelect,as.list(y))$y)
    l<-rbind(l,cbind(max(r),mean(r),sqrt(var(r))))}

scores<-(l[,1]-l[,2])/l[,3]

scores<-l[,1]
rev(order(scores,na.last=FALSE))[1:5]

    

## Data Prep
t0f<-loadFlatData("motifMatrix.txt")
#t0f<-read.table("motifMatrix.txt",header=TRUE,check.names=FALSE,comment.char = "", nrows=70000,colClasses=c("character",rep("numeric",601)) )
rownames<-rownames(t0f)
colnames<-colnames(t0f)
mots<-unique(unlist(lapply(rownames,function(x) strsplit(x,"-")[[1]][1])))
bend<-apply(t0f,2,mean)

noMidBend<-mapply(floorWidth,bend,x)


v<-selectValue(t0f,mots[3],mots[10])#,norm=bend)
hist(v$y)
#ggplot(v,aes(x,y))+geom_step()
t0f[1,0]
t0m<-as.matrix(t0f)

x<-as.numeric(colnames(t0m))

#region<-apply(t0mstrsplit(rownames(t0f[1,]),"-")
motifLengths<-unlist(lapply(mots,nchar))
regions<-lapply(motifLengths,function(l1)lapply(motifLengths,function(l2) rep(cbind(l1,l2),7)))


#ports<-cbind(c(1,1),c(2,1),c(1,2),c(2,2),c(1,3),c(2,3))
for (i in seq(length(mots))){
    a<-motifLengths[i]
    for (j in seq(length(mots))){
        #for (k in 6){
        k<-1               
        v<-selectValue(t0f,mots[i],mots[j],loc=k)
        v$y<-floorWidth(v$y,v$x,-nchar(mots[i]),nchar(mots[j]))
        m<-which.max(v$y)
        name<-paste(i,j,mots[i],mots[j],"(",as.character(v$x[m]),"(",sep="-")
        cols<-rep(1,length(v$x))
        n=c(-nchar(mots[i]),nchar(mots[j]))
        w=unlist(lapply(n,function(r)which(v$x == r)))
        cols[w]=3
        n=c((-nchar(mots[i])-1):-25,(nchar(mots[j])+1):25 )
        w=unlist(lapply(n,function(r)which(v$x == r)))
        cols[w]=4
        cols[m]=2
        plot(v$x,v$y,xlab=name,col=cols)
                    Sys.sleep(1)
        #}
        #b<-motifLengths[j]
        #nonMid<-(x>b | x<(-a))
        #for (k in seq(7)){
        #    sum(y[nonMid])})
#}
}
}
            



size<-apply(t0m,1,sum)

size<-apply(t0m,1,
            function(y){                                
                a<-10
                b<-10
                nonMid<-(x>b | x<(-a))
                sum(y[nonMid])})

dense<-size>1200

for (i in which(dense)[1:4])
    plot(as.numeric(t0f[i,]))

v<-selectValue(t0f,mots[3],mots[10],range="top",norm=bend)
qplot(v$x,as.numeric(cdf(v$y)))+xlab("Post Normalization")+ylab("Events")

selectValue<-function(t0f,m1,m2,loc=1,
                      name=paste(m1,m2,loc,sep="-"),
                      tip=t0f[name,],
                      range="none",norm=rep(1,length(tip))){
    x<-as.numeric(colnames(t0f))
    r<-switch(range,
              top=(x>=0),
              bottom=(x<=0),
              mid=mapply(all,x>=-100,x<=100),
              none=rep(TRUE,length(tip)))
    data.frame(y=as.numeric(tip)[r]/norm[r],x=x[r])
}


visNorm(mots,t0f,mots[1],mots[10])

#' @export
visNorm<-function(mots,
                  t0f,
                  m1="TACCTC",
                  m2="TAKNNMABKN",
                  loc=1,
                  name=paste(m1,m2,loc,sep="-"),
                  tip=t0f[name,]){
    x<-as.numeric(colnames(t0f))
    top<-x>=0
    bottom<-x<=0
    midReg<-mapply(all,x>=-100,x<=100)
    if(! all(sapply(c(m1,m2),function(x) x %in% mots )))
        stop()
    y<-as.numeric(tip)
    y2<-y/bend
    #png("heights_comp.png")
    pushViewport(viewport(layout=grid.layout(3,2)))
    p1<-qplot(y[outsideWidth(x[midReg])])+xlab(name)
    p2<-qplot(y2[outsideWidth(x[midReg])])+xlab(name)
    p3<-locHist2(y[midReg],seq(-100,100),limits=c(-100,100))
    p4<-locHist2(y2[midReg],seq(-100,100),limits=c(-100,100))
    p5<-qplot(x[top],as.numeric(cdf(y[top])))+xlab("Pre Normalization")+ylab("Events")
    p6<-qplot(x[top],as.numeric(cdf(y2[top])))+xlab("Post Normalization")+ylab("Events")
    print(p1,vp=viewport(layout.pos.row=1,layout.pos.col=1))
    print(p2,vp=viewport(layout.pos.row=1,layout.pos.col=2))
    print(p3,vp=viewport(layout.pos.row=2,layout.pos.col=1))
    print(p4,vp=viewport(layout.pos.row=2,layout.pos.col=2))
    print(p5,vp=viewport(layout.pos.row=3,layout.pos.col=1))
    print(p6,vp=viewport(layout.pos.row=3,layout.pos.col=2))
    #dev.off()
}






partialSigMatrix(Sequences,mots,reg,n=c(4,4))
Rprof("speed.test")
#saveData("test",mode="all",append=FALSE,h,seq(1,7),mots[1:10],range(unlist(h)))
saveData("test",mode="perMot",append=FALSE,col.names=TRUE,row.names=TRUE,h,mots[1:10],mots[1:10],seq(1,7),c(-300,300))
Rprof(NULL)
summaryRprof("speed.test")


    len<-m*n[1]*n[2]
    bp<-1000
    
    r<-modulous(len,bp)
    d<-floor(len/bp)
    append(lapply(seq(d),function(x)seq(((x-1)*bp+1),x*bp)),list(seq((d*bp+1),(d*bp+r))))
    combinations()

n<-120
a<-10
b<-6
c<-2

list(value=101,a=1,b=4,c=2)
list(value=49,a=9,b=4,c=1)

testFun1<-function(x){
    floorDiff<-function(x,m){floor(x/m)}
    c(a=(modulous(x,a)+1),modulous(floorDiff(x,a),b),floorDiff(floorDiff(x,a),b ))
}


# 1+   Eryt/K562  1-
# 1-   T-All      1+
# 3+   ECFC       3-
# 3+5+ Other      3-5-, 3+5- , 3-5+
# 3+7- HSPC       3-7-, 3+,7+, 3-,7+
# 3+7+ MEKA       3-,7-, 3-,7+ 3+,7-
# 7+   Eryt/Meka  7-

l<-mapply(function(a,b,c){
    motifsFromPCA(Sequences,data,a,b,title=c)
},
list(1,1,3,c(3,5),c(3,7),c(3,7),7),
list("top","bottom","top",c("top","top"),c("top","bottom"),c("top","bottom"),"top"),
paste("~/Masters/mulcal/inst/data/",list("normal","abnormal","ecfc","other","hspc","meka","diff"),".pwm",sep=""))

l<-motifsFromPCA(Sequences ,data,c(1,2,3),c("top","top","bottom"))


m<-motifs2View("GATAAG","CATTCT",reg[,1],Sequences,TRUE)

m<-motifs2View("GATAAG","CATTCT",reg[,1],Sequences)

mn1<-motifs2View("GATAAG","CANNTG",reg[,1],Sequences)
mn2<-motifs2View("GATAAG","CCTAAT",reg[,1],Sequences)

m1<-mn1[mn1>0]
x<-cdf(getHeights(m1,seq(0,max(m1))))
start<-which(x>0)[1]
y1<-mean(x[start:start+20])
x1<-((start+20)+start)/2
out<-c()
for(i in (start+21):length(x)){
    x2<-i
    y2<-x[i]
    m<-(y2-y1)/(x2-x1)
    b<-y1-x1*m
    vim<-sum(unlist(lapply((start+20):i,function(i) x[i]- m*i+b)))/i^2
    out<-c(out,vim)
}
plot(out)

out
#find the linear range of x+2%

n<-10
y<-unlist(lapply((start+1):(length(x)-n),function(i){cor(i:(i+n),x[i:(i+n)])}))
z<-pdf(c(rep(0,start),pdf(y)))
plot(z)

a<-which(y>sqrt(min(unlist(Map(function(x)var(c(y[x:(x+20)])),c(1,20,40,60)))))*10)

plot(pdf(c(rep(0,start),y)))

#a[a>20]

uc<-92+6
lc<-111+6
m2<-m[apply(cbind(m>-lc,m<uc),1,all)]

png("lin-_k563Eryt_gattag_canntg.png")
plot(y)
dev.off()







x<-do.call(seq,as.list(range(unlist(h))))



       
       scores<-motifGetScore(h,1,mots,addmotifs)

plot(hclust(as.dist(scores[[1]]),method="average"),hang=-1)

quick2motifs("AACCAC","ACCGGM",mots[1:17])


rev(order(unlist(lapply(seq(1:800),function(x) scoreSystem1(data[x,])[1]))))[1:10]

n<-48      
locHist2(as.numeric(data[n,]),as.numeric(colnames(data[n,])),rownames(data[n,]),limits=c(-200,200))
       
colnames(data[n,])[which.max(data[n,])]
       
       
scoreSystem1<-function(t0){
    data<-as.numeric(t0)
    x<-as.numeric(colnames(t0))
    dets<-breakDownName(t0)
    width<-c(-nchar(dets[2]),nchar(dets[1]))
    zData<-zeroWithinWidth(data,width,x)
    prob<-1-mean(zData)/sqrt(var(zData))
    r<-(1-prob)/prob*mean(zData)
    #c(sum(zData),max(zData),length(zData),mean(zData),var(zData))
    #sum(zData)
    list(prob=prob,r=r,l=sum(zData))
}
       
n<-10
       
       d<-scoreSystem1(data[n,])
       
       heightHist(as.numeric(data[n,]))
       
       locHist2(as.numeric(data[n,]),as.numeric(colnames(data[n,])),"",c(-100,100))
       
       with(d,{
       plot(0:50,unlist(lapply(0:50,dnbinom,size=r,prob=(prob)))*l )})
       
lapply(a,function(x)paste(strsplit(as.character(x),"")[[1]][1:2],collapse="")=="#$")
       

motifGetScore<-function(h,m,mots,addmotifs,upt=TRUE){
    n<-length(h)
    mots<-c(mots[1:(n-length(addmotifs))],addmotifs)
    if(upt)
        scores<-lapply(seq(m),
                       function(m){
                           v<-apply(combn(n,2),2,
                                    function(x){
                                        v<-h[[x[1]]][[x[2]]][[m]]
                                        t1<-t0<-getHeights(v)
                                        if(length(t0)>0)
                                        while(abs(order(t0,decreasing = TRUE)[1]+min(v)) > max(length(mots[x[1]]),length(mots[x[2]]) )){
                                            
                                            t0[order(t0)]<-0
                                        }
                                        
                                        motifScoreFunction(v,t1,max(t0))})
                           y<-matrix(0,n,n)
                           y[upper.tri(y)]<-v
                           y[lower.tri(y)]<-v
                           colnames(y)<-mots
                           rownames(y)<-mots
                           y
                       })
    else
        scores<-lapply(seq(m),
                       function(n){
                           ma<-matrix(unlist(lapply(h,function(x) lapply(x,function(y)motifScoreFunction(y[[n]]) ))),nrow=10)
                           colnames(ma)<-mots
                           rownames(ma)<-mots
                           ma})
    scores
}






