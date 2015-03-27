
    

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




