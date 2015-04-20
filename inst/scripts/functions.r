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

outsideWidth<-function(y,a=6,b=-a,out=0){
    out<-c()
    for(x in y)
        if(all(x<a,x>b)) {out<-c(out,FALSE)}
        else {out<-c(out,TRUE)}
    out
}

normTest<-function(v){
    r1<-dnorm(0)
    l1<-dnorm(0)
    if(sum(v$y[v$x>6])>50){
        cutOff<-linearRange(cdf(v$y[v$x>0]))
        if(cutOff>20){
            r2<-unlist(lapply(6:cutOff,function(x) linearity(v$y[v$x>0],c(6,cutOff),x)))
            r1<-dnorm(ztest(r2))}}
    if(sum(v$y[v$x<-6])>50){
        cutOff<-linearRange(cdf(v$y[v$x<0]))
        if(cutOff>20){
            r2<-unlist(lapply(6:cutOff,function(x) linearity(rev(v$y[v$x<0]),c(6,cutOff),x)))
            l1<-dnorm(ztest(r2))}}
    list(r1,l1)
}

zvalue<-function(i,j){
    v<-quickSelect(i,j)
    normTest(v)}

ztest<-function(r2)
    (max(r2)-mean(r2))/sqrt(var(r2))


linearRange<-function(vip,top=270,bottom=20,lit=1){
    lr<-unlist(lapply(lit+bottom:top,function(x)cor(vip[lit:(x+bottom)],lit:(x+bottom))))
    locLr<-which.max(lr)
    minN<-lr[locLr]*0.999
    maxN<-lr[locLr]
    a<-which.min(lr[lr>minN])
    if(a>locLr)
        a
    else
        locLr
}

linearity<-function(v,window,points=c(),bend=rep(1,length(v)))
{
    a<-seq(window[2]-window[1]-length(points) +1)
    region<-which(! unlist(lapply(window[1]:window[2], function(x) x%in% points)))+window[1]-1
    b<-v[region]/bend[region]
    b<-cdf(b)
    cor(a,b)^2
}

plotPort<-function(i,j, k=1){       
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
        plot(v$x,v$y,xlab=name,col=cols,xlim=c(-200,200))
    v}
plotCDF<-function(i,j,reg="bottom"){
        k<-1
        v<-selectValue(t0f,mots[i],mots[j],loc=k,range=reg)
        if(reg=="bottom")
            plot(cdf(rev(v$y)))
        else if (reg=="top")
            plot(cdf(v$y))
        v
    }

quickSelect<-function(x,y,loc=1)
{
    m1<-mots[x]
    m2<-mots[y]
    name=paste(m1,m2,loc,sep="-")
    tip=t0f[name,]
    x<-as.numeric(colnames(t0f))
    data.frame(y=as.numeric(tip),x=x)
}

plotPause<-function(y,i){
    plotPort(y[i,1],y[i,2],2)
    title(paste(y[i,1],y[i,2],1,sep="-"))
    readKey()
}

readKey<-function(){
    cat("Level of interest: ")
    l<-readline()
    l
}

generateRandomSelects<-function(n,start=1,end=301,comb=2){
    x<-t(combn(start:end,comb))
    y<-x[sample(dim(x)[1],n),]
    r<-c()
    for (i in seq(n)){
        r<-c(r,paste(y[i,1],y[i,2],1,sep="-"))
    }
}

multiSplit<-function(data,n){
    maximum<-max(data)
    minimum<-min(data)
    step<-(maximum-minimum)/n
    seq<-seq(n)*step+minimum
    cbind(seq,apply(data,2,equaSplit,n,maximum,minimum))
}

equaSplit<-function(x,n,maximum=max(x),minimum=min(x)){
    step<-(maximum-minimum)/n
    seq<-seq(n)*step+minimum
    unlist(lapply(seq,function(n){length(which(x>=(n-step) &x<n))}))
}

interq<-function(val4){
    gam<-do.call("-",as.list(quantile(val4)[c(4,2)]))/2
    if(gam==0){
        #l<-length(val4)
        #t<-abs(val4)
        #mean(abs(t[order(t)[1:(0.9*l)]]))
        gam<-(weighted.sum(abs(x$y[300:450]))+weighted.sum(abs(x$y[150:300])))/2
        #mean(abs(val4[order(val4)][((length(val4)*(.5-0.45))):(length(val4)*(.5+0.45))]))
    }
    if(gam==0){
        gam<-0.0001}
    gam
}

primeBetween<-function(a,b){
    seq(a,b)[unlist(lapply(seq(a,b),is.prime))]
}


weighted.sum<-function(x,fil=FALSE )
    {
        if(fil==FALSE){
            d<-seq(floor(length(x)/2))
            if(is.even(length(x))){
                fil<-c(d,rev(d))}
            else{
                fil<-c(d,ceiling(length(x)/2),rev(d))}
            fil<-fil/sum(fil)
        }
        sum(x*fil)
    }

stretchInts<-function(x){
    ints<-do.call(seq,as.list(range(x)))
    unlist(lapply(ints,function(y){
        l<-length(which(x==y))
        (seq(l)/l-0.5)+y}))
}

compareModels<-function(val4){
    val5<-getHeights(val4)
    x<-do.call(seq,as.list(range(val4)))
    y<-val5/sum(val5)    
    plot(x, y)
    points(x,unlist(lapply(x,dnorm,mean(val4),sqrt(var(val4)))),col=2)
    points(x,unlist(lapply(x,function(x2) dgeom(abs(x2) ,max(y)))),col=3)
    cx<-median(val4)
    gam<-interq(val4)
    points(x,unlist(lapply(x,dcauchy,cx,gam)),col=4)
    c(cor(y,unlist(lapply(x,dnorm,mean(val4),sqrt(var(val4))))),
      cor(y,unlist(lapply(x,function(x2) dgeom(abs(x2) ,max(y))))),
      cor(y,unlist(lapply(x,dcauchy,cx,gam))))
}

projectDiff<-function(n1,n2){
    lM<-intersect(intersect(locationsM[[n1]],locationsM[[n2]]),which(reg[,1]))
    lC<-intersect(intersect(locationsC[[n1]],locationsC[[n2]]),which(reg[,1]))
    if(length(lM)>0)
        tbM<-list(bM[[n1]][lM],bM[[n2]][lM])
    if(length(lC)>0)
        tbC<-list(bC[[n1]][lC],bM[[n2]][lC])
    vala<-do.call(rbind,Map(function(x,y){
        x<-unlist(x)
        y<-unlist(y)
        c<-length(x)*length(y)
        a<-c/length(x)
        b<-c/length(y)
        cbind(rep(x,a),rep(y,b))},
                            tbM[[1]],tbM[[2]]))
    valb<-do.call(rbind,Map(function(x,y){
        x<-unlist(x)
        y<-unlist(y)
        c<-length(x)*length(y)
        a<-c/length(x)
        b<-c/length(y)
        cbind(rep(x,a),rep(y,b))},
                            tbC[[1]],tbC[[2]]))
    val1<-rbind(vala,valb)
    val2<-val1%*%matrix(c(-1,1,0,0),ncol=2)
    val3<-val1%*%matrix(c(1,1,0,0),ncol=2)
    val4<-getHeights(val2[,1],c(-300,300))-getHeights(val3[,1],c(0,600))
    val4
}

cauchyTest<-function(val4,gam=interq(val4),cx=median(val4)){
   x<-do.call(seq,as.list(range(val4)))    
   dcauchy(max(x),cx,gam)
}
