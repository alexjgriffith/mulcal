mang<-function(x){
    sqrt(sum(sapply(x,"^",2)))
    }

circle<-function(x,y,r,...){
    n<-1000
    li<-seq(n)/n*2
    xx<-x+r*cos(li*pi)
    yy<-y+r*sin(li*pi)
    polygon(xx,yy,...)
}

lensAux<-function(x,y,r,ran){
    xd<-ran[2]-ran[1]
    yd<-ran[4]-ran[3]
    print(m<-apply(matrix(ran,ncol=2),1, function(x) x/mang(x))%*%c(x,y))
    xd2<-(ran[2]+ran[1])/2
    yd2<-(ran[4]+ran[3])/2
    #print(c(xd,yd))
    maxn<-atan2((yd2-y),(xd2-x))
    d<-mang(c(xd,yd))    
    diff<-acos(1-d^2/(2*r^2))
    n<-10000
    print(unlist(lapply(list(maxn,diff),function(x) c(x,180/pi*x))))
    li<-(seq(n)/n-0.5)*diff+maxn
    li
}

lens<-function(x1,y1,r1,x2,y2,r2,...){
    ran<-intersectCircle(x1,y1,r1,x2,y2,r2)
    li<-lensAux(x1,y1,r1,ran)
    li2<-lensAux(x2,y2,r2,ran)
    xx1<-x1+r1*cos(li)
    xx2<-x2+r2*cos(li2)
    yy1<-y1+r1*sin(li)
    yy2<-y2+r2*sin(li2)
    #li2<-li2[1:5000]
    xx<-c(xx2,xx1) 
    yy<-c(yy2,yy1)
    polygon(xx,yy,...)    
}

intersectCircle<-function(x1,y1,r1,x2,y2,r2){
    l1<-x2-x1
    l2<-y2-y1
    #print(c(l1,l2))
    pheta<-atan2(l2,l1)#if(l1==0){pi*(if(l2<0){1}else{3})/2}
           #else{atan(l2/l1)}+(if(l1>0){0}else{pi})+(if(l2==0) {pi}else {0})    
    d<-sqrt(l1^2 +l2^2)
    xm<-(-r2^2+r1^2 +d^2)/(2*d)
    ym<-sqrt(r1^2-xm^2)
    theta<-if(xm==0) {pi/2} else atan(ym/xm)
    xmp<-x1+r1*cos(pheta+c(theta,-theta)) 
    ymp<-y1+r1*sin(pheta+c(theta,-theta))
    #circle(xmp[1],ymp[1],0.1,col="green")
    #circle(xmp[2],ymp[2],0.07,col="blue")
    c(xmp,ymp)
    }

ven2<-function(x1,y1,r1,x2,y2,r2){
    y<-c(y1,y2)
    x<-c(x1,x2)
    r<-c(r1,r2)
    
    minum<-function(x,r,cont,pos)
        pos(do.call(cont,list(x,r)))

    radiusRange<-function(x,r)
        mapply(function(cont,pos)minum(x,r,cont,pos), c("-","+"), c(min,max))
    
    hy<-radiusRange(y,r*1.1)
    hx<-radiusRange(x,r*1.1)
    print(list(hy,hx))
    plot(hx[2]*2,hy[2]*2,ylim=hy,xlim=hx)
    circle(x1,y1,r1,col="yellow",border="black")
    circle(x2,y2,r2,col="blue",border="black")
    lens(x1,y1,r1,x2,y2,r2,col="red",border="red")
}

ven3<-function(x1,y1,r1,x2,y2,r2,x3,y3,r3){
        plot(10,10,ylim=c(-4,9),xlim=c(-4,9))
        circle(x1,y1,r1)#,col="red",border="black")
        circle(x2,y2,r2)#,col="blue",border="black")
        circle(x3,y3,r3)
        ran1<-intersectCircle(x1,y1,r1,x2,y2,r2)
        ran2<-intersectCircle(x1,y1,r1,x3,y3,r3)
        ran3<-intersectCircle(x3,y3,r3,x2,y2,r2)
        li1<-lensAux(x1,y1,r1,c(ran1[2],ran2[2],ran1[4],ran2[4]))
        li3<-lensAux(x3,y3,r3,c(ran2[2],ran3[1],ran2[4],ran3[3]))
        li2<-lensAux(x2,y2,r2,c(ran3[1],ran1[2],ran3[3],ran1[4]))
        xx1<-x1+r1*cos(li1)
        yy1<-y1+r1*sin(li1)
        xx3<-x3+r3*cos(li3)
        yy3<-y3+r3*sin(li3)
        xx2<-x2+r2*cos(li2)
        yy2<-y2+r2*sin(li2)
        c(ran1[2],ran2[2],ran1[4],ran2[4])
        polygon(c(xx1,rev(xx2),xx3),c(yy1,rev(yy2),yy3),col="red")

    }


