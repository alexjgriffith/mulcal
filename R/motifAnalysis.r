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

### Non Finalized Functions
#' Score List
#'
#' The preliminary investigation score mechanism
#' scoreList(lapply(h,function(x) lapply(x,function(y)y[[2]])),mots,addmotifs)
#' @export
scoreList<-function(k,mot,addmotifs){
    which(unlist(mapply(function(k,mot){
        if(length(k)>0){
            t0<-unlist(lapply(seq(min(k),max(k)),function(i)length(which(k==i))))
            vs<-unlist(lapply(min(t0):max(t0),function(x)length(which(x==t0))))
            t<-(vs[length(vs)]> (vs[1]*2^(-length(vs)/sqrt(var(t0)))))
            if(is.na(vs[1]*2^(-length(vs)/sqrt(var(t0)))))               {
                    return(FALSE)
                }
            else if(t)
                if(all(c(vs[1]*2^(-length(vs)/sqrt(var(t0)))<1,length(vs)>5,nchar(addmotifs)<abs(which.max(t0)+min(k)))))
                    return(TRUE)
                else
                    return(FALSE)
        }
        return(FALSE)},h,motifs )))
}

#' Prelim socre function
#' 
#' @export
motifScoreFunction<-function(x) {
    if (length(x)<2){return(0)} else{ if(var(x)==0){ 0} else {t0<-getHeights(x);length(t0)*2^(-max(t0)/sqrt(var(t0)))}}}

#' Get Heights
#'
#' Get the heights of the motif comparison data
#'
#' @export
getHeights<-function(h,range=c(min(h),max(h))){
    rep<-rep(0,(range[2]-range[1]+1))
    for(i in h)
        rep[i-range[1]+1]<-rep[i-range[1]+1]+1
    rep
}
    #unlist(lapply(range,function(i)length(which(h==i))))}

#' geometric score value
#'
#' @export
geometricScore<-function(h){
    n<-length(h)/600
    l<-seq(min(h),max(h))
    t0<-getHeights(h)
    t1<-unlist(lapply(l,function(x) 1+n*(abs(x)/300)/6))
    p<-1/(1+mean(t0))
    test<-function(x,p){(1-p)^x*p}
    plot(unlist(lapply(t0,test,p))*t1)
    plot(unlist(t0)*t1)
}

#' Intersect Count
#'
#' Determine the number of cases which 2 motifs intersect
#'
#' @export
numberIntersect<-function(motifA,motifB,Sequences,reg=rep(TRUE,length(Sequences))){
length(
    intersect(
    
    intersect(grep(unlist(IUPACtoBase(motifA)),Sequences),grep(unlist(IUPACtoBase(motifB)),Sequences))
  , which(reg))
    )
}


### Application Functions
#' Motifs Two View
#' 
#' Visualize the similarity between two motifs,
#' Very ineficent,
#'
#' @export
#' @examples
#' FastaFile<-"~/Dropbox/UTX-Alex/jan/combined.fasta"
#' Sequences <- readDNAStringSet(FastaFile, "fasta")
#' data<-loadHeightFile("~/Dropbox/UTX-Alex/jan/combined_heights.bed")$data
#' reg<-as.matrix(ascore(data,c(1,1,1),c("top","bottom","middle"),c(1,1,1)))
#' m<-motifs2View("CANNTG","TGACCT",reg[,3],Sequences)
#' k<-motifs2View("CANNTG","GATAAG",reg[,1],Sequences)
#' t0<-unlist(lapply(seq(min(k),max(k)),function(i)length(which(k==i))))
motifs2View<-function(m1,m2,reg,Sequences,nearHeights=FALSE){
    m12<-c(m1, m2)
    mList<-unlist(lapply(m12,IUPACtoBase))
    cList<-lapply(unlist(lapply(m12,IUPACtoBase)),compliment)
    locationsM<-lapply(mList,grep,Sequences)
    locationsC<-lapply(cList,grep,Sequences)
    if(nearHeights==FALSE){
    h<-list(motifHist(Sequences,mList,cList,locationsM,locationsC,1,2,reg))
    histVisualize(h,m1,m2)
    h<-unlist(h)
    t0<-unlist(lapply(seq(min(h),max(h)),function(i)length(which(h==i))))}
    else{
        h<-nearSummit(Sequences,mList,cList,locationsM,locationsC,1,reg)}
    #c((max(t0)-mean(t0))/sqrt(var(t0)),scoreFunction(t0))
    h
}

#' Histogram Wrapper
#'
#' @export
#' @examples
#' FastaFile<-"~/Dropbox/UTX-Alex/jan/combined.fasta"
#' data<-loadHeightFile("~/Dropbox/UTX-Alex/jan/combined_heights.bed")$data
#' reg<-as.matrix(ascore(data,c(1,1,1),c("top","bottom","middle"),c(1,1,1)))
#' Sequences <- readDNAStringSet(FastaFile, "fasta")
#' motifFile<-"inst/data/abnormal_normal.pwm"
#' motifs<-homerWrapper(Sequences,reg[,1],reg[,2],"~/Masters/mulcal/inst/lib/homer-4.7/cpp/homer2",motifFile)
#' motifs<-loadPWM(motifFile)
#' mots<-unlist(lapply(unlist(motifs[,1]),function(x)gsub(">","",x)))
#' motifs<-loadPWM("inst/data/jaspar_motifs.txt","jaspar")
#' mots<-c(mots,unlist(lapply(1:length(motifs[,2]),function(i)PWMtoCons(motifs[i,2]$data))))
#' addmotifs<-c("CGNNGC")
#' mList<-unlist(lapply(c(mots,addmotifs),IUPACtoBase))
#' cList<-unlist(lapply(lapply(c(mots,addmotifs),IUPACtoBase),compliment))
#' locationsM<-lapply(mList,grep,Sequences)
#' locationsC<-lapply(cList,grep,Sequences)
#' l<-length(cList)
#' motifHist(Sequences,mList,cList,locationsM,locationsC,4,l,reg)
#' 
#' h<-lapply(seq(3),function(k)lapply(seq(1:l-1),function(i)motifHist(Sequences,mList,cList,locationsM,locationsC,i,l,reg[,k])))
#' 
#' h<-lapply(seq(10),function(k)lapply(seq(10),function(i)lapply(seq(3),function(j)motifHist(Sequences,mList,cList,locationsM,locationsC,i,k,reg[,j]))))
#' lapply(h,scoreList,c(mots,addmotifs)[1:10],addmotifs)
#' histWrapper(lapply(h,function(x) lapply(x,function(y)y[[2]])),c(mots,addmotifs), names=c("NonName"),c(1),c(30))
histWrapper<-function(h,motifNames,names=c("NoName"),...){
    BuildernestedList<-function(h,locs){
        fun<-function(i,h,...){NULL}
        function(funbody,...){
            if(is.expression(funbody))
                body(fun)<-funbody
            else
                fun<-funbody
            subset.nestedList(h,locs,fun,...)}
    }
    D2<-BuildernestedList(h,list(...))
    vals<-D2(expression(h[[i[1]]][[i[2]]][[i[3]]]))
    mot<-D2(function(i,h,motifNames) motifNames[i],motifNames=motifNames)
    motx<-Map(function(x)x[1],mot)
    moty<-Map(function(x)x[2],mot)
    histVisualize(
        vals,
        motx,moty,prod(rapply(list(...),length)),names)
}


### Utility Functions
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
    IUPACCharacters<-list("A","C","G","T",c("A","G"),c("C","T"),c("C","G"),c("A","T"),c("G","T"),c("A","C"),c("C","G","T"),c("A","G","T"),c("A","C","T"),c("A","C","G"),c("A","C","G","T"))
    names(IUPACCharacters)<-c("A","C","G","T","R","Y","S","W","K","M","B","D","H","V","N")

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
 motifHist<-function(data,mList,cList,locationsM,locationsC,n1,n2,reg,width=c(-30,30)){
    lM<-intersect(intersect(locationsM[[n1]],locationsM[[n2]]),which(reg))
    lC<-intersect(intersect(locationsC[[n1]],locationsC[[n2]]),which(reg))
    bM<-c()
    bC<-c()
    if(length(lM)>0)
        bM<-lapply(c(mList[n1],mList[n2]),function(x) lapply(gregexpr(x, data[lM]),as.numeric))
    if(length(lC)>0)        
        bC<-lapply(c(cList[n1],cList[n2]),function(x) lapply(gregexpr(x, data[lC]),as.numeric))
    #h<-c(getDistance(bM[[1]],bM[[2]]),-getDistance(bC[[1]],bC[[2]]))
    #print(c(mList[n1],consensusIUPAC(mList[n2]))
    #print(nchar(mList[n1])-nchar(mList[n2]))
    print(nchar(consenusIUPAC(mList[n1]))-nchar(consenusIUPAC(mList[n2])))
    h<-c(getDistance(bM[[1]],bM[[2]]),-(nchar(consenusIUPAC(mList[n1]))-nchar(consenusIUPAC(mList[n2]))+getDistance(bC[[1]],bC[[2]])))
    #h<-getDistance(bM[[1]],bM[[2]])#getDistance(bM[[1]],bM[[2]])#,-getDistance(bC[[1]],bC[[2]]))
    #if(length(h)>0)
       # hist(h,breaks=1000,xlim=width,xlab=paste(mList[n1],mList[n2],sep=" - "))
    h}


#'@export
nearSummit<-function(data,mList,cList,locationsM,locationsC,n1,reg,width=150)
    {
    lM<-intersect(locationsM[[n1]],which(reg))
    lC<-intersect(locationsC[[n1]],which(reg))
    bM<-c()
    bC<-c()
    if(length(lM)>0){
        bM<-lapply(gregexpr(mList[n1], data[lM]),as.numeric)
        #print(mList[n1])
    }
    if(length(lC)>0)        {
        bC<-lapply(gregexpr(cList[n1], data[lC]),as.numeric)
        #print("b")
    }
    h<-c(getDistance(unlist(bM),width),getDistance(unlist(bC),width))
    h}

#' @export
getDistance<-function(x,y,one=FALSE){
    as.numeric(
        unlist(mapply(function(x,y){
        temp<-outer(x, y,"-")
        temp<-temp[upper.tri(x=temp,diag=TRUE)]
        #print(str(temp))
        if(one ){temp[which.max(abs(temp))]}
        else  {temp}
    },x,y)))}

#' Homer Wrapper
#'
#' This is only a wrapper for the DENOVO functionality of homer2, homer must be installed on the system before this can be used.
#' 
#' @export
homerWrapper<-function(sequences,foreground,background,homerLocation,motifsFile=FALSE,opts="-S 25 -len 6"){
    if (!is.character(motifsFile))
        motifsFile<-tempfile()
    treatmentFile<-tempfile()
    controlFile<-tempfile()    
    writeXStringSet(sequences[foreground],treatmentFile)
    writeXStringSet(sequences[background],controlFile)
    cmd<-paste(homerLocation, "denovo -i ",treatmentFile," -b ",controlFile,opts," > ",motifsFile,sep=" ")
    system(cmd)
    loadPWM(motifsFile,"homer")
}

#' PWM to Consensus Motif
#'
#' A caller for ConsesnusIUPAC, Recives PWM and returns a consusensus motif
#' The PWMs should be loaded through the loadPWM utility.
#'
#' @seealso \link{\code{loadPWM()}}
#' 
#' @export
PWMtoCons<-function(x)
    consenusIUPAC(motifString(x))
#' Consensus IUPAC
#'
#' 
#' @export
#' 
consenusIUPAC<-function(mstring){
        IUPACCharacters<-list("A","C","G","T",c("A","G"),c("C","T"),c("C","G"),c("A","T"),c("G","T"),c("A","C"),c("C","G","T"),c("A","G","T"),c("A","C","T"),c("A","C","G"),c("A","C","G","T"))
        names(IUPACCharacters)<-c("A","C","G","T","R","Y","S","W","K","M","B","D","H","V","N")
        IUPACc<-unlist(lapply(IUPACCharacters,paste,collapse=""))
        x<-splitBlist(mstring)
        paste(lapply(x,function(x)names(which(x==IUPACc))),collapse="")}

#' Split by List
#'
#' Applied to DNA nucleotides (ACGT). Returns a list with [] removed
#' 
#' @export
#' @examples
#' > consensus<-"AGCT[AGCT]G"
#' > splitBlist(consensus)
#'      buf buf buf buf buf    buf
#' [1,] "A" "G" "C" "T" "AGCT" "G"
#'
splitBlist<-function(mstring){
    ret<-c()
    buf<-""
    data<-strsplit(mstring,"")[[1]]
    i<-1
    while (i <= length(data)){
        buf<-""
        if (data[i]=="["){
            while(data[i]!="]"){
                if (data[i] %in% c("A","C","G","T"))
                    buf<-paste(buf,data[i],sep="")
                i<-i+1}}
        else if (data[i] %in% c("A","C","G","T"))
            buf<-data[i]
        ret<-cbind(ret,buf)
        i<-i+1
    }
    ret
}

#' Motif String
#'
#' A scoring mechanism to convert PWM to strings
#' 
#' @export
motifString<-function(x){
    paste(apply(x,2,function(x){
        DNA<-c("A","C","G","T")
        if(sum(x)==0)
            return (paste(c("[",sort(DNA),"]"),collapse=""))
        x<-x/sum(x)
        or<-order(x,decreasing = TRUE)
        t<-or
        if(x[or[1]]>=.6){return(DNA[or[1]])}
        else if(sum(x[or[1:2]])>=0.8){t<-or[1:2]}
        else if((sum(x[or[1:3]])>=0.95)){t<-or[1:3]}
        paste(c("[",sort(DNA[t]),"]"),collapse="")
    }
                ),collapse="")
}

### Plotting Functions

#' @title Visualize batches of height data
#' @export
histVisualize<-function(h,m1,m2,n=1,name=c("NoName")){
    pushViewport(viewport(layout=grid.layout(n,2)))
    if(n==1){
        h<-h[[1]]
        t0<-unlist(lapply(seq(min(h),max(h)),
                          function(i)length(which(h==i))))
        p1<-heightHist(t0,which.max(t0)+min(h)-1)
        p2<-locHist(h,paste(m1,m2,sep="-"))
        print(p1,vp=viewport(layout.pos.row=1,layout.pos.col=1))
        print(p2,vp=viewport(layout.pos.row=1,layout.pos.col=2))
    }
    else{
        for( k in seq(n)){
            t0<-unlist(lapply(seq(min(h[[k]]),max(h[[k]])),
                              function(i)length(which(h[[k]]==i))))
            p1<-heightHist(t0,which.max(t0)+min(h[[k]])-1)+ylab(name[k])
                tm2<-m2[k]
                tm1<-m1[k]
            p2<-locHist(h[[k]],paste(tm1,tm2,sep="-"))+ylab("")
            print(p1,vp=viewport(layout.pos.row=k,layout.pos.col=1))
            print(p2,vp=viewport(layout.pos.row=k,layout.pos.col=2))
        }
    }
}

#' @title Height Histrogram
#'  A plotting tool for histograms, Requires GGPLOT2
#' @export
heightHist<-function(t0,xlab="Histogram"){
    p<-ggplot(as.data.frame(t0),aes(x=t0))+geom_histogram(binwidth=1,xlab=xlab)
    if((max(t0)-min(t0))<20)
        p<-p+scale_x_continuous(breaks=seq(min(t0), max(t0),1))
    p<-p+stat_function(fun=function(x) length(t0)*2^(-x/sqrt(var(t0))),colour="blue" )
    p<-p+xlab(xlab)+ylim(c(0,length(which(t0==0))))
    p}

#' @title Location Histogram
#' A plotting tool for histograms, Requires GGPLOT2
#' @export
locHist<-function(t0,xlab="Histogram",limits=c(-32,32)){
    x<-seq(min(t0),max(t0))
    t1<-unlist(lapply(x,function(x)length(which(x==t0))))
    locHist2(t1,x,xlab,limits)
    #p<-qplot(x,t1,geom="step",ylab="frequency",xlab=xlab)+stat_function(fun=function(x){0})#+xlim(limits)
    #n<-4
    #m<-ceiling(log2(abs(limits[2]-limits[1])/2))
    #s<-seq(max(m-3,1),m)
    #b<-c(-sapply(s,function(x)2^x),sapply(s,function(x)2^x))
    #p+scale_x_continuous(breaks=b,limits=limits)
}

#' @export
locHist2<-function(t1,x,xlab="Histogram",limits=c(-32,32)){
    p<-qplot(x,t1,geom="step",ylab="frequency",xlab=xlab)+stat_function(fun=function(x){0})#+xlim(limits)
    n<-4
    m<-ceiling(log2(abs(limits[2]-limits[1])/2))
    s<-seq(max(m-3,1),m)
    b<-c(-sapply(s,function(x)2^x),sapply(s,function(x)2^x))
    p+scale_x_continuous(breaks=b,limits=limits)
}

### Import Functions

#' @title Load PWM
#' This utiltiy is designed to load several PWM formats and create a uniform layout within R for analysis
#'@export
loadPWM<-function(fileLocation,version="homer")
{
    loadDataBuilder<-function(splitfun,header=FALSE,skip=1,id=">",split=" ",region=1){
        function(fileLocation){
            data<-read.delim(fileLocation,header=header,skip=skip,sep="\n")
            out<-list()
            n<-0
            name<-c()
            info<-c()
            for(i in seq(length(t(data)))) {
                d<-strsplit(as.character(data[i,]),split)[[1]]
                l<-length(d)
                if(  id == strsplit(d[1],"")[[1]][region])
                    {
                        if (! n==0){out<-append(out,list(matrix(box,4)))}
                        box<-c()
                        name<-c(name,d[1])
                        info<-c(info,d[2:l])
                        n<-n+1}
                else{
                    box<-c(box,splitfun(d,l))
                }}
            cbind(name=name,info=info,data=append(out,list(matrix(box,4))) )}}    
    loadFunctions<-list(
        homer=c(
            skip=0,
            split="\t",
            splitfun=function(d,l){
                as.numeric(unlist(strsplit(as.character(d),"\t")))}),
        jaspar=c(skip=0,
            splitfun=function(d,l){
                na.omit(as.numeric(unlist(strsplit(unlist(d), "[^0-9]+"))))
        }))
    if (version %in% names(loadFunctions))
        do.call(loadDataBuilder,loadFunctions[[version]])(fileLocation)
    else
        warning(paste(version,"is not an option."))
}
