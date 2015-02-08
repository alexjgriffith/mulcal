#!/data/binaries/R-3.1.0/bin/Rscript
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
#sets<-c("normal_not_normal.pwm", "abnormal_not_abnormal.pwm")
#shortNames<-c("n_not_a","a_not_n")
#fileLocation<-"~/Dropbox/temp-data/"
#pScoreFile<-"/home/agriffith/Dropbox/temp-data/p-score_table_1.tab"
#imageDirectory="/home/agriffith/Dropbox/temp-data/"
#printFasta(fileLocation,sets,shortNames,imageDirectory)
#write(htmlGenerateMain(fileLocation,sets,shortNames,imageDirectory,n=2),
#       "a-score_other.html")
#
######################################################################
######################################################################


######################################################################
######################################################################

library('getopt')
#loadProject("../src/r/",c("aux.r","html.r"))

loadData<-function(fileLocation){
    data<-read.delim(fileLocation,header=FALSE,sep="\n")    
    out<-list()
    n<-0
    name<-c()
    for(i in seq(length(t(data)))){
        if(  ">" %in% strsplit(as.character(data[i,]),"")[[1]])
            {if (! n==0){out<-append(out,list(matrix(box,dim(box))))}
             box<-c()
             name<-c(name,strsplit(as.character(data[i,]),"\t")[[1]][1])
             n<-n+1}
        else{
            box<-rbind(box,(as.numeric(unlist(strsplit(as.character(data[i,]),"\t")))))}}
    cbind(name=name,data=append(out,list(matrix(box,dim(box)))))}



loadData<-function(fileLocation){
    data<-unlist(read.delim(fileLocation,skip=1,header=FALSE,sep="\n"))
    out<-list()
    n<-0
    name<-c()
    for(i in seq(length(data))){
        d<-unlist(strsplit(as.character(data[i])," "))
        if( ">"==unlist(strsplit(as.character(data[i]),""))[1])
            {if (! n==0){out<-append(out,list(t(matrix(box,dim(box)))))}

             box<-c()
             name<-c(name,d[1])
             n<-n+1}
        else{
            l<-length(d)
            box<-rbind(box,as.numeric(d[seq(3,l-1)]))}}
    cbind(name=name,data=append(out,list(t(matrix(box,dim(box))))))}



palindrome<-function(data){
    makePWM(t(apply(apply(data@pwm,2,rev),1,rev)))}

getScore<-function(pScore,motif,name){
    as.numeric(unlist(pScore[name]))[which(as.character(pScore$motif)==motif)[1]]}

imageList<-function(fileLocation,sets,x,shortNames,imageDirectory,front="motif_"){
    data<-loadData(paste(fileLocation,sets[x],sep=""))
    name<-shortNames[x]
    apply(data,1,function(x){htmlImage(paste(imageDirectory,front,name,"_",strsplit(x[1]$name,">")[[1]][2],".png",sep=""),"style","width:101px;height:50px")})}

scoreList<-function(i,fileLocation,sets,j,shortNames,pScoreFile){
    data<-read.delim(pScoreFile)
    print(length(data))
    pScore<-read.table(pScoreFile,header=1)
    data<-loadData(paste(fileLocation,sets[j],sep=""))
    apply(data,1,function(x){getScore(pScore,x[1]$name,shortNames[i])})}

getMotifs<-function(data){
    apply(data,1,function(x){strsplit(x$name,">")[[1]][2]})}

prepareTable<-function(loadData,shortNames,fileLocation,sets,pScoreFile,n,x){
    temp<-sapply(seq(n),scoreList,fileLocation,sets,x,shortNames,pScoreFile)
    apply(
        cbind(
            getMotifs(loadData(paste(fileLocation,sets[x],sep=""))),
            temp),
        1,function(x){htmlTable(matrix(collapse(as.character(x))),anotations=buildAnotations("style","font-size:10px"))})}

getFrameChar<-function(lis,val){
    as.character(lis[val][[1]])}

sequenceGen<-function(... ,fn, shortNames,n=3){
    I<-lapply(seq(n),function(x) {fn( x=x,shortNames=shortNames,...)})
    ma<-max(sapply(I,length))
    I<-mapply( function(x,y){c(y,rep("",x))},ma-sapply(I,length),I)
    I<-data.frame(I)
    colnames(I)<-shortNames
    I}

printFasta<-function(fileLocation,sets,shortNames,imageDirectory)
    for(j in seq(length(sets))){
        data<-loadData(paste(fileLocation,sets[j],sep=""))
        a<-lapply(seq(length(data[,1])), function(i){makePWM(t(as.matrix(data[i,2]$data)))})
        for( i in seq(length(a))) {
            name<-shortNames[j]
            motif<-strsplit(data[i,1]$name,">")[[1]][2]
            png(paste(imageDirectory,"motif_p_",name,"_",motif,".png",sep=""))
            seqLogo(palindrome(a[[i]]),ic.scale=TRUE,xaxis=FALSE,yaxis=FALSE)
            dev.off()
            png(paste(imageDirectory,"motif_",name,"_",motif,".png",sep=""))
            seqLogo(a[[i]],ic.scale=TRUE,xaxis=FALSE,yaxis=FALSE)
            dev.off()}}

htmlGenerateMain<-function(fileLocation,sets,shortNames,imageDirectory,n=3,pScoreFile){
    I<-sequenceGen(fn=imageList,n=n,shortNames=shortNames,fileLocation,sets,imageDirectory)
    IP<-sequenceGen(fn=imageList,n=n,shortNames=shortNames,fileLocation,sets,imageDirectory,front="motif_p_")
    M<-sequenceGen(fn=prepareTable,n=n,shortNames=shortNames,loadData,fileLocation,sets,pScoreFile,n)
    mat<-matrix(unlist(lapply(seq(n),function(x){lapply(list(I,IP,M),getFrameChar,x)})),ncol=3*n)
    htmlDoc(htmlTags("head",htmlTags("title", "Test Images")),htmlTags("body",c(htmlTable(mat))))}

htmlGenerateHelp <- function(){
    output<-c("Usage:\n\n",
              "htmlGenerate.r \\\n",
              "\t-f ~/Dropbox/temp-data/ \\\n",
              "\t-i /home/agriffith/Dropbox/temp-data/ \\\n",
              "\t-p /home/agriffith/Dropbox/temp-data/p-score_table_1.tab \\\n",
              "\t-s normal_not_normal.pwm,abnormal_not_abnormal.pwm \\\n",
              "\t-n n_not_n,a_not_a -o test.html\n\n",
              "The location of the sets\t-f --filelocation \n",
              "The location of the images\t-i --imageDirectory \n")
    write(lcollapse(output),stdout())
    quit()}

main<-function(){
    spec = matrix(c(
        'fileLocation','f', 1,"character",
        'imageDirectory','i', 1,"character",
        'pScoreFile','p', 1,"character",
        'sets','s', 1,"character",
        'outFile','o', 1,"character",
        'shortNames','n',1,"character",
        'help','h',0,"logical",
        'generateLogos',"g",0,"logical"
    ),byrow=TRUE,ncol=4)
    args=getopt(spec)
    if(! is.null(args$help)){htmlGenerateHelp()}
    if(is.null(args$generateLogos)){args$generateLogos=FALSE}
    with(as.list(args),{
        sets<-strsplit(sets,",")[[1]]
        shortNames<-strsplit(shortNames,",")[[1]]
        n=length(sets)
        if(generateLogos)
            {library("seqLogo")
             printFasta(fileLocation,sets,shortNames,imageDirectory)}
        write(
            htmlGenerateMain(fileLocation,sets,shortNames,imageDirectory,n,pScoreFile)
           ,outFile)})}

main()

