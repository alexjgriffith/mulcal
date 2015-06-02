#!/usr/bin/env R
#
# This file is part of peakAnalysis,
# http://github.com/alexjgriffith/alpha-score/, 
# and is Copyright (C) University of Ottawa, 2015. It is Licensed under 
# the three-clause BSD License; see LICENSE.txt.
# Author : Alexander Griffith
# Contact: griffitaj@gmail.com
#

#' @export
loadHeightFile<-function(file="~/masters/normal-abnormal/single_heights.bed",n=3,norm=TRUE){  
  cdata<-read.table(file)
  stats<-cdata[seq(n)]
  l=length(cdata)
  if(norm==TRUE)
      data<-as.matrix(apply(cdata[(n+1):l],2, function(x) as.vector((unlist(x/(stats[3]-stats[2]))))))
  else
      data<-as.matrix(apply(cdata[(n+1):l],2, function(x) as.vector((unlist(x)))))
  list(stats=as.matrix(stats),
            data=data)}

#' @export
qn <-function(data){
    shape<-dim(data)
    sequence<-apply(data,2,order)
    reverseSequence<-unlist(apply(sequence,2,order))
    ranks<-apply(matrix(unlist(lapply(seq(shape[2]),function(i,x,y) x[y[,i],i],data,sequence)),ncol=shape[2]),1,sum)/shape[2]
    apply(reverseSequence,2,function(x) ranks[x])}

#' @export
#' 
standPCAPrep <-function(data,v="colQn"){
  switch(v,
         rowSumOne=t(apply(data,1, function(x) {x/sum(x)})),
         colSumOne=apply(data,2, function(x) x/sum(x)),
         row=t(apply(data,1, function(x) (x-mean(x))/var(x))),
         col=apply(data,2, function(x) (x-mean(x))/var(x)),
         rows1=t(apply(data,1, function(x) x/var(x))),
         cols1=apply(data,2, function(x) x/var(x)),
         colQn=qn(data),
         non=t(data))}

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

#' @export
plotPCs<-function(pcs,pos,data,cats,lab=c("xlabel","ylable","Title")){
    if("rotation" %in% names(pcs))
        x<-t(as.matrix(pcs$rotation)) %*% as.matrix(data)
    else
        x<-t(as.matrix(pcs)) %*% as.matrix(data)
    d1<-data.frame(x[pos[1],])
    d2<-data.frame(x[pos[2],])
    plot(t(d1),t(d2),xlab=lab[1],ylab=lab[2])
    title(main=lab[3])
    text(x[pos[1],],x[pos[2],],labels=cats,cex=0.7,pos=3)}

#' @export
interClass<-function(pca,data,pos,C=3,r=13){
  sm<-data %*%  pca$rotation[,pos]
  m<-mean(sm)
  v<-var(sm)
  S<-which(unlist(lapply(sm,function(x) x< m-v*C)))
  G<-which(unlist(lapply(sm,function(x) x> m+v*C)))
  data.frame(data=unlist(lapply(seq(r),function(x) mean(data[G,x])/mean(data[S,x]))),row.names=cats)}

#' Principle Componenet Analysis
#' @export
#' @examples
#'    peakFile <- "sorted_peaks.bed"
#'    catFile <- "catagories"
#'    pcs<-pcaAnalysis(peakFile,catFile)
#'    pc1=3
#'    pc2=6    
#'    plotPCs(pcs,cbind(pc1,pc2),normData,cats,c( paste("PC:",as.character(pc1),sep="" ),paste("PC:",as.character(pc2),sep="" ),mock[i]))
pcaAnalysis<-function(heightFile,catFile,normMethod="colQn")
    {
    values<-loadHeightFile( heightFile)
    catagories<-t(read.table(catFile))
    stats<-values$stats
    data<-values$data
    colnames(data)<-catagories
    normData<-standPCAPrep(data,normMethod)
    list(pcs=prcomp(t(normData)),data=normData,cats=catagories,stats=stats)
    }


ascoreWeighting<-function(pcs,scores,weights){
    t(matrix(weights))%*%t(as.matrix(pcs$rotation[,scores]))}

#' The Processing of PC values
#' A utility designed to help seperate a vector based on the mean and standard deviation of the data. For the most part ascoreSeparation will return a vector of logicals but funs such as "value" and "rank" return numerics.
#' @template authorTemplate
#' @param vect A vector of numerical values.
#' @param fun A function which is used to sepearate memebers of the vector. 
#' \describe{
#' \item{\strong{"top"}}{TRUE if greater than mean + n* standard deviation.}
#' \item{\strong{"bottom"}}{TRUE if less than mean - n* standard deviation.}
#' \item{\strong{"nottop"}}{TRUE if less than mean + n* standard deviation.}
#' \item{\strong{"notbottom"}}{TRUE if greater than mean - n* standard deviation.}
#' \item{\strong{"middle"}}{TRUE if within mean +/- n* standard deviation.}
#' \item{\strong{"value"}}{Returns the principle component of interest.}
#' \item{\strong{"rank"}}{Returns the ranks of the principle components.}
#' }
#' @param n The value which is used to adjust the standard deviation cut off.
#' @return A vector of equal length to the parameter vect.
#' @export
ascoreSeperation<-function(vect,fun,n=1){
    loopfuns<-c("top","bottom","nottop","notbottom","middle","value")
    if(fun %in% loopfuns){
        ap<-switch(fun,
                   top=function(x,m,sd,n){x>m+sd*n},
                   bottom=function(x,m,sd,n){x<m-sd*n},
                   nottop=function(x,m,sd,n){x<=m+sd*n},
                   notbottom=function(x,m,sd,n){x>=m-sd*n},
                   middle=function(x,m,sd,n){all(x>=m-sd*n,x<=m+sd*n)},
                   value=function(x,m,sd,n){x})
        m<-mean(vect)
        sd<-sqrt(var(vect))
        sapply(vect,ap,m,sd,n)}
    else if(fun =="ranked"){
        order(vect)}
    else
        vect}


#' Batch alpha score seperation
#' Apply a serise of options to ascoreSeperation and return a matrix. Used by \code{\link{ascore()}}.
#' @seealso \code{\link{ascoreSeperation}}
#' @seealso \code{\link{ascore}}
#' @export 
batchscore<-function(pcs,scores,funs,ns){
    mapply(function(x,y,z){ascoreSeperation(pcs$rotation[,x],y,z)},
           scores,funs,ns)    }

#' Alpha Score
#'
#' Returns a seperating vector built from principle components. returns a matrix with a length of the data set inputed in heightFile and a width of the lenght of the number of options. As such a wide variety of seperating mechanisms can be called at the same time automating the process of isolating regions of data.
#' 
#' @param heightFile Can be either a file name, a data.frame, or a matrix
#' @param pc The principle component of interest, note the length of pc, funs and ns must be the same.
#' @param funs funs=c("top","bottom","nottop","notbottom","middle","ranked","value")
#' @param ns An integer to adjust the SD in the fun
#' @param normMethod Default="colQn", use "non" if no normalization is desired
#' @seealso For more infromation on callable functions \code{\link{ascoreSeperation}}
#' @return An nxm matrix where n=length of data and m=number of optoins.
#' @export
#' @examples
#' 
#' # load a bed file of format
#' # chr start end h1 ... hn
#' heightFile<-"some_bed_file.bed"
#' 
#' # loadHeightFile splits the first three columns into
#' # stats and the remainder into data
#' # You can use any utiltity you wish to load data so
#' # long as it clearly seperates the data (heights, read pile up,
#' # etc) from the stats (chromosome inforomation, data source, genes etc)
#' values<-loadHeightFile(heightFile)
#' data<-values$data
#' stats<-values$stats
#' 
#' pc<-1
#' function<-"top" # x>mean+sd*n
#' n<-3
#' reg<-ascore(data,pc,function,n) 
#' stats[reg]
#' 
#' # Determining two values at the same time
#' reg2<-ascore(data,c(1,2),c("top","top"),c(3,1))
#' regions<-apply(reg2,1,all)
#' stats[regions]
#'
#' # Saving the first 3 principle components
#' pcs<-ascore(data,c(1,2,3,4))
#'
#' reg4<-ascoreSeperation(pcs[1,],"top",3)
#' 
#' @template authorTemplate
ascore<-function(heightFile,pc,funs=rep("value",length(pc)),ns=rep(0,length(pc)),normMethod="colQn"){
    if(is.character(heightFile)){
        values<-loadHeightFile(heightFile)
        data<-values$data}
    else if(is.data.frame(heightFile))
        data<-as.matrix(heightFile)
    else if(! is.matrix(heightFile))
        stop("heightFile must be either a filename data.frame or matrix.")
    normData<-standPCAPrep(data,normMethod)
    pcs=prcomp(t(normData))
    batchscore(pcs,pc,funs,ns)}

#' Principle Componenet Analysis Test
#'
#' @export
pcaAnalysisTest<- function(pc1=1, pc2=3)
{
    root_file<-"~/Dropbox/UTX-Alex/jan/"
    mock<-cbind("combined")
    i=1    
    peakFile<- paste(root_file,mock[i],"_heights.bed",sep="" )
    catFile<-paste(root_file,"catagories",sep="" )
    pcsData<-pcaAnalysis(peakFile,catFile)
    pcs<-pcsData$pcs
    normData<-pcsData$data
    cats<-unlist(pcsData$cats)
    p3<-ascoreSeperation(pcs$rotation[,3],"top",3)
    p1<-ascoreSeperation(pcs$rotation[,1],"top",3)
    batch<-batchscore(pcs,c(1,3),c("top","top"),c(3,3))
    print(str(batch))
    length(which(apply(batch,2,all)==TRUE))
    pc1=1
    pc2=3
    #length(which(ascoreSeperation(vect,{x<m-sd*3})==TRUE))
    #plotPCs(pcs,cbind(pc1,pc2),normData,cats,c( paste("PC:",as.character(pc1),sep="" ),paste("PC:",as.character(pc2),sep="" ),mock[i]))}
}

#' some hierarchichal clustring  Tests
#'
#' @export
hcTest<-function(){

    cats2<-data.frame(tall_p1="Prima5",
                      tall_p2_1="Prima2",
                      tall_p2_2="Prima2",
                      tall_p3_1="Prima5",
                      tall_p3_2="Prima5",
                      jurk_sandar_1="Jurkat",
                      jurk_sandar="Jurkat",
                      jurk="Jurakt",
                      rpmi_1="RPMI",
                      rpmi_2="RPMI",
                      cem_1="CEM",
                      cem_2="CEM",
                      "ecfc-tsa"="ECFC",
                      meka="Meka",
                      cd133="HSC",
                      cd34="HSC",
                      cd34_new="HSC",
                      eryt="Erythroid",
                      eryt_f="Erythroid",
                      eryt_a="Erythroid",
                      k562_1="K562",
                      k562_2="K562")
    
    root_file<-"~/Dropbox/UTX-Alex/jan/"
    mock<-cbind("single")
    i=1
    values<-loadHeightFile( paste(root_file,mock[i],"_heights.bed",sep="" ))
    cats<-t(read.table(paste(root_file,"catagories",sep="" )))
    cats3<-unlist(lapply(cats,function(x){as.character(cats2[[gsub("-",".",x)]])}))
    data<-values$data
    temp<-cor(data)
    rownames(temp)<-cats3
    colnames(temp)<-cats3
    #png("~/Dropbox/dend_corr_average_combiend.png")
    dend=hclust(dist(temp),method="average")
    plot(dend,hang=-1,main="Combinded Mock Dendogram",xlab="Cell Type")
    require(dendextend)
    d1<-as.dendrogram(hclust(dist(1-temp),method="single"))
    c1<-color_branches(d1,k=4,col=c("green","red","blue","orange"))
    plot(c1)
    #dev.off()
    #png("~/Dropbox/heatmap.png")
    heatmap(1-log(cor(values$data)) ,hclustfun=function(x)hclust(x,method="complete") ,distfun=function(x)as.dist(x),scale="column")
    #dev.off()
}

