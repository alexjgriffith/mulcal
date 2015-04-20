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


#' Exons from Genome
#' 
#' Generates a bed file format for exons from a RefSeqGene set.
#' 
#' @export
#' @examples
#' genome<-read.table("hg19.RefSeqGenes.csv",header=TRUE)
#' exonList<-exonsFromGenome(genomes)
#' write.table(exonLocs,"hg19.exonLocs",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
exonsfromGenome<-function(genome,chroms=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)){
    exons=as.data.frame(t(matrix(unlist(apply(genome[,c('chrom','exonStarts','exonEnds','exonCount','name' )],1,function(x){
    lapply(seq(as.numeric(x[4])),function(y){
           c(x[1],strsplit(as.character(x[2]), ",")[[1]][y],strsplit(as.character(x[3]), ",")[[1]][y],x[5])})  })),nrow=4 )))
    chroms<-unlist(lapply(chroms,function(x){paste("chr",as.character(x),sep="")}))
    exons[as.character(exons[,1]) %in% chroms,]
    }

# work in progress
rnaSeqAnalysis<-function(){
    geneomes<-read.table("inst/data/hg19.RefSeqGenes.csv",header=TRUE,comment.char = "")
    vals<-loadHeightFile("inst/data/rna_heights3.bed",n=4)
    cats<-read.table("inst/data/rna_cats")
    data<-do.call(rbind,mapply(function(genes,chro){
        levs<-as.character(genes)
        reg1<-vals$stats[,1]==chro
        vets<-as.numeric(factor(vals$stats[reg1,4],levels=unique(levs) ))
        bets<-as.numeric(factor(levs,levels=unique(levs) ))
        blata<-vals$data[reg1,]
        if(length(blata)>0)
            t(mapply(
                function(x){
                    reg<-x==vets
                    data<-blata[reg,]
                    le<-sum(c(reg,TRUE))
                    if(le==2)
                        data
                    else if (le>2)
                        colSums(data)/(le-1)
                    else
                        rep(0,6)
                },
                bets)) },
                               split(geneomes$name,geneomes$chrom),levels(geneomes$chrom)))

    normData<-standPCAPrep(data,"colQn")
    pcs<-prcomp(t(normData))
    plotPCs(pcs,c(1,2),normData,Map(as.character,cats)[[1]])

    pc=3
    ord<-order(pcs$rotation[,pc])
    n=15
    ##plot(hclust(dist(rbind(data[ord[1:n],],data[rev(ord)[1:n],]))))
    d<-1-cor(t(rbind(data[ord[1:n],],data[rev(ord)[1:n],])))^2
    rownames(d)<-rbind(as.character(geneomes$name2[ord[1:n]]),as.character(geneomes$name2[rev(ord)[1:n]]))
    colnames(d)<-c(rep(0,n),rep(1,n))
    #heatmap(d)
    d<-as.dist(d)
    plot(hclust(d))
    a<-unique(geneomes[order(pcs$rotation[,1]),"name2"][1:50])
    b<-unique(geneomes[order(pcs$rotation[,1],decreasing = TRUE),"name2"][1:50])
    #plot(seq(length(pcs$rotation[,1])),10^-pcs$rotation[order(pcs$rotation[,1]),1])
    #plot(pcs$rotation[,c(5,1)])
    #plotBox(pcs,5,data,Map(as.character,cats)[[1]])    
    opts<-list(list("awt",c(1,2),c("bottom","bottom"),c(2,2)),
               list("akd",c(1,2),c("bottom","top"),c(2,2)))           
    #opts<-list(list("awt",3,"top",3),list("akd",3,"bottom",3))
           
    a<-lapply(opts,function(x)importantGenes(vals,geneomes,x[[2]],x[[3]],x[[4]]))
    #a<-vals$stats[which((function(x,n,fun,s=sqrt(var(x)),m=mean(x)){fun(x,m,s,n)})(pcs$rotation[,3],1,function(x,m,s,n){x>m+s*n}))]
    #a<-vals$stats[intersect(
    #which((function(x,n,fun,s=sqrt(var(x)),m=mean(x)){fun(x,m,s,n)})(pcs$rotation[,3],1,function(x,m,s,n){x<m-s*n})),
    #which((function(x,n,fun,s=sqrt(var(x)),m=mean(x)){fun(x,m,s,n)})(pcs$rotation[,2],0,function(x,m,s,n){x<m-s*n}))),]
    #cat(unique(unlist(apply(a,1,function(x)findGene(x,geneomes)))),sep="\n")      
    findGene<-function(x,geneome,get="name2"){
        a<-geneome[
                   which(apply(cbind(geneomes$chrom==x[1],as.numeric(geneomes$txStart)<=as.numeric(x[2]),as.numeric(geneomes$txEnd)>=as.numeric(x[3]) ),1,all))
                  ,]
        unlist( Map(as.character,a[get] )[[1]])
    }

    importantGenes<-function(vals,geneome,...){
        test<-ascore(vals$data,...)
        t2<-vals$stats[apply(test,1,all),]
        unique(unlist(apply(t2,1,findGene,geneome)))
    }
    a
}
