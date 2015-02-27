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
