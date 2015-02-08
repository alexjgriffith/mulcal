#!/data/binaries/R-3.1.0/bin/Rscript
library('getopt')
source("gene-assoc.r")

spec = matrix(c(
    'fileLocation','f', 1,"character"
    ),byrow=TRUE,ncol=4)
args=getopt(spec)

fileLocation=args$fileLocation 
bedData<-read.delim(fileLocation,header=0)
geneList<-read.delim("../info/hg19.RefSeqGenes.csv")
genes<-geneAssociation(bedData,geneList,  c(50000,0,0,0),n=15)

data<-cbind(bedData,sapply(genes, function(x){x[1]}))

write.table(data,stdout(),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")

