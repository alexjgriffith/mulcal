library(mulcal)
library(Biostrings)

heightFile<-"~/Dropbox/UTX-Alex/jan/combined_heights.bed"
fileLocation<-"~/Dropbox/UTX-Alex/jan/"
geneList<-read.delim(paste(fileLocation,"hg19.RefSeqGenes.csv",sep=""))
bedData<-read.delim(paste(fileLocation,"combined_sorted.bed",sep=""),header=0)
data<-loadHeightFile(heightFile)$data



## Generate a gene list using Stanford great
chrom<-as.character(geneList$chrom)
tss<-as.numeric(geneList$txStart)
a<-genomicRegions(chrom,tss,5000,1000,50000) # takes about 3 min


inlist<-paste( as.character(geneList$name),
      as.character(geneList$name2), sep="\t")


for(n in c(1,3,6))
{
reg<-mapply(function(pc,loc)buildRegions(data,pc,loc,ns=n)[,1],
            list(1,1,3,c(3,5),c(3,7),c(3,7),7),
            list("top","bottom","top",c("top","top"),c("top","bottom"),c("top","bottom"),"top"))
filenames<-lapply(cbind("erythroid","t-all","ecfc","other","hspc","meka","diff"),paste, "-refseq-", n, "-genes.txt",sep="")
print(filenames)
for (i in seq(7)){
   genes<-peakGeneRegions(bedData[reg[,i],],a,inlist)
   write.table(genes,filenames[[i]],quote=FALSE,col.names=FALSE,row.names=FALSE)}
}

## Generate a gene list based on proximity
## this method is gene centric and usefull for finding
## the binding ocupancy frequency of each cell type
a<-geneAssocFlipAux(geneList,bedData,c(50000,500,0,0))
names<-mapply(function(a,b) paste(a,b,sep="-"), as.character(geneList$name),as.character(geneList$name2)) 
for (i in seq(7)){
   print(filenames[[i]])
   filename<-filenames[[i]]
   sequence<-reg[,i]
   b<-lapply(a,function(x,y) list(x[[1]],x[[2]][which(x[[2]]%in%y)]) ,which(sequence))
   c<-nearestGene(b,sequence,as.numeric(geneList$txStart),names,bedData)
   write.table(unique(unlist(c)),filename ,col.names = FALSE,row.names=FALSE,quote=FALSE)
}

## most basic (and slow) gene assoiation
## has the possibility of parallel alignement
genes<-geneAssociation(bedData,geneList,  c(50000,0,0,0))
