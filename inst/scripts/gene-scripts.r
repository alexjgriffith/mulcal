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


#############################################################
### Using the GREAT2.0 list of genes


heightFile<-"~/Dropbox/UTX-Alex/jan/combined_heights.bed"
fileLocation<-"~/Dropbox/UTX-Alex/jan/"
geneList<-read.delim(paste("inst/data/","hg19.great2.csv",sep=""))
#bedData<-read.delim(paste(fileLocation,"combined_sorted.bed",sep=""),header=0)
v<-loadHeightFile(heightFile)
data<-v$data
bedData<-data.frame(chro=v$stats[,1],start=as.numeric(v$stats[,2]), end= as.numeric(v$stats[,3]))


## Generate a gene list using Stanford great
chrom<-as.character(geneList$chrom)
tss<-as.numeric(geneList$txStart)
strand<-rep(1,length(tss))
strand[geneList$strand=="-"]=-1
inlist<-as.character(geneList$name2)

## Find Subsets of the unified peak set which are biologicaly relevent
n<-1
i<-1
reg<-mapply(function(pc,loc)buildRegions(data,pc,loc,ns=n)[,1],
            list(1,1,3,c(3,5),c(3,7),c(3,7),7),
            list("top","bottom","top",c("top","top"),c("top","top"),c("top","bottom"),"top"))

## Generate genomic Regions of Interest
a<-genomeDomains(chrom,tss,strand,5000,1000,1000000) # takes about 3 min
b<-genomicRegions(a,chrom)

## Find Genes of Interest
genes<-peakGeneRegions(bedData[reg[,2],],b,inlist)



c1<-as.data.frame(do.call(rbind,lapply(b,function(x) if(x[1]=="chr1")as.numeric(x[c(2,3)]))))
colnames(c1)<-c("start","end")
c2<-bedData[bedData[,1]=="chr1",c(2,3)]
colnames(c2)<-c("start","end")

outer((c1$end+c1$start)/2,c2$start,function(x,y) x+y)


tests<-function(a,b){
max(unlist(lapply(b,function(x) -1* do.call("-",as.list(as.numeric(x[c(2,3)]))))))

max(a[,2]-a[,1])
plot(genes,seq(1:length(genes)),type="l")
}


v<-loadHeightFile("~/Dropbox/UTX-Alex/jan/combined_heights.bed")
data<-v$data
stats<-v$stats
formatOut<-apply(sub("   ","",stats[,1:3]), 1,paste,collapse="\t")

FastaFile<-"~/Dropbox/UTX-Alex/jan/combined.fasta"
Sequences <- readDNAStringSet(FastaFile, "fasta")
smad<-rep(FALSE,length(Sequences))
smad[grep(IUPACtoBase("AGNCAGAC"),Sequences)]=TRUE
smad[grep(compliment(IUPACtoBase("AGNCAGAC")),Sequences)]=TRUE
tallSmad<-(smad&reg[,2])
erytSmad<-(smad&reg[,1])
hspcSmad<-(smad&reg[,5])

write.table(formatOut[smad&reg[,1]],"tall-smad_peaks.bed",quote=FALSE,row.names=FALSE,col.names=FALSE)


great<-unlist(levels(read.table("inst/data/tallSmadGreatGenes.txt")$V1))

tallSmadGenes<-peakGeneRegions(bedData[tallSmad,],b,inlist)
erytSmadGenes<-peakGeneRegions(bedData[erytSmad,],b,inlist)
hspcSmadGenes<-peakGeneRegions(bedData[hspcSmad,],b,inlist)

buildCompare<-function(x,inlist){
    ret<-rep(FALSE,length(inlist))
    ts<-unlist(lapply(x,function(x) which(inlist==x)))
    ret[ts]=TRUE
    ret
}

tallLG<-buildCompare(tallSmadGenes,inlist)
erytLG<-buildCompare(erytSmadGenes,inlist)
hspcLG<-buildCompare(hspcSmadGenes,inlist)

write.table(peakGeneRegions(bedData[tallSmad,],b,inlist),"~/Dropbox/UTX-Alex/smad/tall-smad-code-genes.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

write.table(peakGeneRegions(bedData[erytSmad,],b,inlist),"~/Dropbox/UTX-Alex/smad/erty-smad-code-genes.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

write.table(peakGeneRegions(bedData[hspcSmad,],b,inlist),"~/Dropbox/UTX-Alex/smad/hspc-smad-code-genes.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

write.table(peakGeneRegions(bedData[(hspcSmad | tallSmad),],b,inlist),"~/Dropbox/UTX-Alex/smad/hspc_union_tall-smad-code-genes.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)


write.table(inlist[union(which( tallLG), which( erytLG))],"~/Dropbox/UTX-Alex/smad/erty_union_tall-smad-code-genes.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

great[which(unlist(lapply(great,function(i)  (i %in% inlist &  ! i %in% tallSmadGenes))))]

tallSmadGenes[which(unlist(lapply(tallSmadGenes,function(i)  (i %in% inlist &  ! i %in% great))))]

length(great[unlist(lapply(great,function(i)    i %in% tallSmadGenes))])/
length(great[unlist(lapply(great,function(i)    i %in% inlist))])



write.table(formatOut[smad],"smad_peaks.bed",quote=FALSE,row.names=FALSE,col.names=FALSE)


gata3=which(as.character(bedData[,1])=="chr10" &  as.numeric(bedData[,2])==8103372 )

bedData[gata3,]

which(unlist(lapply(b,function(x) x[1]=="chr10" & as.numeric(x[2])<8103372 & as.numeric(x[3])>8103372 )))


which(unlist(lapply(b,function(x) 1838 %in% x[4:length(x)])))
