hostname<-"/mnt/brand01-00/mbrand_analysis/peaks/october/"
#catfile<-"/home/griffita/Dropbox/UTX-Alex/jan/catagories"
catfile<-"/home/agriffith/Dropbox/UTX-Alex/jan/catagories"
cats<-read.table(catfile)

fileList<-paste(hostname,apply(cats,1,as.character),"/combined_mock_peaks.xls",sep="")
    prefix<-"/mnt/brand01-00/mbrand_analysis/data_sets/"
    suffix<-"_sorted.bed"


macsCutoffs<-c(5,10,20)
for(i in macsCutoffs){
    fnam<-paste("combined_peaks_",i,".bed",sep="")
    fnam2<-paste("combined_heights_",i,".bed",sep="")
    a<-tags(fileList,cats,300,150,i)
    write.table(a,fnam,quote=FALSE,col.names=FALSE,row.names=FALSE)
    rawdata<-apply(cats,1,function(x){paste(prefix,x,"/",x,suffix,sep="")})
    b<-loadBedFile(fnam)
    data<-hg19Sort(b)
    score<-pileUp(data,rawdata,n=22)
    write.table(cbind(data,score),fnam2,quote=FALSE,col.names=FALSE,row.names=FALSE)
    temp<-cor(score)    
    rownames(temp)<-t(cats)
    colnames(temp)<-t(cats)
    png(paste("~/dend-test",i,".png"))
    plot(hclust(dist(temp),method="average"),hang=-1)
    dev.off()}

score<-loadHeightFile("combined_peaks_1.bed")
temp<-cor(score)    
rownames(temp)<-t(cats)
colnames(temp)<-t(cats)
temp["cem_1","cem_2"]

#FastaFile<-"~/Dropbox/UTX-Alex/jan/combined.fasta"


peakFile <- "~/Dropbox/overlap/combined_heights_10.bed"
#catFile <- "catagories"
pcs<-pcaAnalysis(peakFile,catfile)

pc1=1
pc2=3
#png(paste("~/pcs-test-20",pc1,"-",pc2,".png",sep=""))
plotPCs(pcs[[1]]$rotation,cbind(pc1,pc2),pcs[[2]],apply(cats,1,as.character))
#dev.off()


temp<-cor(loadHeightFile("~/Dropbox/overlap/combined_heights_5.bed")$data)
colnames(temp)<-catagories
rownames(temp)<-catagories
temp["cem_1","k562_1"]


       #,c( paste("PC:",as.character(pc1),sep="" ),paste("PC:",as.character(pc2),sep="" )))


data<-loadHeightFile("~/Dropbox/UTX-Alex/jan/combined_heights.bed")$data

# 1+   Eryt/K562  1-
# 1-   T-All      1+
# 3+   ECFC       3-
# 3+5+ Other      3-5-, 3+5- , 3-5+
# 3+7- HSPC       3-7-, 3+,7+, 3-,7+
# 3+7+ MEKA       3-,7-, 3-,7+ 3+,7-
# 7+   Eryt/Meka  7-

reg<-mapply(function(pc,loc)buildRegions(data,pc,loc)[,1],
            list(1,1,3,c(3,5),c(3,7),c(3,7),7),
            list("top","bottom","top",c("top","top"),c("top","bottom"),c("top","bottom"),"top"))

Sequences <- readDNAStringSet(FastaFile, "fasta")
motifFile<-"inst/data/abnormal_normal.pwm"
 motifs<-homerWrapper(Sequences,reg[,1],reg[,2],"~/Masters/mulcal/inst/lib/homer-4.7/cpp/homer2",motifFile)
motifs<-loadPWM(motifFile)

#tags<-read.table("~/Dropbox/UTX-Alex/jan/combined_sorted.bed")[,4]
#catagories<-unique(unlist(lapply(levels(tags),strsplit,"-")))
#values<-loadHeightFile("~/Dropbox/UTX-Alex/jan/combined_heights.bed")
tags<-read.table("~/Dropbox/overlap/combined_peaks_5.bed")[,4]
#catagories<-unique(unlist(lapply(levels(tags),strsplit,"-")))

catagories<-unlist(lapply(cats,as.character))
values<-loadHeightFile("~/Dropbox/overlap/combined_heights_5.bed")
data<-values$data
info<-values$info
pcs<-as.matrix(mulcal::ascore(data,c(1:3)))

cnames<-c("tall_p1",
"tall_p2_1",
"tall_p2_2",
"tall_p3_1",
"tall_p3_2",
"jurk_sandar_1",
"jurk_sandar",
"jurk",
"rpmi_1",
"rpmi_2",
"cem_1",
"cem_2",
"ecfc-tsa",
"meka",
"cd133",
"cd34",
"cd34_new",
"eryt",
"eryt_f",
"eryt_a",
"k562_1",
"k562_2")

cnames<-list(list("tall_p1","tall_p2_1","tall_p2_2","tall_p3_1","tall_p3_2","jurk_sandar_1","jurk_sandar","jurk","rpmi_1","rpmi_2","cem_1","cem_2"),
list("ecfc-tsa"),list("meka"),list("cd133","cd34","cd34_new"),list("eryt","eryt_f","eryt_a","k562_1","k562_2"))

cnames<-list(list("tall_p1","tall_p2_1|tall_p2_2","tall_p3_1|tall_p3_2","jurk_sandar_1|jurk_sandar|jurk","rpmi_1|rpmi_2","cem_1|cem_2"),
list("ecfc-tsa"),list("meka"),list("cd133","cd34|cd34_new"),list("eryt|eryt_f|eryt_a","k562_1|k562_2"))

cnames<-list(list("tall_p1|tall_p2_1|tall_p2_2","tall_p3_1|tall_p3_2|jurk_sandar_1|jurk_sandar|jurk|rpmi_1|rpmi_2|cem_1|cem_2"),
list("ecfc-tsa"),list("meka"),list("cd133|cd34|cd34_new"),list("eryt|eryt_f|eryt_a|k562_1|k562_2"))

names<-c("t-all","endotheilel","meka","stem","eryt")

strplit("tall_p1 tall_p2 tall_p2 tall_p3 tall_p3 jurk jurk jurk rpmi rpmi cem cem ecfc meka cd133 cd34 cd34 eryt eryt eryt k562 k562",sep=" ")

keyframe<-data.frame(s="stem",t="t-all",e="Endotheiliel",d="Eryt",m="meka")

swCat<-as.character(unlist(sapply(strsplit("t,t,t,t,t,t,t,t,t,t,t,t,e,m,s,s,s,d,d,d,d,d",",")[[1]],function(x) keyframe[x])))

names(swCat)<-cnames

#swCat<-c("normal","normal","abnormal","abnormal","abnormal","abnormal","abnormal","abnormal","meka","ecfc","stem","stem")
#names(swCat)<-c("eryt","k562","jurk","cem","rpmi","tall_p1","tall_p2","tall_p3","meka","ecfc","cd133","cd34")
tSwCat<-as.character(unlist(lapply(catagories,function(x) swCat[x])))

png("~/contrib-10-3.png")
stackedContrib2(pcs[,3],tags,names,cnames,steps=100,n=7,inord=FALSE,sum=TRUE)
dev.off()
