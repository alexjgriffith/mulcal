##########################################

# 1+   Eryt/K562  1-
# 1-   T-All      1+
# 3+   ECFC       3-
# 3+5+ Other      3-5-, 3+5- , 3-5+
# 3+7- HSPC       3-7-, 3+,7+, 3-,7+
# 3+7+ MEKA       3-,7-, 3-,7+ 3+,7-
# 7+   Eryt/Meka  7-


library(mulcal)
library(Biostrings)
library(ggplot2)

denovoMotifs<-paste("~/Masters/mulcal/inst/data/",list("normal","abnormal","ecfc","other","hspc","meka","diff"),".pwm",sep="")
FastaFile<-"~/Dropbox/UTX-Alex/jan/combined.fasta"
heightFile<-"~/Dropbox/UTX-Alex/jan/combined_heights.bed"
jasparFile<-"inst/data/jaspar_motifs.txt"


denovoMotifs<-paste("~/prefered/",list("normal","abnormal","ecfc","other","hspc","meka","diff"),".pwm",sep="")
FastaFile<-"~/prefered/combined.fasta"
heightFile<-"~/prefered/combined_heights.bed"
jasparFile<-"~/prefered/jaspar_motifs.txt"



data<-loadHeightFile(heightFile)$data
reg<-mapply(function(pc,loc)buildRegions(data,pc,loc)[,1],
            list(1,1,3,c(3,5),c(3,7),c(3,7),7),
            list("top","bottom","top",c("top","top"),c("top","bottom"),c("top","bottom"),"top"))

Sequences <- readDNAStringSet(FastaFile, "fasta")

len<-length(l)

### jaspar specific
jasparNames<-unlist(lapply(strsplit(grep(">",Map(as.character,read.csv(jasparFile,header=FALSE))[[1]],value=TRUE)," "),"[",2))
motifs<-loadPWM(jasparFile,"jaspar")
strseq<-unlist(Map(function(i)mulcal::PWMtoCons(motifs[i,2][[1]]),seq(dim(motifs)[1])))
jasparVals<-data.frame(motif=strseq,names=jasparNames)

addmotifs<-c("CANNTG","TACCTC")#,"GATAAG")
motifs<-loadPWM(jasparFile,"jaspar")
motifs[,1]<-apply(motifs,1,function(x) PWMtoCons(x[2][[1]]))
motifs<-do.call(rbind,append(Map(loadPWM,denovoMotifs),list(motifs)))
       
mots<-unique(c(unlist(lapply(unlist(motifs[,1]),function(x)gsub(">","",x))),addmotifs))

#h1<-partialSigMatrix(Sequences,mots,reg,c(-300,300))
#h<-sigMatrix(Sequences,mots,reg,addmotifs,10)

#######################################################

#n=length()
mots<-c(mots[1:(n-length(addmotifs))],addmotifs)
mList<-unlist(lapply(mots,IUPACtoBase))
cList<-unlist(lapply(lapply(mots,IUPACtoBase),compliment))
locationsM<-lapply(mList,grep,Sequences)
locationsC<-lapply(cList,grep,Sequences)

bM<-lapply(mList,function(x) lapply(gregexpr(x, Sequences),as.numeric))
bC<-lapply(cList,function(x) lapply(gregexpr(x, Sequences),as.numeric))


opt<-t(combn(n,2))



m1<-"GATAAG"
m2<-"CANNTG"

viewMotifs<-function(m1,m2,c,mots,bM,bC){


    
}

start<-proc.time()
co<-lapply(seq(dim(opt)[1]),function(x)compareModels(do.call(projectDiff,as.list(opt[x,]))))
proc.time()-start

coM<-do.call(rbind,co)

a<-do.call(cbind,lapply(c(mean,var),function(x)apply(coM,2,x)))
x<-seq(800,1000)/1000
y<-lapply(seq(3), function(i) unlist(lapply(x,dnorm,a[i,1],a[i,2])))
plot(x,y[[1]],type="l",col=2,ylim=c(0,max(unlist(y))))
lines(x,y[[2]],type="l",col=3)
lines(x,y[[3]],type="l",col=4)


b<-c()
for(i in seq(20)){
    x<-data.frame(y=projectDiff(opt[i,1],opt[i,2]),x=seq(-300,300))
    b<-c(b,max(-log10(windowing(function(x) dcauchy(x$y[length(x$y)],median(x$y),interq(x$y)) , 100)(rev(x)))))
}
