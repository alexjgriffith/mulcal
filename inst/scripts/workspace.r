data<-loadHeightFile("~/Dropbox/UTX-Alex/jan/combined_heights.bed")$data
reg<-as.matrix(ascore(data,c(1,1,1),c("top","bottom","middle"),c(1,1,1)))
FastaFile<-"~/Dropbox/UTX-Alex/jan/combined.fasta"
Sequences <- readDNAStringSet(FastaFile, "fasta")
motifs<-rbind(loadPWM("inst/data/abnormal_normal.pwm"),loadPWM("inst/data/normal_abnormal.pwm"),loadPWM("inst/data/jaspar_motifs.txt","jaspar"))
mots<-unlist(lapply(unlist(motifs[,1]),function(x)gsub(">","",x)))
addmotifs<-c("CANNGT","TACCTC")

n=10
m=3
mots<-c(mots[1:(n-length(addmotifs)]),addmotifs)
mList<-unlist(lapply(c(mots,addmotifs),IUPACtoBase))
cList<-unlist(lapply(lapply(c(mots,addmotifs),IUPACtoBase),compliment))
locationsM<-lapply(mList,grep,Sequences)
locationsC<-lapply(cList,grep,Sequences)

h<-lapply(seq(n),function(k)lapply(seq(n),function(i)lapply(seq(m),function(j)motifHist(Sequences,mList,cList,locationsM,locationsC,i,k,reg[,j]))))



scores<-lapply(seq(m),
       function(n)
           matrix(unlist(lapply(h,function(x) lapply(x,function(y)motifScoreFunction(y[[n]]) ))),nrow=10))
