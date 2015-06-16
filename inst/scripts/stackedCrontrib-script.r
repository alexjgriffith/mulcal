tags<-read.table("~/Dropbox/UTX-Alex/jan/combined_sorted.bed")[,4]
catagories<-unique(unlist(lapply(levels(tags),strsplit,"-")))
values<-loadHeightFile("~/Dropbox/UTX-Alex/jan/combined_heights.bed")
data<-values$data
info<-values$info
pcs<-as.matrix(mulcal::ascore(data,c(1:10)))


swCat<-c("normal","normal","abnormal","abnormal","abnormal","abnormal","abnormal","abnormal","meka","ecfc","stem","stem")
#swCat<-c("Erythroid","K562","Jurkat","CEM","RPMI","Prima5(Brand)","Prima2","Prima5","Megakaryocite","ECFC","CD133","CD34")
names(swCat)<-c("eryt","k562","jurk","cem","rpmi","tall_p1","tall_p2","tall_p3","meka","ecfc","cd133","cd34")
tSwCat<-as.character(unlist(lapply(catagories,function(x) swCat[x])))

stackedContrib(pcs[,5],tags,catagories,swCat,steps=41,n=3,inord=c(1,3,5,2,4),sum=FALSE)

saveSC<-function(n,name){
png(name)
stackedContrib(pcs[,n],tags,catagories,tSwCat,steps=41,n=3,inord=c(1,3,5,2,4),sum=TRUE)
dev.off()
}

pos<-c(1,3,5,7)
filenames<-paste("~/Dropbox/pcs-",pos,"-contrib-s3.png",sep="")

for(i in seq(length(pos))){
    n<-pos[i]
name<-filenames[i]

    png(name)
    stackedContrib(pcs[,n],tags,catagories,tSwCat,steps=41,n=3,inord=c(1,3,5,2,4),sum=FALSE)
    dev.off()    

}



png("~/Dropbox/pcs-1-contrib-s3.png" )
stackedContrib(pcs[,1],tags,catagories,tSwCat,steps=41,n=3,inord=c(1,3,5,2,4),sum=TRUE)
dev.off()

png("~/Dropbox/pcs-1-contrib-stacked-s3.png" )
stackedContrib(pcs[,1],tags,catagories,tSwCat,steps=41,n=3,inord=c(1,3,5,2,4),sum=FALSE)
dev.off()

png("~/Dropbox/pcs-3-contrib-s3.png" )
stackedContrib(pcs[,3],tags,catagories,tSwCat,steps=41,n=3,inord=c(1,3,5,2,4),sum=TRUE)
dev.off()

png("~/Dropbox/pcs-3-contrib-stacked-s3.png" )
stackedContrib(pcs[,3],tags,catagories,tSwCat,steps=41,n=3,inord=c(1,3,5,2,4),sum=FALSE)
dev.off()

n=2000
plot(matrix(rnorm(n,0,1), ncol=2)%*%matrix(c(1,1,2,1),ncol=2),c="r")

