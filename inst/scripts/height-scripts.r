mainProf<-function(file="/home/agriffith/Dropbox/UTX-Alex/jan/combined_sorted.bed",
                   rawdata=rep(file,22),
                   n=0){
    Rprof("temp.prof")
    data<-hg19Sort(loadBedFile(file))
    score<-pileUp(data,rawdata,n=n)
    print(str(score))
    Rprof(NULL)
    summaryRprof("temp.prof")}


testHeights<-function(){
    file<-"/home/griffita/Dropbox/UTX-Alex/jan/combined_sorted.bed"

    cats<-read.table("/home/griffita/Dropbox/UTX-Alex/jan/catagories")
    prefix<-"/mnt/brand01-00/mbrand_analysis/data_sets/"
    suffix<-"_sorted.bed"
    rawdata<-apply(cats,1,function(x){paste(prefix,x,"/",x,suffix,sep="")})

    data<-hg19Sort(loadBedFile(file))
    score<-pileUp(data,rawdata,n=22)
    temp<-cor(score)
    rownames(temp)<-t(cats)
    colnames(temp)<-t(cats)
    pdf("test.pdf")
    plot(hclust(dist(temp)),hang=-1)
    dev.off()}
