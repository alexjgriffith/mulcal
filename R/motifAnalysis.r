#library(Biostrings)

#' @export
motifTest<-function(fastaFile="~/Dropbox/UTX-Alex/jan/combined.fasta",
                    heightFile="~/Dropbox/UTX-Alex/jan/combined_heights.bed",
                    motifs=c("CGGGCG","GGATAA")){
    values<-loadHeightFile(heightFile)
    data<-values$data
    reg<-ascore(data,1,"top",3)
    test<-readDNAStringSet(fastaFile,use.names = TRUE)[reg]
    ma<-as.matrix(test)
    a<-lapply(motifs,peaksWithMotifs,test)
    locations<-intersect(a[[1]],a[[2]])
    b<-lapply(motifs,motifLocations , ma,locations)
    h<-getDistance(b[[1]],b[[2]])
    sqrt(var(h))
    hist(h,breaks=20)}


peaksWithMotifs<-function(motif,data){
    l<-nchar(data[1])
    lm<-nchar(motif)
    sort(unlist(lapply(seq(l-lm),function(i){which(motif==subseq(data,i,width=lm))})))}

motifLocations<-function(motif,ma,locations){
    seq<-unlist(strsplit(motif,""))
    apply(ma[locations,],1,function(vect){
        l<-length(vect)
        l2<-length(seq)
        v1<-as.matrix(sapply(seq(l2),function(x){vect[x:(l-l2+x)]}))
        which(apply(v1,1,function(v1){all(v1==seq)}))})}

getDistance<-function(x,y)
    as.numeric(mapply(function(x,y){
        temp<-outer(x, y,"-")
        temp<-temp[upper.tri(x=temp,diag=TRUE)]
        temp[which.min(abs(temp))]},x,y))



