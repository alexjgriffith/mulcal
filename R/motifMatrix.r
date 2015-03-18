# The portions of saving, motif hist and sigMatrix required to
# save the data in chunks
# motif_1-motif_2-enviroment
# option to save every x motif_2

#' @export
partialSigMatrix<-function(Sequences,mots,reg,ranges,
                           fields=unlist(Map(as.character,1:dim(reg)[2])),
                           addmotifs=c(),
                           n=rep((length(mots)+length(addmotifs)),2),
                           m=dim(reg)[2],filename="test"){
    mots<-c(mots[1:(n[2]-length(addmotifs))],addmotifs)
    mList<-unlist(lapply(c(mots,addmotifs),IUPACtoBase))
    cList<-unlist(lapply(lapply(c(mots,addmotifs),IUPACtoBase),compliment))
    locationsM<-lapply(mList,grep,Sequences)
    locationsC<-lapply(cList,grep,Sequences)
    # Using for loops in place of the apply family
    # in order to provide options to dynamicly
    # modify the printed output.
    # An alternative method will be implimented which
    # will not permit for such flexabiluty but which
    # may be parralized easily.
    temp<-0
    h<-list()
    bp<-4
    m1<-c()
    m2<-c()
    append=FALSE
    for(k in seq(n[1])){
        print(mots[k])
        for(i in seq(n[2])){
         h<-append(h,lapply(seq(m),function(j)motifHist(Sequences,mList,cList,locationsM,locationsC,i,k,reg[,j])))
         m1<-c(m1,mots[k])
         m2<-c(m2,mots[i])
         if(! temp<bp){
             print("saving")             
             saveData(filename,mode="perMot2",append=append,h,m1,m2,fields,ranges)
             append=TRUE
             m1<-c()
             h<-list()
             m2<-c()
             temp<-0
         }
         else
             temp<-temp+1
     }
    }
    if(temp!=0){
        print("saving")             
        saveData(filename,mode="perMot2",append=append,h,m1,m2,fields,ranges)
    }
    h<-read.table(filename,header=TRUE,check.names = FALSE)
    h
}

#' @export
anotatePaste<-function(anot,value,start="#$"){
    paste(c(paste(c(start,anot),collapse=""),value),collapse=" ")
}

#' @export
flattenHeights<-function(h,range=range(unlist(h))){
       x<-range#do.call(seq,as.list(range))
       do.call(rbind,lapply(h,function(y1) do.call(rbind,lapply(y1,function(y) do.call(rbind,lapply(y,getHeights,x))))))}
       #do.call(rbind,lapply(h,function(y1) do.call(rbind,lapply(y1,function(y) do.call(rbind,lapply(y,"+",x))))))}

#' @export
saveData<-function(filename,mode=FALSE,append=FALSE,...){
    generateHeader<-function(mots,fields,ranges){
        c(anotatePaste("fields",fields),
                   anotatePaste("motifs",mots),
                   anotatePaste("range",ranges))
    }
    generateData<-function(h,ranges=range(unlist(h))){
        a<-flattenHeights(h,ranges)
        b<-apply(a,1,paste,collapse=" ")
        b
    }
    generateAll<-function(h,mots,fields,ranges){
        c(generateHeader(mots,fields,ranges),generateData(h))
    }
    perMot2<-function(h,m1,m2,fields,ranges){
        colNames<-do.call(seq,as.list(ranges))
        rowNames<-unlist(lapply(combinations(list(paste(m1,m2,sep="-"),fields)),paste,collapse="-"))        
        output<-do.call(rbind,lapply(h,getHeights,ranges))
        colnames(output)<-colNames
        rownames(output)<-rowNames
        output
    }
    perMot1<-function(h,m1,m2,fields,ranges){
        colNames<-do.call(seq,as.list(ranges))
        rowNames<-unlist(lapply(combinations(list(m1,m2,fields)),paste,collapse="-"))
        output<-flattenHeights(h,ranges)
        colnames(output)<-colNames
        rownames(output)<-rowNames
        output
    }
    f<-list(row.names=FALSE,
            col.names=FALSE,
            quote=FALSE)
    n<-list(row.names=TRUE,
            col.names=(! append),
            quote=FALSE)    
    output<-switch(mode,
                   body=append(list(data=generateData(...)),f),
                   head=append(list(data=generateHeader(...)),f),
                   all=append(list(data=generateAll(...)),f),
                   perMot=append(list(data=perMot1(...)),n),
                   perMot2=append(list(data=perMot2(...)),n))
    with(output,
         write.table(data,filename,
                     row.names=row.names,
                     col.names=col.names,
                     quote=quote,
                     append=append)
         )
}

#' @export
loadData<-function(filename){    
    a<-unlist(as.character(read.delim(filename,header=FALSE,comment.char="")[,1]))
    coms<-strsplit(a[grep("#",a)]," ")
    keys<-lapply(coms,"[[",1)
    values<-lapply(coms,function(com){
        l<-length(com)
        com[2:l]})               
    x<-do.call(seq,as.list(as.numeric(values[keys=="#$range"][[1]]) ))
    motifs=as.list(values[keys=="#$motifs"][[1]])
    fields=as.list(values[keys=="#$fields"][[1]])
    data<-read.table(filename)
    colnames(data)<-x
    rownames(data)<-unlist(lapply(motifs,function(m1)lapply(motifs,function(m2)lapply(fields,function(f)paste(m1,m2,f,sep="-") ))))
    data
}

#' @export
buildRegions<-function(data,pc,location,directory="~/Masters/mulcal/inst/data/",title=paste(directory,paste(apply(cbind(pc,location),1,function(x) paste(x,collapse="_") ),collapse="_"),sep="")){
    testCombs<-function(x,data){
        data<-as.matrix(data)
        m<-do.call(ascore,append(list(data),lapply(x[2:4],unlist)))
        xy<-dim(x[[5]])
        crossFun<-function(data,pos,fun){
            if(length(pos)>1)
                apply(data[,pos],1,fun)
            else
                data[,pos]
        }
        fg<-crossFun(m,x[[5]][1,],"all")
        if(xy[1]>2){
            bg1<-apply(x[[5]][2:xy[1],],1,function(x)crossFun(m,x,"all"))
            bg<-crossFun(bg1,seq(xy[1]-1),"any")}
        else{
            bg<-crossFun(m,x[[5]][2,],"all")
        }
        cbind(fg,bg)
    }    
    fg<-mapply(function(x,loc){
        y<-1
        if(loc=="bottom"){y=2}
        (x-1)*2+y
    },seq(length(pc)),location)
    temp<-do.call(rbind,combinations(lapply(1:length(pc),function(x) list((x-1)*2+1,(x-1)*2+2) )))
    loc<-apply( temp,1,function(x) all(x==fg ))
    bg<-matrix(temp[! loc],ncol=length(pc))
    l<-list(title,unlist(lapply(pc,function(x)rep(x,2))),rep(c("top","bottom"),length(pc)), rep(1,length(pc)*2),rbind(fg,bg))
    reg<-do.call(testCombs,list(l,list(data)))
    reg
}

#' motifsFromPCA
#'
#' @examples
#' l<-mapply(function(a,b,c){
#'    motifsFromPCA(Sequences,data,a,b,title=c)
#' },
#' list(1,1,3,c(3,5),c(3,7),c(3,7),7),
#' list("top","bottom","top",c("top","top"),c("top","bottom"),c("top","bottom"),"top"),
#' paste("~/Masters/mulcal/inst/data/",list("normal","abnormal","ecfc","other","hspc","meka","diff"),".pwm",sep=""))
#' 
#'l<-motifsFromPCA(Sequences ,data,c(1,2,3),c("top","top","bottom"))
#' @export
motifsFromPCA<-function(Sequences,data,pc,location,directory="~/Masters/mulcal/inst/data/",title=paste(directory,paste(apply(cbind(pc,location),1,function(x) paste(x,collapse="_") ),collapse="_"),sep="")){
    reg<-buildRegions(data,pc,location)
    motifs<-lapply(6:8,function(x)
        homerWrapper(Sequences,reg[,1],reg[,2],"~/Masters/mulcal/inst/lib/homer-4.7/cpp/homer2",opts=paste("-S 25 -len",as.character(x))))
    mots<-unlist(lapply(motifs,function(x) gsub(">","",x[,1])))
    fg<-reg[,1]
    all<-rep(TRUE,length(fg))
    fgLength<-length(which(fg)) 
    x<-lapply(unlist(lapply(mots,IUPACtoBase)), function(y)lapply(list(fg,all), function(reg) length(grep(y, Sequences[reg]))))    
    allMotifs<-data.frame(motifs=mots,indep=unlist(lapply(x,function(x)x[[1]]/x[[2]]*length(fg)/fgLength)),num=unlist(lapply(x,function(x)x[[1]])))    
    topMotifs<-apply(rbind(allMotifs$num/length(which(fg))>0.05,allMotifs$indep>1.1),2,all)
    saveMotifFile(cbind(paste(">",as.character(allMotifs[topMotifs,"motifs"]),sep=""),do.call(c,lapply(motifs,function(x) x[,2]))[topMotifs]),title)
    loadPWM(title)
    }

#' @export
saveMotifFile<-function(l,filename){
    output<-unlist(apply(l,1,function(x) {c(x[1][[1]],apply(t(x[2][[1]]),1,paste,collapse="\t"))}))
    write.table(output,filename,row.names=FALSE,col.names=FALSE,quote=FALSE)}

#' @export
breakDownName<-function(t0){
    strsplit(as.character(rownames(t0)),"-")[[1]]
}

#' @export
zeroWithinWidth<-function(t0,width,x){
    ma<-order(t0,decreasing = TRUE)[1]+x[1]
    while( all(ma>width[1],ma<width[2]) ){                                            
        t0[order(t0)]<-0
        ma<-order(t0,decreasing = TRUE)[1]+x[1]
    }
    t0
}

#' @export
quick2motifs<-function(motifA,motifB,motifs,names="None",n=c(1,2,3)){
    L1<- which(motifA==motifs)
    L2<- which(motifB==motifs)
    print(L1)
    print(L2)
    histWrapper(h,motifs,names="None",L1,L2,n)}

#' @export
sigMatrix<-function(Sequences,mots,reg,addmotifs=c(),
           n=(length(mots)+length(addmotifs)),
           m=dim(reg)[2]){
    mots<-c(mots[1:(n-length(addmotifs))],addmotifs)
    mList<-unlist(lapply(c(mots,addmotifs),IUPACtoBase))
    cList<-unlist(lapply(lapply(c(mots,addmotifs),IUPACtoBase),compliment))
    locationsM<-lapply(mList,grep,Sequences)
    locationsC<-lapply(cList,grep,Sequences)
    h<-lapply(seq(n),function(k)lapply(seq(n),function(i)lapply(seq(m),function(j)motifHist(Sequences,mList,cList,locationsM,locationsC,i,k,reg[,j]))))
    #motifGetScore(m,h,n)
    h
}
