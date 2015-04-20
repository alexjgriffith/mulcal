#!/usr/bin/env R
#
# This file is part of peakAnalysis
# , http://github.com/alexjgriffith/alpha-score/, 
# and is Copyright (C) University of Ottawa, 2015. It is Licensed under 
# the three-clause BSD License; see LICENSE.txt.
# Author : Alexander Griffith
# Contact: griffitaj@gmail.com


#' member
#' @export
member<-function(posibList,mem,test="default"){
    fun<-switch(test,default=function(x,y){all(y %in% x)})
    sapply(posibList,fun,mem)}

#' Member Wrapper
#' Ment for bed data where the member is found in the fourth column
#' @export
memberWrapper<-function(bedData,lis ,n=4,sep="-"){
    member(strsplit(as.character(bedData[,n]),sep),lis)}

#' @export
cdf<-function(m,range){
    x<-abs(m)    
    out<-c(x[1])
    for(i in seq(length(x)-1))
        out<-c(out,out[i]+x[i+1])
    out     
}

#' @export
pdf<-function(x){
 l<-length(x)
 x[2:(l-0)]-x[1:(l-1)]
}

#' @export
windowing.data.frame<-function(x,fun,window,sliding=FALSE) {
        n<-dim(x)[1]
        if(sliding)
            unlist(lapply(seq(floor(n/window)),function(i) do.call(fun,list(x[((i-1)*window+1):(i*window),]))))
        else
            unlist(lapply(seq(n-window),function(i) do.call(fun,list(x[i:(i+window),])))) 
}

#' @export
windowing.vector<-function(x,fun,window,sliding=FALSE){
        n<-length(x)
        if(sliding)
            unlist(lapply(seq(floor(n/window)),function(i) do.call(fun,list(x[((i-1)*window+1):(i*window)]))))
        else
            unlist(lapply(seq(n-window),function(i) do.call(fun,list(x[i:(i+window)]))))       
        }

#' @export
windowing<-function(x,...){
    if(is.null(attr(x,"class"))){
        args<-list(x,...)
        if(is.numeric(x))
            do.call(windowing.vector,args)
        else if(is.matrix(x))
            do.call(windowing.data.frame,args)
        }
    else UseMethod("windowing",x)
}



# aux.r -> utils.r
#' @export
collapse<-function(...,sep=""){paste(...,sep=sep)}
#' @export
lcollapse<-function(x,sep=""){br<-"";for(i in x){br<-paste(br,i,sep=sep)};br}
#' @export
collect<-function(x,fn,...){lcollapse(sapply(x,fn,...))}
#' @export
modulous<-function(x,m)
    {t1<-floor(x/m)
     (x-t1*m)}

#' @export
is.even<-function(x){
    modulous(x,2)==0
}

#' @export
is.odd<-function(x){
    modulous(x,2)!=0
}

#' @export
is.prime<-function(x){
    length(which( modulous(x,seq(1:ceiling(x/2)))==0))==1
}



#' Subset nestedList
#'
#' subsets a nested list
#' 
#' @export
subset.nestedList<-function(h,locs,fun,...){
    lapply(combinations(locs),fun,h,...)           
}

#' Combinations
#'
#' Generates combinations along a set of lists
#' 
#' @export
combinations<-function(vars){
    v<-vars[[1]]
    for(i in vars[-1]){
        t<-c()
        for(k in v)
            for(j in i){
                t<-c(t,list(c(k,j)))
            }
        v<-t
    }
    return(v)
}
