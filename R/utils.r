#!/usr/bin/env R
#
# This file is part of peakAnalysis
# , http://github.com/alexjgriffith/alpha-score/, 
# and is Copyright (C) University of Ottawa, 2015. It is Licensed under 
# the three-clause BSD License; see LICENSE.txt.
# Author : Alexander Griffith
# Contact: griffitaj@gmail.com

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
