fun<-function(x){-log10(do.call(max,normTest(x)))}

windowing<-function(fun,window){
    function(x){
        n<-length(x$y)
        n2<-n-window
        unlist(lapply(seq(n2),function(i) do.call(fun,list(x[i:(i+window),]))))
    }
}

nonSlidingWindowing<-function(fun,window){
    function(x){
        n<-length(x$y)        
        n2<-floor(n/window)
        unlist(lapply(seq(n2),function(i) do.call(fun,list(x[((i-1)*window+1):(i*window),]))))
    }
}

twoDNorm<-function(translation=c(0,0),scale=c(1,1),rpheta=0){
    rv<-function(n=1){rnorm(n,0,1)}
    pheta=rpheta/pi*180
    rotation<-matrix(c(cos(pheta),-sin(pheta),sin(pheta),cos(pheta)),ncol=2)
    function(n=1)t(apply(as.matrix(cbind(rv(n)*scale[1],rv(n)*scale[2]))%*%rotation,1,function(x)c(translation[1]+x[1],translation[2]+x[2]))) }
