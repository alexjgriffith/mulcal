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

