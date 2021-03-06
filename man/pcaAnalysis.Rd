% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/pcaAnalysis.r
\name{pcaAnalysis}
\alias{pcaAnalysis}
\title{Principle Componenet Analysis}
\usage{
pcaAnalysis(heightFile, catFile, normMethod = "colQn")
}
\description{
Principle Componenet Analysis
}
\examples{
peakFile <- "sorted_peaks.bed"
   catFile <- "catagories"
   pcs<-pcaAnalysis(peakFile,catFile)
   pc1=3
   pc2=6
   plotPCs(pcs,cbind(pc1,pc2),normData,cats,c( paste("PC:",as.character(pc1),sep="" ),paste("PC:",as.character(pc2),sep="" ),mock[i]))
}

