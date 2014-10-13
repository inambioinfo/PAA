\name{plotFeaturesHeatmap}
\alias{plotFeaturesHeatmap}



\title{
Plot feature intensities as a heatmap.
}



\description{
Plots intensities of given features as a heatmap.  
}



\usage{
plotFeaturesHeatmap(features = NULL, elist = NULL, n1 = NULL, n2 = NULL, 
 output.path = NULL, description=FALSE)
}



\arguments{
  \item{features}{vector containing the selected features as \code{"BRC"}-IDs
  (mandatory).}
  \item{elist}{\code{EListRaw} or \code{EList} object containing all intensity
  data in log scale (mandatory).}
  \item{n1}{integer indicating the sample size of group 1 (mandatory).}
  \item{n2}{integer indicating the sample size of group 2 (mandatory).}
  \item{output.path}{path for saving the heatmap as a tiff file (default:
  NULL).}
  \item{description}{if TRUE, features will be described via protein names
  instead of uniprot accessions (default: FALSE).}
}



\details{
Plots intensities of all features given in the vector \code{features} via their
corresponding \code{"BRC"}-IDs as a heatmap. If \code{description} is TRUE
(default: FALSE), features will be described via protein names instead of
uniprot accessions. Furthermore, if \code{output.path} is not NULL, the heatmap
will be saved as a tiff file in \code{output.path}. This function can be used to
check whether the selected features are differential.
}



\value{
No value is returned.
}



\author{
Michael Turewicz, \email{michael.turewicz@rub.de}
}



\examples{
cwd <- system.file(package="PAA")
load(paste(cwd, "/extdata/Alzheimer.RData", sep=""))
elist <- elist[elist$genes$Block < 10,]

c1 <- paste(rep("AD",20), 1:20, sep="")
c2 <- paste(rep("NDC",20), 1:20, sep="")

pre.sel.results <- preselect(elist=elist, columns1=c1, columns2=c2, label1="AD",
 label2="NDC", discard.threshold=0.1, fold.thresh=2, discard.features=TRUE,
 method="tTest")
elist <- elist[-pre.sel.results$discard,]

selectFeatures.results <- selectFeatures(elist,n1=20,n2=20,label1="AD",
 label2="NDC",selection.method="rf.rfe",preselection.method="none",subruns=2,
 k=2,candidate.number=20,method="frequency")
 
plotFeaturesHeatmap(features=selectFeatures.results$features, elist=elist,
 n1=20, n2=20, description=TRUE)
}



\keyword{ Feature selection }