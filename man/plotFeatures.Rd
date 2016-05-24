\name{plotFeatures}
\alias{plotFeatures}



\title{
Plot intensities of features.
}



\description{
Plots intensities of all given features (one sub-plot per feature) in group-
specific colors.  
}



\usage{
plotFeatures(features = NULL, elist = NULL, n1 = NULL, n2 = NULL, 
 group1 = "group1", group2 = "group2", output.path = NULL)
}



\arguments{
  \item{features}{vector containing \code{"BRC"}-IDs (mandatory).}
  \item{elist}{\code{EListRaw} or \code{EList} object containing all intensity
  data in log2 scale (mandatory).}
  \item{n1}{integer indicating the sample size of group 1 (mandatory).}
  \item{n2}{integer indicating the sample size of group 2 (mandatory).}
  \item{group1}{class label of group 1.}
  \item{group2}{class label of group 2.}
  \item{output.path}{string indicating the folder where the figure will be saved
  (optional).}
}



\details{
Plots intensities of given features (e.g., selected by the function
\code{selectFeatures()}) in group-specific colors (one sub-plot per feature).
All sub-plots are aggregated to one figure. When the argument
\code{output.path} is not NULL this figure will be saved in a tiff file in
\code{output.path}. This function can be used to check whether the selected
features are differential.
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
#elist <- elist[elist$genes$Block < 10,]

#c1 <- paste(rep("AD",20), 1:20, sep="")
#c2 <- paste(rep("NDC",20), 1:20, sep="")

#pre.sel.results <- preselect(elist=elist, columns1=c1, columns2=c2, label1="AD",
# label2="NDC", discard.threshold=0.1, fold.thresh=1.9, discard.features=TRUE,
# method="tTest")
#elist <- elist[-pre.sel.results$discard,]

#selectFeatures.results <- selectFeatures(elist,n1=20,n2=20,label1="AD",
# label2="NDC",selection.method="rf.rfe",preselection.method="none",subruns=2,
# k=2,candidate.number=20,method="frequency")

load(paste(cwd, "/extdata/selectFeaturesResultsFreq.RData", sep="")) 
plotFeatures(features=selectFeatures.results$features, elist=elist, n1=20,
 n2=20, group1="AD", group2="NDC")
}



\keyword{ Feature selection }