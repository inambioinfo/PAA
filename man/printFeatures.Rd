\name{printFeatures}
\alias{printFeatures}



\title{
Print biomarker candidate panel into a table. 
}



\description{
Creates a table containing the selected biomarker candidate panel.  
}



\usage{
printFeatures(features = NULL, elist = NULL, output.path = NULL)
}



\arguments{
  \item{features}{ vector containing the selected features as \code{"BRC"}-IDs
  (mandatory). }
  \item{elist}{ \code{EListRaw} or \code{EList} object containing all intensity
  data in log scale (mandatory). }
  \item{output.path}{ string indicating the folder where the table will be saved
  as a txt file (optional).}
}



\details{
Creates a table containing the selected biomarker candidate panel as well as
additional information. When \code{output.path} is defined this table will be
saved in a txt file (\code{"candidates.txt"}).
}



\value{
Table containing the selected biomarker panel.
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
 
printFeatures(features=selectFeatures.results$features, elist=elist)
}



\keyword{ Feature selection }