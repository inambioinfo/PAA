\name{shuffleData}
\alias{shuffleData}



\title{
Shuffles class labels to obtain random groups.
}



\description{
Shuffles class labels of an \code{EList} or \code{EListRaw} object randomly to
obtain two random groups (e.g. "A" and "B").
}



\usage{
shuffleData(elist=NULL, n1=NULL, n2=NULL, label1="A", label2="B")
}



\arguments{
  \item{elist}{ \code{EList} or \code{EListRaw} object (mandatory). }
  \item{n1}{ sample size of random group 1 (mandatory). }
  \item{n2}{ sample size of random group 2 (mandatory). }
  \item{label1}{ class label of random group 1 (default: \code{"A"}). }
  \item{label2}{ class label of random group 2 (default: \code{"B"}). }
}



\details{
Shuffles class labels of an \code{EList} or \code{EListRaw} object randomly to
obtain two random groups (e.g. "A" and "B").     
}



\value{\code{EList} or \code{EListRaw} object with random groups.}



\author{
Michael Turewicz, \email{michael.turewicz@rub.de}
}



\examples{
cwd <- system.file(package="PAA")
load(paste(cwd, "/extdata/Alzheimer.RData", sep=""))
shuffleData(elist=elist, n1=20, n2=20, label1="A", label2="B")
}



\keyword{ Pre-processing }