\name{pvaluePlot}                   
\alias{pvaluePlot}



\title{
Draw a p-value plot.
}



\description{
Draws a p-value plot to visualize the p-values for all features stored in a
\code{EList} or \code{EListRaw} object.
}



\usage{
pvaluePlot(elist=NULL, group1=NULL, group2=NULL, method="tTest",
output.path=NULL, tag="", mMs.matrix1=NULL, mMs.matrix2=NULL, above=1500,
between=400, adjust=FALSE)
}



\arguments{
  \item{elist}{ \code{EList} or \code{EListRaw} object (mandatory). }
  \item{group1}{ vector of column names for group 1 (mandatory). }
  \item{group2}{ vector of column names for group 2 (mandatory). }
  \item{method}{ method for p-value computation: \code{"tTest"} or \code{"mMs"}.
  Default is \code{"tTest"}. }
  \item{output.path}{ string indicating a path for saving the plot (optional). }
  \item{tag}{ string that can be used for tagging the saved plot (optional). }
  \item{mMs.matrix1}{ pre-computed M score reference matrix (see
  \code{mMsMatrix()}) for group 1 (mandatory when method = "mMs"). }
  \item{mMs.matrix2}{ pre-computed M score reference matrix (see
  \code{mMsMatrix()}) for group 2 (mandatory when method = "mMs"). }
  \item{above}{ M score above parameter (integer). Default is \code{"1500"}. }
  \item{between}{ M score between parameter (integer). Default is
  \code{"400"}. }
  \item{adjust}{ Boolean indicating whether p-values should be adjusted. Default
  is \code{FALSE}. }
}



\details{
This function takes an \code{EList} or \code{EListRaw} object and the
corresponding column name vectors to draw a plot of p-values for all features
stored in \code{elist} (sorted in increasing order and in log2 scale). The
p-value computation method (\code{"tTest"} or \code{"mMs"}) can be set via the
argument \code{method}. Furthermore, when \code{adjust=TRUE} adjusted p-values 
(method: Benjamini & Hochberg, 1995, computed via \code{p.adjust()}) will be
used. When an output path is defined (via \code{output.path}) the plot will be
saved as a tiff file.     
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
pvaluePlot(elist=elist, group1=c1, group2=c2, method="tTest", tag="_tTest",
 adjust=FALSE)
}



\keyword{ Differential analysis }
