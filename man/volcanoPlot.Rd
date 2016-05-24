\name{volcanoPlot}                   
\alias{volcanoPlot}



\title{
Draw a volcano plot.
}



\description{
Draws a volcano plot to visualize differential features.
}



\usage{
volcanoPlot(elist=NULL, group1=NULL, group2=NULL, log=NULL, method="tTest",
p.thresh=NULL, fold.thresh=NULL, output.path=NULL, tag="", mMs.matrix1=NULL,
mMs.matrix2=NULL, above=1500, between=400)
}



\arguments{
  \item{elist}{ \code{EList} or \code{EListRaw} object (mandatory). }
  \item{group1}{ vector of column names for group 1 (mandatory). }
  \item{group2}{ vector of column names for group 2 (mandatory). }
  \item{log}{ indicates whether the data is in log scale (mandatory; note: if
  TRUE log2 scale is expected; mandatory). }
  \item{method}{ method for p-value computation: \code{"tTest"} or \code{"mMs"}.
  Default is \code{"tTest"}. }  
  \item{p.thresh}{ positive float number between 0 and 1 indicating the maximum
  p-value for features to be considered as differential (e.g., \code{"0.5"}).
  This argument is optional. }
  \item{fold.thresh}{ float number indicating the minimum fold change for
  features to be considered as differential (e.g., \code{"1.5"}). This argument
  is optional. }
  \item{output.path}{ string indicating a path for saving the plot (optional). }
  \item{tag}{ string that can be used for tagging the saved plot (optional). }
  \item{mMs.matrix1}{ a precomputed M score reference matrix (see
  \code{mMsMatrix()}) for group 1 (mandatory when method = "mMs"). }
  \item{mMs.matrix2}{ a precomputed M score reference matrix (see
  \code{mMsMatrix()}) for group 2 (mandatory when method = "mMs"). }
  \item{above}{ M score above parameter (integer). Default is \code{"1500"}. }
  \item{between}{ M score between parameter (integer). Default is
  \code{"400"}. }
}



\details{
This function takes an \code{EList} or \code{EListRaw} object and the
corresponding column name vectors to draw a volcano plot. To visualize
differential features, thresholds for p-values and fold changes can be defined.
Furthermore, the p-value computation method (\code{"mMs"} or \code{"tTest"}) can
be set. When an output path is defined (via \code{output.path}) the plot will be
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
volcanoPlot(elist=elist, group1=c1, group2=c2, log=FALSE, method="tTest",
 p.thresh=0.01, fold.thresh=2)
}



\keyword{ Differential analysis }