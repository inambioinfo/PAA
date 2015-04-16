\name{batchFilter}                   
\alias{batchFilter}



\title{
Remove differential features regarding array batches/lots.
}



\description{
Finds differential features regarding array batches/lots and removes them.
}



\usage{
batchFilter(elist=NULL, lot1=NULL, lot2=NULL, p.thresh=0.05, fold.thresh=1.5,
output.path=NULL)
}



\arguments{
  \item{elist}{\code{EList} or \code{EListRaw} object (mandatory).}
  \item{lot1}{vector of column names for group 1 (mandatory).}
  \item{lot2}{vector of column names for group 2 (mandatory).}
  \item{p.thresh}{positive float number between 0 and 1 indicating the maximum
  Student's t-test p-value for features to be considered as differential (e.g.,
  \code{"0.5"}).}
  \item{fold.thresh}{float number indicating the minimum fold change for
  features to be considered as differential (e.g., \code{"1.5"}).}
  \item{output.path}{string indicating a path for saving results (optional).}
}



\details{
This function takes an \code{EList} or \code{EListRaw} object (see \code{limma}
documentation) and the batch-specific column name vectors \code{lot1} and
\code{lot2} to find differential features regarding batches/lots. For this
purpose, thresholds for p-values (Student's t-test) and fold changes can be
defined. To visualize the differential features a volcano plot is drawn. Then,
differential features are removed and the remaining data are returned. When an
output path is defined (via \code{output.path}) volcano plots and result files
are saved on the hard disk.     
}



\value{
An \code{EList} or \code{EListRaw} object without differential features
regarding array batches/lots.
}



\author{
Michael Turewicz, \email{michael.turewicz@rub.de}
}



\examples{
cwd <- system.file(package="PAA")
load(paste(cwd, "/extdata/Alzheimer.RData", sep=""))
elist <- elist[elist$genes$Block < 10,]
lot1 <- elist$targets[elist$targets$Batch=='Batch1','ArrayID']
lot2 <- elist$targets[elist$targets$Batch=='Batch2','ArrayID']
elist <- batchFilter(elist=elist, lot1=lot1, lot2=lot2, p.thresh=0.001,
  fold.thresh=3)
}



\keyword{ Pre-processing }