\name{batchAdjust}
\alias{batchAdjust}



\title{
Adjust microarray data for batch effects.
}



\description{
Adjusts \code{EListRaw} or \code{EList} data for batch/lot effects.
}



\usage{
batchAdjust(elist=NULL, log=NULL)
}



\arguments{
  \item{elist}{\code{EList} or \code{EListRaw} object containing the data to be
  adjusted (mandatory).}
  \item{log}{logical indicating whether the data is in log scale (mandatory;
  note: if TRUE log2 scale is expected).}
}



\details{
This is a wrapper to \code{sva}'s function \code{ComBat()} for batch adjustment
using the empirical Bayes approach. To use batchAdjust the targets information
of the \code{EList} or \code{EListRaw} object must contain the columns
\code{"Batch"} (containing batch/lot information for each particular array) and
\code{"Group"} (containing experimental group information for each particular
array). 
}



\value{
An \code{EListRaw} or \code{EList} object with the adjusted data in log scale is
returned.
}



\references{
The package \code{sva} by Jeffrey T. Leek et al. can be downloaded from
Bioconductor (\url{http://www.bioconductor.org/}).

Johnson WE, Li C, and Rabinovic A (2007) Adjusting batch effects in microarray
expression data using empirical Bayes methods. Biostatistics 8:118-27.
}



\author{
Michael Turewicz, \email{michael.turewicz@rub.de}
}



\note{
The targets information of the \code{EListRaw} or \code{EList} object must
contain the columns \code{"Batch"} and \code{"Group"}. 
}



\examples{
cwd <- system.file(package="PAA")
load(paste(cwd, "/extdata/Alzheimer.RData", sep=""))
elist <- elist[elist$genes$Block < 10,]
elist <- batchAdjust(elist=elist, log=FALSE)
}



\keyword{ Preprocessing }