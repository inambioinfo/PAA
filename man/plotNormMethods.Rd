\name{plotNormMethods}
\alias{plotNormMethods}



\title{
Check normalization results with boxplots.
}



\description{
Draws boxplots (one boxplot per sample) of raw data and data after all kinds of
normalization provided by PAA.  
}



\usage{
plotNormMethods(elist = NULL, include.rlm=FALSE, controls="internal",
output.path = NULL)
}



\arguments{
  \item{elist}{ \code{EListRaw} object containing raw data (mandatory).}
  \item{include.rlm}{ boolean indicating whether RLM normalization should be
  included (for ProtoArrays only, deafault: FALSE).}
  \item{controls}{ sring indicating the ProtoArray controls for \code{rlm}
  normalization (optional). Valid options are \code{"internal"} (default),
  \code{"external"}, \code{"both"} or a regular expression defining a specific
  control or a specific set of controls.}
  \item{output.path}{ string indicating a folder for saving the boxplots as tiff
  files (optional).}
}



\details{
For each normalization approach sample-wise boxplots are created. All boxplots
can be saved as high-quality tiff files (when an output path has been specified
via the argument \code{output.path}). The resulting boxplots can be used to
compare the results of different normalization methods.
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
elist <- elist[elist$genes$Block == 1,]
plotNormMethods(elist=elist)
}



\keyword{ Pre-processing }