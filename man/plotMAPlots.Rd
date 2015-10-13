\name{plotMAPlots}
\alias{plotMAPlots}



\title{
Check normalization results with MA plots.
}



\description{
Draws MA plots of raw data and data after all kinds of normalization provided by
PAA.
}



\usage{
plotMAPlots(elist = NULL, idx="all", include.rlm=FALSE, controls="internal",
output.path = NULL)
}



\arguments{
  \item{elist}{ \code{EListRaw} object containing raw data (mandatory).}
  \item{idx}{ integer indicating the column index of the sample for drawing MA
  plots or the string 'all' for drawing MA plots for all samples
  (default: all).}
  \item{include.rlm}{ boolean indicating whether RLM normalization should be
  included (for ProtoArrays only, deafault: FALSE).}
  \item{controls}{ sring indicating the ProtoArray controls for \code{rlm}
  normalization (optional). Valid options are \code{"internal"} (default),
  \code{"external"}, \code{"both"} or a regular expression defining a specific
  control or a specific set of controls.}
  \item{output.path}{ string indicating the folder where the tiff files will be
  saved (mandatory when idx='all').}
}



\details{
When \code{idx="all"} (default) for each microarray a tiff file containing MA
plots for raw data, cyclicoess normalized data, quantile normalized data and vsn
normalized data (and, optionally, for ProtoArrays, rlm normalized data) will be
created. When \code{idx} is an integer indicating the column index of a
particular sample, MA plots only for this sample will be created. For A and
M value computation the artificial median array is used as reference signal.
All figures can be saved in \code{output.path} (mandatory when
\code{idx="all"}). The resulting MA plots can be used to compare the results of
the different normalization methods.
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
plotMAPlots(elist=elist, idx=1)
}



\keyword{ Pre-processing }