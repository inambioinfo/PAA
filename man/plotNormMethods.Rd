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
plotNormMethods(elist = NULL, include.rlm=FALSE, controls.elist=NULL,
gpr.path=NULL, targets.path=NULL, contr.names=NULL, output.path = NULL)
}



\arguments{
  \item{elist}{ \code{EListRaw} object containing raw data (mandatory). }
  \item{include.rlm}{ boolean indicating whether RLM normalization should be
  included (for ProtoArrays only, deafault: FALSE). }
  \item{controls.elist}{ \code{EListRaw} object containing intensity values of
  control features for RLM normalization for ProtoArrays (optional). }
  \item{gpr.path}{ when controls.elist is not defined: path to original gpr
  files for loading control feature intensities for ProtoArray RLM normalization
  (optional). }
  \item{targets.path}{ when controls.elist is not defined: path to original
  targets file for loading control feature intensities for ProtoArray RLM
  normalization (optional). }
  \item{contr.names}{ character vector containing the names of the control
  features for ProtoArray RLM normalization (optional). }
  \item{output.path}{ string indicating a folder for saving the boxplots as tiff
  files (optional). }
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
elist <- elist[elist$genes$Block < 10,]
plotNormMethods(elist=elist)
}



\keyword{ Pre-processing }