\name{normalizeArrays}
\alias{normalizeArrays}



\title{
Normalize microarray data.
}



\description{
Normalizes \code{EListRaw} data and returns an \code{EList} object containing
normalized data in log scale.
}



\usage{
normalizeArrays(elist = NULL, method = "quantile", cyclicloess.method = "pairs", 
group1 = NULL, group2 = NULL, controls.elist=NULL, gpr.path=NULL,
targets.path=NULL, contr.names=NULL, output.path=NULL)
}



\arguments{
  \item{elist}{ \code{EListRaw} object containing raw data to be normalized
  (mandatory).}
  \item{method}{ string indicating the normalization method
  (\code{"cyclicloess"}, \code{"quantile"} or \code{"vsn"}) to be used
  (mandatory).}
  \item{cyclicloess.method}{ string indicating which type of cyclicloess
  normalization (\code{"pairs"}, \code{"fast"}, \code{"affy"}) should be
  performed (optional).}
  \item{group1}{ vector of integers (column indices) indicating all group 1
  samples (optional).}
  \item{group2}{ vector of integers (column indices) indicating all group 2
  samples (optional).}
  \item{controls.elist}{ \code{EListRaw} object containing intensity values of
  control features for RLM normalization for ProtoArrays (optional).}
  \item{gpr.path}{ when controls.elist is not defined: path to original gpr
  files for loading control feature intensities for ProtoArray RLM normalization
  (optional).}
  \item{targets.path}{ when controls.elist is not defined: path to original
  targets file for loading control feature intensities for ProtoArray RLM
  normalization (optional).}
  \item{contr.names}{ character vector containing the names of the control
  features for ProtoArray RLM normalization (optional).}
  \item{output.path}{ output.path for ProtoArray RLM normalization (optional).}
}



\details{
This function is partially a wrapper to \code{limma}'s function
\code{normalizeBetweenArrays()} for inter-array normalization featuring optional
groupwise normalization when the arguments \code{group1} AND \code{group2} are
assigned. For more information on \code{"cyclicloess"}, \code{"quantile"} or
\code{"vsn"} see the documentation of the \code{limma} package. Furthermore, for
ProtoArrays robust linear normalization (\code{"rlm"}, see Sboner A. et al.) is
provided. 
}



\value{
An EList object with the normalized data in log scale is returned.
}



\references{
The package \code{limma} by Gordon Smyth et al. can be downloaded from
Bioconductor (\url{http://www.bioconductor.org/}).

Smyth, G. K. (2005). Limma: linear models for microarray data. In:
Bioinformatics and Computational Biology Solutions using R and Bioconductor,
R. Gentleman, V. Carey, S. Dudoit, R. Irizarry, W. Huber (eds.), Springer,
New York, pages 397-420.

Sboner A. et al., Robust-linear-model normalization to reduce technical
variability in functional protein microarrays. J Proteome Res 2009,
8(12):5451-5464.
}



\author{
Michael Turewicz, \email{michael.turewicz@rub.de}
}



\note{
Groupwise normalization will be performed when the arguments \code{group1} AND
\code{group2} are assigned. 
}



\examples{
cwd <- system.file(package="PAA")
load(paste(cwd, "/extdata/Alzheimer.RData", sep=""))
elist <- elist[elist$genes$Block < 10,]
normalized.elist <- normalizeArrays(elist=elist, method="quantile")
}



\keyword{ Pre-processing }