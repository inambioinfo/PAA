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
controls="internal", group1 = NULL, group2 = NULL, output.path=NULL)
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
  \item{controls}{ sring indicating the ProtoArray controls for \code{rlm}
  normalization (optional). Valid options are \code{"internal"} (default),
  \code{"external"}, \code{"both"} or a regular expression defining a specific
  control or a specific set of controls.}
  \item{group1}{ vector of integers (column indices) indicating all group 1
  samples (optional).}
  \item{group2}{ vector of integers (column indices) indicating all group 2
  samples (optional).}
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

For \code{rlm} normalization (\code{method = "rlm"}) the additional argument
\code{controls} needs to be specified in order to select a set of controls used
for normalization. Valid options are \code{"internal"} (default),
\code{"external"} and \code{"both"} which refer to the following sets of
ProtoArray controls:
\itemize{
    \item{internal}{ The set of all internal controls spotted on the ProtoArray.
    The human-IgG series and anti-human-IgG series, which respond to serum and
    secondary antibodies. }
    \item{external}{ The V5-CMK1 series spotted on the ProtoArray which responds
    to exogenously added anti-V5 antibody (external control). }
    \item{both}{ The combined set of both the internal and the external
    controls (i.e., the human-IgG and anti-human-IgG series and the V5-CMK1
    series). }
}

Moreover, via \code{controls} a regular expression can be passed in order to
select a more specific group of controls. Please check the column \code{"Name"}
in your gpr files in order to obtain the complete list of names of all controls
spotted on the ProtoArray. In the following some examples of valid regular
expressions are given:
\itemize{
    \item{\code{"^HumanIg"}} {Only human IgGs and IgAs are selected (esp.,
    no anti-human Igs).}
    \item{\code{"Anti-HumanIgA"} {Only anti-human-IgAs are selected (esp.,
    no human IgGs and IgAs).}}
    \item{\code{"(Anti-HumanIg|^V5control|BSA|ERa)"} {Only anti-human IgGs and
    anti-human IgAs, the V5-CMK1 series, BSA and ERa are selected.}}
    \item{\code{"HumanIgG"} {Only human IgGs and anti-human IgGs are selected}.}
    \item{\code{"V5control"} {Only the V5-CMK1 series is selected.}}
}    
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