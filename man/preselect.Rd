\name{preselect}                   
\alias{preselect}



\title{
Score and preselect features.
}



\description{
Iterates all features to score them via mMs, Student's t-test, or mRMR.
Optionally, a list of not informative features can be obtained (for discarding
them).
}



\usage{
preselect(elist=NULL, columns1=NULL, columns2=NULL, label1="A", label2="B",
 log=NULL, discard.threshold=0.5, fold.thresh=1.5, discard.features=TRUE,
 mMs.above=1500, mMs.between=400, mMs.matrix1=NULL,
 mMs.matrix2=NULL, method=NULL)
}



\arguments{
  \item{elist}{ \code{EListRaw} or \code{EList} object (mandatory). }
  \item{columns1}{ column name vector (string vector) of group 1 (mandatory). }
  \item{columns2}{ column name vector (string vector) of group 2 (mandatory). }
  \item{label1}{ class label of group 1.}
  \item{label2}{ class label of group 2.}
  \item{log}{ indicates whether the data is in log scale (mandatory; note: if
  TRUE log2 scale is expected). }
  \item{discard.threshold}{ positive numeric between 0 and 1 indicating the
  maximum mMs or, respectively, the maximum t-test p-value for features to be
  included for further analysis. Default is \code{"0.5"}. }
  \item{fold.thresh}{ numeric indicating the minimum fold change for
  features to be included for further analysis. Default is \code{"1.5"}. }
  \item{discard.features}{ boolean indicating whether merely feature scores
  (i.e., mMs or t-test p-values) (=\code{"FALSE"}) or feature scores and a
  discard list (=\code{"TRUE"}) should be returned. Default is \code{"TRUE"}. }
  \item{mMs.above}{ mMs above parameter (integer). Default is \code{"1500"}. }
  \item{mMs.between}{ mMs between parameter (integer). Default is
  \code{"400"}. }
  \item{mMs.matrix1}{ precomputed mMs reference matrix (see \code{mMsMatrix()})
  for group 1 (mandatory). }
  \item{mMs.matrix2}{ precomputed mMs reference matrix (see \code{mMsMatrix()})
  for group 2 (mandatory). }
  \item{method}{preselection method ( \code{"mMs"}, \code{"tTest"},
  \code{"mrmr"}). Default is \code{"mMs"}. }
}



\details{
This function takes an \code{EListRaw} or \code{EList} object and group-specific
column vectors. Furthermore, the class labels of group 1 and group 2 are needed.
If \code{discard.features} is \code{"TRUE"} (default), all features that are
considered as not differential will be collected and returned for discarding.

If \code{method = "mMs"}, additionally precomputed mMs reference matrices (see
\code{mMsMatrix()}) for group 1 and group 2 will be needed to compute mMs values
(see Love B.) as scoring method. All mMs parameters (\code{mMs.above} and
\code{mMs.between}) can be set. The defaults are \code{"1500"} for
\code{mMs.above} and \code{"400"} for \code{mMs.between}.  Features having an
mMs value larger than \code{discard.threshold} (here: numeric between 0.0 and
1.0) or do not satisfy the minimal absolute fold change \code{fold.thresh} are
considered as not differential. 

If \code{method = "tTest"}, Student's t-test will be used as scoring method.
Features having a p-value larger than \code{discard.threshold} (here: numeric
between 0.0 and 1.0) or do not satisfy the minimal absolute fold change
\code{fold.thresh} are considered as not differential.

If \code{method = "mrmr"}, mRMR scores for all features will be computed as
scoring method (using the function \code{mRMR.classic()} of the CRAN R package
\code{mRMRe}). Features that are not the \code{discard.threshold} (here: integer
indicating a number of features) best features regarding their mRMR score are
considered as not differential.    
}



\value{
If \code{discard.features} is \code{"FALSE"}: matrix containing metadata,
feature scores and intensity values for the whole data set.

If \code{discard.features} is \code{"TRUE"}, a list containing:
  \item{results}{matrix containing metadata, feature scores and intensity values
  for the whole data set.}
  \item{discard}{vector containing row indices (= features) for discarding
  features considered as not differential.}
}



\references{
Love B: The Analysis of Protein Arrays. In: Functional Protein Microarrays in
Drug Discovery. CRC Press; 2007: 381-402.

The software \code{"Prospector"} for ProtoArray analysis can be downloaded from
the Thermo Fisher Scientific web page (\url{https://www.thermofisher.com}).

The R package mRMRe can be downloaded from CRAN.
See also: De Jay N, Papillon-Cavanagh S, Olsen C, El-Hachem N, Bontempi G,
Haibe-Kains B. mRMRe: an R package for parallelized mRMR ensemble feature
selection. Bioinformatics 2013.

The package \code{limma} by Gordon Smyth et al. can be downloaded from
Bioconductor (\url{https://www.bioconductor.org}).

Smyth, G. K. (2005). Limma: linear models for microarray data. In:
Bioinformatics and Computational Biology Solutions using R and Bioconductor,
R. Gentleman, V. Carey, S. Dudoit, R. Irizarry, W. Huber (eds.), Springer,
New York, pages 397-420.
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
preselect(elist, columns1=c1, columns2=c2, label1="AD", label2="NDC", log=FALSE,
 discard.threshold=0.5, fold.thresh=1.5, discard.features=TRUE, method="tTest")
}



\keyword{ Feature selection }