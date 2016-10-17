\name{diffAnalysis}                   
\alias{diffAnalysis}



\title{
Differential analysis.
}



\description{
Performs a univariate differential analysis.
}



\usage{
diffAnalysis(input=NULL, label1=NULL, label2=NULL, class1=NULL, class2=NULL,
output.path=NULL, mMs.matrix1=NULL, mMs.matrix2=NULL, above=1500,
between=400, features=NULL, feature.names=NULL)
}



\arguments{
  \item{input}{ \code{EList$E}- or \code{EListRaw$E}-matrix extended by row
  names comprising BRC-IDs of the corresponding features (mandatory; note: it is
  expected that this matrix is in original scale and not in log2 scale). }
  \item{label1}{ vector of column names for group 1 (mandatory). }
  \item{label2}{ vector of column names for group 2 (mandatory). }
  \item{class1}{ label of group 1 (mandatory). }
  \item{class2}{ label of group 2 (mandatory). }
  \item{output.path}{ string indicating a path for saving the results
  (optionally). }
  \item{mMs.matrix1}{ precomputed mMs reference matrix (see \code{mMsMatrix()})
  for group 1 (mandatory). }
  \item{mMs.matrix2}{ precomputed mMs reference matrix (see \code{mMsMatrix()})
  for group 2 (mandatory). }
  \item{above}{ mMs above parameter (integer). Default is \code{"1500"}. }
  \item{between}{ mMs between parameter (integer). Default is \code{"400"}. }
  \item{features}{ vector of row indices (optional). }
  \item{feature.names}{ vector of corresponding feature names (additionally to
  \code{features}). }
}



\details{
This function takes an \code{EList$E}- or \code{EListRaw$E}-matrix (e.g.,
\code{temp <- elist$E}) extended by row names comprising BRC-IDs of the
corresponding features. The BRC-IDs can be created via:\cr
\code{brc <- paste(elist$genes[,1], elist$genes[,3], elist.$genes[,2])}.\cr
The BRC-row names can be defined as follows: \code{rownames(temp) <- brc}.
Furthermore, the corresponding column name vectors, group labels and
mMs-parameters are needed to perform the univariate differential analysis. This
analysis covers inter alia p-value computation, p-value adjustment (method:
Benjamini & Hochberg, 1995), and fold change computation. Since the results
table is usually large, a path for saving the results can be defined via
\code{output.path}. Optionally, a vector of row indices (\code{features}) and
additionally (not mandatory for subset analysis) a vector of corresponding
feature names (\code{feature.names}) can be forwarded to perform the analysis
for a feature subset.     
}



\value{
A matrix containing the analysis results is returned.
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
mMs.matrix1 <- mMs.matrix2 <- mMsMatrix(x=20, y=20)
temp <- elist$E
rownames(temp) <- paste(elist$genes[,1], elist$genes[,3], elist$genes[,2])
diffAnalysis(input=temp, label1=c1, label2=c2, class1="AD", class2="NDC",
 mMs.matrix1=mMs.matrix1, mMs.matrix2=mMs.matrix2, above=1500,
 between=400)
}



\keyword{ Differential analysis }
