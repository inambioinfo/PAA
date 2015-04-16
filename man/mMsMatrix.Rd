\name{mMsMatrix}
\alias{mMsMatrix}



\title{
Compute a reference minimum M statistic (n1 x n2)-matrix.
}



\description{
Computes a reference minimum M statistic (n1 x n2)-matrix (mMs matrix).
}



\usage{
mMsMatrix(x, y)
}



\arguments{
  \item{x}{integer, first dimension (i.e., number of samples in group 1) of the
  mMs matrix to be computed (mandatory).}
  \item{y}{integer, second dimension (i.e., number of samples in group 2) of the
  mMs matrix to be computed (mandatory).}
}



\details{
For feature pre-selection the "minimum M Statistic" (mMs) proposed by Love B. is
used. The mMs is a univariate measure that is sensitive to population subgroups. 
To avoid redundant mMs computations for a large number of features (e.g.,
ca. 9500 features on ProtoArray v5) a reference matrix containing all relevant
mMs values can be pre-computed. Therefore, only two parameters are needed: the
number of samples in group 1 (\code{n1}) and the number of samples in group 2
(\code{n2}). According to mMs definition for each matrix element (i,m) a mMs
value (= the probability of) for having m values in group 1 larger than the i-th
largest value in group 2 is computed.  
}



\value{
A  (n1 x n2)-matrix containing all mMs values for group 1 and group 2. 
}



\references{
Love B: The Analysis of Protein Arrays. In: Functional Protein Microarrays in
Drug Discovery. CRC Press; 2007: 381-402.
}



\author{
Michael Turewicz, \email{michael.turewicz@rub.de}
}



\note{
To check whether a feature is more prevalent in group 1 or in group 2, PAA needs
both the mMs for having m  values in group 1 larger than the i-th largest
element in group 2 as well as the mMs for having m values in group 2 larger than
the i-th largest element in group 1. Hence, always both needs be computed
\code{mMsMatrix(n1,n2)} and  \code{mMsMatrix(n2,n1)}.  
}



\examples{
#exemplary computation for a group 1 comprising 10 arrays and a group 2
#comprising 12 arrays
mMs.matrix1 <- mMsMatrix(x=10, y=12)
mMs.matrix2 <- mMsMatrix(x=12, y=10)
}



\keyword{ mMs }