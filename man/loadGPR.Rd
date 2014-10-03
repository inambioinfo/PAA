\name{loadGPR}
\alias{loadGPR}



\title{
A function to read raw data from gpr files.
}



\description{
Reads an EListRaw from a set of gpr files containing ProtoArray data or other
one-color data.
}



\usage{
loadGPR(gpr.path = NULL, targets.path = NULL, array.type = "ProtoArray", 
 protoarray.aggregation = "min", array.columns = list(E = "F635 Median", 
 Eb = "B635 Median"), 
 array.annotation = c("Block", "Column", "Row", "Description", "Name", "ID"))
}



\arguments{
  \item{gpr.path}{string indicating the path to a folder containing gpr files
  (mandatory).}
  \item{targets.path}{string indicating the path to targets file (see limma,
  mandatory).}
  \item{array.type}{string indicating whether the gpr files contain ProtoArray
  data (default is \cr\code{"ProtoArray"}). In any other case no duplicate
  aggregation will be performed (mandatory).}
  \item{protoarray.aggregation}{string indicating which type of ProtoArray spot
  duplicate aggregation should be performed. If \code{"min"} is chosen, the
  value for the corresponding feature will be the minimum of both duplicate
  values. If \code{"mean"} is chosen, the arithmetic mean will be computed.
  The default is \code{"min"} (optional).}
  \item{array.columns}{list containing the column names for foreground
  intensities (E) and background intensities (Eb) in the gpr files that is
  passed to \code{limma}'s \code{"read.maimages"} function (optional).}
  \item{array.annotation}{string vector containing further mandatory column
  names that are passed to limma (optional).}
}



\details{
This function is partially a wrapper to \code{limma}'s function
\code{read.maimages()} featuring optional duplicate aggregation for ProtoArray
data. Paths to a target file and to a folder containing gpr files (all gpr files
in that folder that are listed in the targets file will be read) are mandatory.
The folder \code{"R_HOME/library/PAA/extdata"} contains an exemplary target file
that can be used as a template. If \code{array.type} is set to
\code{"ProtoArray"} (default), duplicate spots will be aggregated. 
The corresponding method (\code{"min"} or \code{"mean"}) can be specified via
\code{protoarray.aggregation}. Arguments \code{array.columns} and
\code{array.annotation} define the columns where \code{read.maimages()} will
find foreground and background intensity values as well as other important
columns.          
}



\value{
An object of class \code{EListRaw} (see the documentation of \code{limma} for
details) is returned.
}



\references{
The package \code{limma} by Gordon Smyth et al. can be downloaded from
Bioconductor (\url{http://www.bioconductor.org/}).

Smyth, G. K. (2005). Limma: linear models for microarray data. In:
Bioinformatics and Computational Biology Solutions using R and Bioconductor,
R. Gentleman, V. Carey, S. Dudoit, R. Irizarry, W. Huber (eds.), Springer,
New York, pages 397-420.
}



\author{
Michael Turewicz, \email{michael.turewicz@rub.de}
}



\note{
Don't forget to check column names in your gpr files. They may differ from the
default settings of \code{loadGPR()}. Note that if \code{array.type} is not
\code{"ProtoArray"}, no aggregation will be done and all computations will be
performed by \code{limma}'s \code{read.maimages()}.  
}



\examples{
gpr <- system.file("extdata", package="PAA") 
targets <- list.files(system.file("extdata", package="PAA"),
 pattern = "dummy_targets", full.names=TRUE)   
elist <- loadGPR(gpr.path=gpr, targets.path=targets)
}



\keyword{Input/output}