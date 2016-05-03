\name{loadGPR}
\alias{loadGPR}



\title{
A function to read raw data from gpr files.
}



\description{
Constructs an \code{EListRaw} object from a set of gpr files containing
ProtoArray data or other one-color microarray data.
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
  Alternatively, no aggregation will be performed, if \code{"none"} is chosen.
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
\code{"ProtoArray"} (default), duplicate spots can be aggregated. The
corresponding method (\code{"min"}, \code{"mean"} or \code{"none"}) can be
specified via \code{protoarray.aggregation}. As another ProtoArray-specific
feature, control spot data and information will be stored in additional
components of the returned object (see below). Arguments \code{array.columns}
and \code{array.annotation} define the columns where \code{read.maimages()} will
find foreground and background intensity values as well as other important
columns. For \code{array.annotation} the default columns \code{"Block"},
\code{"Column"}, \code{"Row"}, \code{"Description"}, \code{"Name"} and
\code{"ID"} are mandatory (exception: the \code{"Description"} column for
ProtoArrays since a very general column can be constructed from the column
\code{"Name"}).         
}



\value{
An extended object of class \code{EListRaw} (see the documentation of
\code{limma} for details) is returned. If \code{array.type} is set to
\code{"ProtoArray"} (default), the object provides additional components for
control spot data: \code{C}, \code{Cb} and \code{cgenes} which are analogous to
the probe spot data \code{E}, \code{Eb} and \code{genes}.
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
default settings of \code{loadGPR()} and should be renamed to the default column
names (see also the exemplary gpr files accompanying PAA as a reference for the
default column names). At worst, important columns in your gpr files may be
completely missing and should be added in order to provide all information
needed by PAA.

Note that if \code{array.type} is not \code{"ProtoArray"}, neither aggregation
will be done nor controls components will be added to the returned object of
class \code{EListRaw}.  
}



\examples{
gpr <- system.file("extdata", package="PAA") 
targets <- list.files(system.file("extdata", package="PAA"),
 pattern = "dummy_targets", full.names=TRUE)   
elist <- loadGPR(gpr.path=gpr, targets.path=targets)
}



\keyword{Input/output}