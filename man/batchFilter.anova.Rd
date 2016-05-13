\name{batchFilter.anova}                   
\alias{batchFilter.anova}



\title{
Remove features which are differential regarding microarray batches / lots in a
multi-batch scenario.
}



\description{
Finds features which are differential regarding at least two  microarray
batches / lots in a multi-batch scenario (i.e., > 2 batches) via one-way
analysis of variance (ANOVA) and removes them.
}



\usage{
batchFilter.anova(elist=NULL, log=NULL, p.thresh=0.05, fold.thresh=1.5,
output.path=NULL)
}



\arguments{
  \item{elist}{\code{EList} or \code{EListRaw} object (mandatory).}
  \item{log}{ logical indicating whether the data is in log scale (mandatory;
  note: if TRUE log2 scale is expected). }
  \item{p.thresh}{positive float number between 0 and 1 indicating the maximum
  Student's t-test p-value for features to be considered as differential (e.g.,
  \code{"0.5"}).}
  \item{fold.thresh}{float number indicating the minimum fold change for
  features to be considered as differential (e.g., \code{"1.5"}).}
  \item{output.path}{string indicating a path for saving results (optional).}
}



\details{
This function takes an \code{EList} or \code{EListRaw} object (see \code{limma}
documentation) to find features which are differential regarding at least two
microarray batches / lots in a multi-batch scenario (i.e., more than two
batches). For this purpose, thresholds for p-values obtained from an one-way
analysis of variance (ANOVA) and fold changes can be defined. To visualize the
differential features a volcano plot is drawn. Then, differential features are
removed and the remaining data are returned. When an output path is defined
(via \code{output.path}) volcano plots and result files are saved on the hard
disk.     
}



\value{
An \code{EList} or \code{EListRaw} object without differential features
regarding at least two microarray batches / lots.
}



\author{
Ivan Grishagin (Rancho BioSciences LLC, San Diego, CA, USA),
John Obenauer (Rancho BioSciences LLC, San Diego, CA, USA)
and Michael Turewicz (Ruhr-University Bochum, Bochum, Germany),
\email{michael.turewicz@rub.de}
}



\examples{
cwd <- system.file(package="PAA")
load(paste(cwd, "/extdata/Alzheimer.RData", sep=""))
elist <- elist[elist$genes$Block < 10,]
elist <- batchFilter.anova(elist=elist, log=FALSE, p.thresh=0.001,
 fold.thresh=3)
}



\keyword{ Preprocessing }