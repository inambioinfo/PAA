\name{plotArray}
\alias{plotArray}



\title{
Plot ProtoArray expression intensities in the original arrangement mimicing the
original scan image.
}



\description{
Uses the \code{"Block"}, \code{"Row"} and \code{"Column"} information of an
\code{EList} or \code{EListRaw} object to resemble the original positions on the
array(s). The resulting plot is similar to the original scan image of the
considered array(s). Thus, this function is a visualization tool that can be
used to visualize ProtoArrays for which the original scan image is not
available. Visual inspection of the spatial expression pattern can then identify
possible local tendencies and strong spatial biases. Moreover, the array can be
inspected at all stages of the pre-processing workflow in order to check the
impact of the particular methods that have been applied.  
}



\usage{
plotArray(elist=NULL, idx=NULL, data.type="fg", log=NULL, normalized=NULL,
  protoarray.aggregation=NULL, colpal="heat.colors", output.path=NULL)
}



\arguments{
  \item{elist}{\code{EList} or \code{EListRaw} object (mandatory).}
  \item{idx}{integer, vector of integers or the string \code{"all"} indicating
    the column indices of the sample(s) for drawing the plot(s) (mandatory).}
  \item{data.type}{string indicating whether the foreground (\code{"fg"}) or
    background (\code{"bg"}) data should be plotted. The default is \code{"fg"} 
   (optional).}
  \item{log}{logical indicating whether the input data is logarithmized. If
    FALSE, the log2 will be plotted (mandatory).}
  \item{normalized}{logical indicating whether \code{elist} was normalized
    (mandatory).}
  \item{protoarray.aggregation}{string indicating whether the data stored in
    \code{elist} has been aggregated and, if this is the case, which method has
    been used by the function \code{loadGPR()}. Possible values are
    \code{"min"}, \code{"mean"} and \code{"none"}(mandatory).}
  \item{colpal}{string indicating the color palette for the plot(s). The default
    is \code{"heat.colors"} (optional).}
  \item{output.path}{string indicating the output path for the plots
    (optional).} 
}



\details{
This function allows plotting of ProtoArray data using the \code{gplots}
function \code{heatmap.2()} for visual quality control. The data obtained 
from an \code{EList} or \code{EListRaw} object is re-ordered and represented 
in the same way the spots are ordered on the actual ProtoArray. Consequently,
the resulting plot is similar to the original scan image of the considered
array. This allows for visual control and assessment of possible patterns in
spatial distribution.

Mandatory arguments are \code{elist}, \code{idx}, \code{log}, \code{normalized}
and \code{protoarray.aggregation}. While \code{elist} specifies the \code{EList}
or \code{EListRaw} object to be used, \code{idx} designates the ProtoArray
column index in \code{elist} to plot a single array from the \code{EList}
object. Alternatively, a vector (e.g., \code{1:5}) or the string \code{"all"}
can be designated to include multiple, respectively, all arrays that were
imported.

Furthermore, \code{data.type} allows for plotting of \code{"fg"}, foreground
data (i.e., \code{elist$E} and \code{elist$C}), which is the default or
\code{"bg"}, background data (i.e., \code{elist$Eb} and \code{elist$Cb}).

The normalization approaches of PAA which comprise also data logarithmization
do not include control data. With \code{normalized=TRUE} it is indicated that
the input data was normalized, so the control data will be logarithmized before
plotting as well. However, since the complete data (foreground and background
values of protein features and control spots) can be logarithmized regardless of
normalization the argument \code{log} states whether the designated data is
already logarithmized.

The parameter \code{protoarray.aggregation} indicates whether the ProtoArray
data has been aggregated by \code{loadGPR()} and, if so, which method has been
used.

Moreover, the parameter \code{colpal} defines the color palette that will
be used for the plot. Some exemplary values are \code{"heat.colors"} (default),
\code{"terrain.colors"}, \code{"topo.colors"}, \code{"greenred"} and
\code{"bluered"}.

Finally, the output path optionally can be specified with the argument
\code{output.path} to save the plot(s). Then, one or more tiff file(s)
containing the corresponding plot(s) are saved into the subfolder "array_plots".
}

\value{
No value is returned.  
}


\references{
The package \code{gplots} by Gregory R. Warnes et al. can be downloaded from
CRAN (\url{http://CRAN.R-project.org/package=gplots}).  

Gregory R. Warnes, Ben Bolker, Lodewijk Bonebakker, Robert Gentleman, 
Wolfgang Huber, Andy Liaw, Thomas Lumley, Martin Maechler, Arni Magnusson, 
Steffen Moeller, Marc Schwartz and Bill Venables (2015). gplots: Various R 
Programming Tools for Plotting Data. R package version 2.17.0. 
http://CRAN.R-project.org/package=gplots
}



\author{
Daniel Bemmerl and Michael Turewicz \email{michael.turewicz@rub.de}
}


\note{
This function is designed to represent ProtoArray data. Please note the
instructions of the PAA function \code{loadGPR()}. Note that the data has to be
imported as \code{array.type = "ProtoArray"} (default) to include controls and
that the data can be imported without aggregation by \code{loadGPR()} (when
\code{protoarray.aggregation="none"}) in order to inspect the array visually
before duplicate aggregation.  
}


\examples{
cwd <- system.file(package="PAA")
load(paste(cwd, "/extdata/BadData.RData", sep=""))
plotArray(elist=bad.elist, idx=1, data.type="bg", log=FALSE, normalized=FALSE,
  protoarray.aggregation="none")
}


\keyword{Pre-processing}