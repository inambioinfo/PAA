\name{selectFeatures}
\alias{selectFeatures}



\title{
Select features using frequency-based or ensemble feature selection.
}



\description{
Performs a multivariate feature selection using frequency-based feature
selection (based on RF-RFE, RJ-RFE or SVM-RFE) or ensemble feature selection
(based on SVM-RFE). 
}



\usage{
selectFeatures(elist = NULL, n1 = NULL, n2 = NULL, label1 = "A", label2 = "B",
 cutoff = 10, selection.method = "rf.rfe", preselection.method = "mMs",
 subruns = 100, k = 10, subsamples = 10, bootstraps = 10,
 candidate.number = 300, above=1500, between=400,
 panel.selection.criterion="accuracy", importance.measure="MDA", ntree = 500,
 mtry = NULL, plot = FALSE, output.path = NULL, verbose = FALSE,
 method = "frequency")
}



\arguments{
  \item{elist}{ \code{EListRaw} or \code{EList} object containing all microarray
  data (mandatory). }
  \item{n1}{ integer indicating the sample number in group 1 (mandatory). }
  \item{n2}{ integer indicating the sample number in group 2 (mandatory). }
  \item{label1}{ class label of group 1 (default: "A"). }
  \item{label2}{ class label of group 2 (default: "B"). }
  \item{cutoff}{ integer indicating how many features will be selected
  (default: 10). }
  \item{selection.method}{ string indicating the feature selection method:
  \code{"rf.rfe"} (default), \code{"svm.rfe"} or \code{"rj.rfe"}. }
  \item{preselection.method}{ string indicating the feature pre-selection
  method: \code{"mMs"} (default), \code{"tTest"}, \code{"mrmr"} or
  \code{"none"}. }
  \item{subruns}{ integer indicating the number of re-sampling repeats to be
  performed (default: 100). }
  \item{k}{ integer indicating the number of k-fold cross validation subsets
  (default: 10, i.e., 10-fold CV). }
  \item{subsamples}{ integer indicating the number of subsamples for ensemble
  feature selection (default: 10). }
  \item{bootstraps}{ integer indicating the number of bootstrap samples for
  ensemble feature selection (default: 10). }
  \item{candidate.number}{ integer indicating how many features shall be pre-
  selected. Default is \code{"300"}. }
  \item{above}{ mMs above parameter (integer). Default is \code{"1500"}. }
  \item{between}{ mMs between parameter (integer). Default is \code{"400"}. }
  \item{panel.selection.criterion}{ string indicating the panel selection
  criterion: \code{"accuracy"} (default), \code{"sensitivity"}\cr or
  \code{"specificity"}. }
  \item{importance.measure}{ string indicating the random forest importance
  measure: \code{"MDA"} (default) or \code{"MDG"}. }
  \item{ntree}{ random forest parameter ntree (default: \code{"500"}). }
  \item{mtry}{ random forest parameter mtry (default: \code{sqrt(p)} where p
  is the number of predictors). }
  \item{plot}{ indicates whether performance plots shall be plotted
  (default: FALSE). }
  \item{output.path}{ string indicating the results output folder (optional). }
  \item{verbose}{ logical, indicating whether additional information shall be
  printed to the console (default: FALSE). }
  \item{method}{ the feature selection method: "frequency" (default) for
  frequency-based and "ensemble" for ensemble feature selection. }
}



\details{
This function takes an \code{EListRaw} or \code{EList} object, group-specific
sample numbers, group labels and parameters choosing and setting up a univariate
feature pre-selection method as well as a multivariate feature selection method
(frequency-based or ensemble feature selection) to select a panel of
differential features. When an output path is defined (via \code{output.path})
results will be saved on the hard disk and when \code{verbose} is TRUE
additional information will be printed to the console.

Frequency-based feature selection (\code{method="frequency"}): The whole data is
splitted in k cross validation training and test set pairs. For each training
set a multivariate feature selection procedure is performed. The resulting k
feature subsets are tested using the corresponding test sets (via
classification). As a result, \code{selectFeatures()} returns the average k-fold
cross validation classification accuracy as well as the selected feature panel
(i.e., the union set of the k particular feature subsets). As multivariate
feature selection methods random forest recursive feature elimination (RF-RFE),
random jungle recursive feature elimination (RJ-RFE) and support vector machine
recursive feature elimination (SVM-RFE) are supported. To reduce running times,
optionally, univariate feature pre-selection can be performed (usage via
\code{preselection.method}). As univariate pre-selection methods mMs
\code{"mMs"}, Student's t-test \code{"tTest"} and mRMR \code{"mrmr"} are
supported. Alternatively, no pre-selection can be chosen \code{"none"}. This
approach is similar to the method proposed in Baek et al.

Ensemble feature selection (\code{method="ensemble"}): From the whole data the
previously defined number of subsamples is drawn defining pairs of training and
test sets. Moreover, for each training set a previously defined number of
bootstrap samples is drawn. Then, for each bootstrap sample SVM-RFE is performed
and a feature ranking is obtained. To obtain a final ranking for a particular
training set, all associated bootstrap rankings are aggregated to a single
ranking. To score the \code{cutoff} best features, for each subsample a
classification of the test set is performed (using a svm trained with the
\code{cutoff} best features from the training set) and the classification
accuracy is determined. Finally, the stability of the subsample-specific panels
is assessed (via Kuncheva index, Kuncheva LI, 2007), all subsample-specific
rankings are aggregated, the top n features (defined by \code{cutoff}) are
selected, the average classification accuracy is computed, and all these results
are returned in a list. This approach has been proposed and is described in
Abeel et al.        
}



\value{
If \code{method} is \code{"frequency"}, the results list contains the following
elements:
  \item{accuracy}{ average k-fold cross validation accuracy. }
  \item{sensitivity}{ average k-fold cross validation sensitivity. }
  \item{specificity}{ average k-fold cross validation specificity. }
  \item{features}{ selected feature panel. }
  \item{all.results}{ complete cross validation results. }

If \code{method} is \code{"ensemble"}, the results list contains the following
elements:
  \item{accuracy}{ average accuracy regarding all subsamples. }
  \item{sensitivity}{ average sensitivity regarding all subsamples. }
  \item{specificity}{ average specificity regarding all subsamples. }
  \item{features}{ selected feature panel. }
  \item{all.results}{ all feature ranking results. }
  \item{stability}{ stability of the feature panel (i.e., Kuncheva index for the
  subrun-specific panels). }
}



\references{
Baek S, Tsai CA, Chen JJ.: Development of biomarker classifiers from high-
dimensional data. Brief Bioinform. 2009 Sep;10(5):537-46.

Abeel T, Helleputte T, Van de Peer Y, Dupont P, Saeys Y: Robust biomarker
identification for cancer diagnosis with ensemble feature selection methods.
Bioinformatics. 2010 Feb 1;26(3):392-8.

Kuncheva, LI: A stability index for feature selection. Proceedings of the IASTED
International Conference on Artificial Intelligence and Applications.
February 12-14, 2007. Pages: 390-395.
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

pre.sel.results <- preselect(elist=elist, columns1=c1, columns2=c2, label1="AD",
 label2="NDC", discard.threshold=0.1, fold.thresh=2, discard.features=TRUE,
 method="tTest")
elist <- elist[-pre.sel.results$discard,]

selectFeatures.results <- selectFeatures(elist,n1=20,n2=20,label1="AD",
 label2="NDC",selection.method="rf.rfe",preselection.method="none",subruns=2,
 k=2,candidate.number=20,method="frequency")
}



\keyword{ Feature selection }