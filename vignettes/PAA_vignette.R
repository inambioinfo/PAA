### R code from vignette source 'C:/R/R-3.2.4revised/library/PAA_DONT_DELETE/vignettes/PAA_vignette.Rnw'

###################################################
### code chunk number 1: style
###################################################
BiocStyle::latex()


###################################################
### code chunk number 2: install (eval = FALSE)
###################################################
## # only if you install a Bioconductor package for the first time
## source("http://www.bioconductor.org/biocLite.R")
## # else
## library("BiocInstaller")
## biocLite("PAA", dependencies=TRUE)


###################################################
### code chunk number 3: start
###################################################
library(PAA)


###################################################
### code chunk number 4: targets
###################################################
targets <- read.table(file=list.files(system.file("extdata", package="PAA"),
 pattern = "^targets", full.names = TRUE), header=TRUE)
print(targets[1:3,])


###################################################
### code chunk number 5: loadGPR
###################################################
gpr <- system.file("extdata", package="PAA")
targets <- list.files(system.file("extdata", package="PAA"),
   pattern = "dummy_targets", full.names=TRUE)
dummy.elist <- loadGPR(gpr.path=gpr, targets.path=targets,
   array.type="ProtoArray")
save(dummy.elist, file=paste(gpr, "/DummyData.RData",
   sep=""), compress="xz")


###################################################
### code chunk number 6: loadGPR
###################################################
targets2 <- list.files(system.file("extdata", package="PAA"),
 pattern = "dummy_no_descr_targets", full.names=TRUE)
elist2 <- loadGPR(gpr.path=gpr, targets.path=targets2, array.type="other",
 description="Name", description.features="^Hs~", description.discard="Empty")


###################################################
### code chunk number 7: load
###################################################
cwd <- system.file(package="PAA")
dir.create(paste(cwd, "/demo/demo_output", sep=""))
output.path <- paste(cwd, "/demo/demo_output",  sep="")
load(paste(cwd, "/extdata/Alzheimer.RData", sep=""))


###################################################
### code chunk number 8: plotArray1
###################################################
plotArray(elist=elist, idx=3, data.type="bg", log=FALSE, normalized=FALSE,
  aggregation="min", colpal="topo.colors")


###################################################
### code chunk number 9: plotArray2
###################################################
plotArray(elist=elist, idx=3, data.type="fg", log=FALSE, normalized=FALSE,
  aggregation="min", colpal="topo.colors")


###################################################
### code chunk number 10: backgroundCorrect
###################################################
library(limma)
elist <- backgroundCorrect(elist, method="normexp",
 normexp.method="saddle")


###################################################
### code chunk number 11: batchFilter
###################################################
lot1 <- elist$targets[elist$targets$Batch=='Batch1','ArrayID']
lot2 <- elist$targets[elist$targets$Batch=='Batch2','ArrayID']
elist.bF <- batchFilter(elist=elist, lot1=lot1, lot2=lot2, log=FALSE,
 p.thresh=0.001, fold.thresh=3)


###################################################
### code chunk number 12: batchFilterAnova
###################################################
elist.bF.a <- batchFilter.anova(elist=elist, log=FALSE, p.thresh=0.001,
 fold.thresh=3)
elist <- elist.bF


###################################################
### code chunk number 13: plotNormMethods (eval = FALSE)
###################################################
## plotNormMethods(elist=elist)


###################################################
### code chunk number 14: plotMAPlots
###################################################
plotMAPlots(elist=elist, idx=10)


###################################################
### code chunk number 15: normalizeArrays
###################################################
elist <- normalizeArrays(elist=elist, method="cyclicloess",
cyclicloess.method="fast")


###################################################
### code chunk number 16: batchAdjust
###################################################
elist <- batchAdjust(elist=elist, log=TRUE)


###################################################
### code chunk number 17: plotArray3
###################################################
plotArray(elist=elist, idx=3, data.type="fg", log=TRUE, normalized=TRUE,
  aggregation="min", colpal="topo.colors")


###################################################
### code chunk number 18: unlog
###################################################
elist.unlog <- elist
elist.unlog$E <- 2^(elist$E)


###################################################
### code chunk number 19: volcanoPlot1
###################################################
c1 <- paste(rep("AD",20), 1:20, sep="")
c2 <- paste(rep("NDC",20), 1:20, sep="")
#volcanoPlot(elist=elist.unlog, group1=c1, group2=c2, method="tTest",
volcanoPlot(elist=elist, group1=c1, group2=c2, log=TRUE, method="tTest",
p.thresh=0.01, fold.thresh=2)


###################################################
### code chunk number 20: volcanoPlot2 (eval = FALSE)
###################################################
## mMs.matrix1 <- mMs.matrix2 <- mMsMatrix(x=20, y=20)
## volcanoPlot(elist=elist.unlog, group1=c1, group2=c2, log=FALSE, method="mMs",
## p.thresh=0.01, fold.thresh=2, mMs.matrix1=mMs.matrix1,
## mMs.matrix2=mMs.matrix2, above=1500, between=400)


###################################################
### code chunk number 21: pvaluePlot1
###################################################
pvaluePlot(elist=elist, group1=c1, group2=c2, log=TRUE, method="tTest")


###################################################
### code chunk number 22: pvaluePlot2 (eval = FALSE)
###################################################
## mMs.matrix1 <- mMs.matrix2 <- mMsMatrix(x=20, y=20)
## pvaluePlot(elist=elist.unlog, group1=c1, group2=c2, log=FALSE, method="mMs",
## mMs.matrix1=mMs.matrix1, mMs.matrix2=mMs.matrix2, above=1500,
## between=400)


###################################################
### code chunk number 23: pvaluePlot3
###################################################
pvaluePlot(elist=elist, group1=c1, group2=c2, log=TRUE, method="tTest",
 adjust=TRUE)


###################################################
### code chunk number 24: pvaluePlot4 (eval = FALSE)
###################################################
## pvaluePlot(elist=elist.unlog, group1=c1, group2=c2, log=FALSE, method="mMs",
## mMs.matrix1=mMs.matrix1, mMs.matrix2=mMs.matrix2, above=1500,
## between=400, adjust=TRUE)


###################################################
### code chunk number 25: diffAnalysis
###################################################
E <- elist.unlog$E
rownames(E) <- paste(elist.unlog$genes[,1], elist.unlog$genes[,3],
    elist.unlog$genes[,2])
write.table(x=cbind(rownames(E),E),
    file=paste(cwd,"/demo/demo_output/data.txt", sep=""), sep="\t", eol="\n",
    row.names=FALSE, quote=FALSE)
mMs.matrix1 <- mMs.matrix2 <- mMsMatrix(x=20, y=20)
diff.analysis.results <- diffAnalysis(input=E, label1=c1, label2=c2,
    class1="AD", class2="NDC", output.path=output.path,
    mMs.matrix1=mMs.matrix1, mMs.matrix2=mMs.matrix2, above=1500,
    between=400)
print(diff.analysis.results[1:10,])


###################################################
### code chunk number 26: preselect
###################################################
mMs.matrix1 <- mMs.matrix2 <- mMsMatrix(x=20, y=20)
pre.sel.results <- preselect(elist=elist.unlog, columns1=c1, columns2=c2,
    label1="AD", label2="NDC", log=FALSE, discard.threshold=0.5,
    fold.thresh=1.5, discard.features=TRUE, mMs.above=1500, mMs.between=400,
    mMs.matrix1=mMs.matrix1, mMs.matrix2=mMs.matrix2,
    method="mMs")
elist <- elist[-pre.sel.results$discard,]


###################################################
### code chunk number 27: selectFeatures1 (eval = FALSE)
###################################################
## selectFeatures.results <- selectFeatures(elist,n1=20,n2=20,label1="AD",
##     label2="NDC",log=TRUE,selection.method="rf.rfe",subruns=2,
##     candidate.number=1000,method="frequency")


###################################################
### code chunk number 28: selectFeatures2 (eval = FALSE)
###################################################
## selectFeatures.results <- selectFeatures(elist,n1=20,n2=20,label1="AD",
##     label2="NDC",log=TRUE,subsamples=10,bootstraps=10,method="ensemble")


###################################################
### code chunk number 29: loadSelectFeaturesResults
###################################################
# results of frequency-based feature selection:
load(paste(cwd, "/extdata/selectFeaturesResultsFreq.RData", sep=""))
# or results of ensemble feature selection:
load(paste(cwd, "/extdata/selectFeaturesResultsEns.RData", sep=""))


###################################################
### code chunk number 30: plotFeatures
###################################################
plotFeatures(features=selectFeatures.results$features, elist=elist, n1=20,
    n2=20, group1="AD", group2="NDC")


###################################################
### code chunk number 31: plotFeaturesHeatmap
###################################################
plotFeaturesHeatmap(features=selectFeatures.results$features, elist=elist,
    n1=20, n2=20, description=TRUE)


###################################################
### code chunk number 32: plotFeaturesHeatmap2
###################################################
plotFeaturesHeatmap.2(features=selectFeatures.results$features, elist=elist,
    n1=20, n2=20, description=TRUE)


###################################################
### code chunk number 33: printFeatures
###################################################
elist$E <- round(elist$E,2)
printFeatures(features=selectFeatures.results$features, elist=elist)[,-2]


