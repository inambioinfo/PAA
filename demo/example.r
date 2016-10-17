library(PAA)

#find R_HOME
cwd <- system.file(package="PAA")
dir.create(paste(cwd, "/demo/demo_output", sep=""))
output.path <- paste(cwd, "/demo/demo_output",  sep="")



#load exemplary raw data and save expression list object
#(limma class EListRaw) for further sessions - DON'T RUN!!!
#elist <- loadGPR(gpr.path=paste(cwd, "/path/to/gpr_files", sep=""),
# targets.path=paste(cwd, "/path/to/targets_file/targets_file.txt", sep=""),
# array.type="ProtoArray", aggregation="min")
#save(elist, file=paste(cwd, "/path/to/expression_list/EListRawObject.RData",
# sep=""))


#load expression lists
# --> loading 'elist'
load(paste(cwd, "/extdata/Alzheimer.RData", sep=""), verbose=TRUE)  
#shuffleData(elist=elist, n1=20, n2=20, label1="A", label2="B")



#diagnostic plots before pre-processing
c1 <- paste(rep("AD",20), 1:20, sep="")
c2 <- paste(rep("NDC",20), 1:20, sep="")
mMs.matrix1 <- mMs.matrix2 <- mMsMatrix(x=20, y=20) #>>>>>>>>> BECAUSE EQUAL!!!
volcanoPlot(elist=elist, group1=c1, group2=c2, log=FALSE, method="mMs",
 p.thresh=0.01, fold.thresh=2, tag="_raw_mMs", mMs.matrix1=mMs.matrix1,
 mMs.matrix2=mMs.matrix2, above=1500, between=400)
volcanoPlot(elist=elist, group1=c1, group2=c2, log=FALSE, method="mMs",
 p.thresh=0.01, fold.thresh=2, output.path=output.path, tag="_raw_mMs",
 mMs.matrix1=mMs.matrix1, mMs.matrix2=mMs.matrix2, above=1500, between=400)
volcanoPlot(elist=elist, group1=c1, group2=c2, log=FALSE, method="tTest",
 p.thresh=0.01, fold.thresh=2, tag="_raw_tTest")
volcanoPlot(elist=elist, group1=c1, group2=c2, log=FALSE, method="tTest",
 p.thresh=0.01, fold.thresh=2, output.path=output.path, tag="_raw_tTest")
pvaluePlot(elist=elist, group1=c1, group2=c2, log=FALSE, method="tTest",
 tag="_lots_raw_tTest")
pvaluePlot(elist=elist, group1=c1, group2=c2, log=FALSE, method="tTest",
 output.path=output.path, tag="_lots_raw_tTest")
pvaluePlot(elist=elist, group1=c1, group2=c2, log=FALSE, method="tTest",
 tag="_lots_raw_tTest", adjust=TRUE)
pvaluePlot(elist=elist, group1=c1, group2=c2, log=FALSE, method="tTest",
 output.path=output.path, tag="_lots_raw_tTest", adjust=TRUE)



#pre-processing
lot1 <- elist$targets[elist$targets$Batch=='Batch1','ArrayID']
lot2 <- elist$targets[elist$targets$Batch=='Batch2','ArrayID']
elist <- batchFilter(elist=elist, lot1=lot1, lot2=lot2, log=FALSE,
 p.thresh=0.001, fold.thresh=3, output.path=output.path)
elist <- limma:::backgroundCorrect(elist, method="normexp",
 normexp.method="saddle")
plotNormMethods(elist=elist, include.rlm=TRUE)
plotNormMethods(elist=elist, include.rlm=TRUE,
  output.path=paste(cwd,"/demo/demo_output",sep=""))
plotMAPlots(elist=elist, idx=10, include.rlm=FALSE)
plotMAPlots(elist=elist, idx="all", include.rlm=TRUE,
 output.path=paste(cwd,"/demo/demo_output",sep=""))
elist <- normalizeArrays(elist=elist, method="cyclicloess",
 cyclicloess.method="pairs")
elist <- batchAdjust(elist=elist, log=TRUE)
elist.unlog <- elist
elist.unlog$E <- 2^(elist$E)



#diagnostic plots after pre-processing
volcanoPlot(elist=elist, group1=c1, group2=c2, log=TRUE, method="mMs",
 p.thresh=0.01, fold.thresh=2, tag="_normalized_mMs", mMs.matrix1=mMs.matrix1,
 mMs.matrix2=mMs.matrix2, above=1500, between=400)
volcanoPlot(elist=elist, group1=c1, group2=c2, log=TRUE, method="mMs",
 p.thresh=0.01, fold.thresh=2, output.path=output.path, tag="_normalized_mMs",
 mMs.matrix1=mMs.matrix1, mMs.matrix2=mMs.matrix2, above=1500, between=400)
volcanoPlot(elist=elist, group1=c1, group2=c2, log=TRUE, method="tTest",
 p.thresh=0.01, fold.thresh=2, tag="_normalized_tTest")
volcanoPlot(elist=elist, group1=c1, group2=c2, log=TRUE, method="tTest",
 p.thresh=0.01, fold.thresh=2, output.path=output.path, tag="_normalized_tTest")
pvaluePlot(elist=elist, group1=c1, group2=c2, log=TRUE, method="tTest",
 tag="_normalized_tTest")
pvaluePlot(elist=elist, group1=c1, group2=c2, log=TRUE, method="tTest",
 output.path=output.path, tag="_normalized_tTest")
pvaluePlot(elist=elist, group1=c1, group2=c2, log=TRUE, method="mMs",
 tag="_normalized_mMs", mMs.matrix1=mMs.matrix1, mMs.matrix2=mMs.matrix2,
 above=1500, between=400)
pvaluePlot(elist=elist, group1=c1, group2=c2, log=TRUE, method="mMs",
 output.path=output.path, tag="_normalized_mMs", mMs.matrix1=mMs.matrix1,
 mMs.matrix2=mMs.matrix2, above=1500, between=400)
pvaluePlot(elist=elist, group1=c1, group2=c2, log=TRUE, method="tTest",
 tag="_normalized_tTest", adjust=TRUE)
pvaluePlot(elist=elist, group1=c1, group2=c2, log=TRUE, method="tTest",
 output.path=output.path, tag="_normalized_tTest", adjust=TRUE)
pvaluePlot(elist=elist, group1=c1, group2=c2, log=TRUE, method="mMs",
 tag="_normalized_mMs", mMs.matrix1=mMs.matrix1, mMs.matrix2=mMs.matrix2,
 above=1500, between=400, adjust=TRUE)
pvaluePlot(elist=elist, group1=c1, group2=c2, log=TRUE, method="mMs",
 output.path=output.path, tag="_normalized_mMs", mMs.matrix1=mMs.matrix1,
 mMs.matrix2=mMs.matrix2, above=1500, between=400, adjust=TRUE)



#differential analysis
output <- elist.unlog$E
rownames(output) <- paste(elist.unlog$genes$Block, elist.unlog$genes$Row,
  elist.unlog$genes$Column)
write.table(x=cbind(rownames(output),output),
 file=paste(cwd, "/demo/demo_output/data.txt",  sep=""), sep="\t", eol="\n",
 row.names=FALSE, quote=FALSE)
diffAnalysis(input=output, label1=c1, label2=c2, class1="AD", class2="NDC",
 mMs.matrix1=mMs.matrix1, mMs.matrix2=mMs.matrix2, above=1500, between=400)
diff.analysis.results <- diffAnalysis(input=output, label1=c1, label2=c2,
 class1="AD", class2="NDC", output.path=output.path, mMs.matrix1=mMs.matrix1,
 mMs.matrix2=mMs.matrix2, above=1500, between=400)
print(diff.analysis.results[1:10,])



#feature preselection
pre.sel.results <- preselect(elist=elist.unlog, columns1=c1, columns2=c2,
 label1="AD", label2="NDC", log=FALSE, discard.threshold=0.5, fold.thresh=1.5,
 discard.features=TRUE, mMs.above=1500, mMs.between=400,
 mMs.matrix1=mMs.matrix1, mMs.matrix2=mMs.matrix2, method="mMs")
elist <- elist[-pre.sel.results$discard,]



#feature selection and final classification - note:
#subruns should be set to 100 (here: subruns=10 beacaue running times would be
#too long for this example)
row.number <- nrow(elist$E)
selectFeatures.results <- selectFeatures(elist,n1=20,n2=20,label1="AD",
 label2="NDC",log=TRUE,selection.method="rf.rfe",subruns=2,
 candidate.number=1000,method="frequency")
selectFeatures.results <- selectFeatures(elist,n1=20,n2=20,label1="AD",
 label2="NDC",log=TRUE,subsamples=10,bootstraps=10,method="ensemble")

plotFeatures(features=selectFeatures.results$features, elist=elist, n1=20,
 n2=20, group1="AD", group2="NDC")
plotFeatures(features=selectFeatures.results$features, elist=elist, n1=20,
 n2=20, group1="AD", group2="NDC", output.path=paste(cwd, "/demo/demo_output",
 sep=""))
 
plotFeaturesHeatmap(features=selectFeatures.results$features, elist=elist,
 n1=20, n2=20, description=TRUE)
plotFeaturesHeatmap(features=selectFeatures.results$features, elist=elist,
 n1=20, n2=20, output.path=paste(cwd, "/demo/demo_output", sep=""),
 description=TRUE)
 
plotFeaturesHeatmap.2(features=selectFeatures.results$features, elist=elist,
 n1=20, n2=20, description=TRUE)
plotFeaturesHeatmap.2(features=selectFeatures.results$features, elist=elist,
 n1=20, n2=20, output.path=paste(cwd, "/demo/demo_output", sep=""),
 description=TRUE)
 
printFeatures(features=selectFeatures.results$features, elist=elist.unlog)
printFeatures(features=selectFeatures.results$features, elist=elist.unlog,
 output.path=paste(cwd, "/demo/demo_output",  sep=""))