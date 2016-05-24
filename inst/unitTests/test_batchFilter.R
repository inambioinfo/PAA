test_batchFilter <- function() {
    cwd <- system.file("extdata", package="PAA")
    load(paste(cwd, "/Alzheimer.RData", sep=""))
    
    lot1 <- elist$targets[elist$targets$Batch=='Batch1','ArrayID']
    lot2 <- elist$targets[elist$targets$Batch=='Batch2','ArrayID']
    elist <- elist[1:100,]
    test.elist <- batchFilter(elist=elist, lot1=lot1, lot2=lot2, log=FALSE,
     p.thresh=0.001, fold.thresh=3)



    checkTrue(is.matrix(test.elist$E))
    checkTrue(nrow(test.elist$E) >= 1)
    checkTrue(is.numeric(test.elist$E[1,1]))
    checkTrue(is.matrix(test.elist$Eb))
    checkTrue(nrow(test.elist$Eb) >= 1)
    checkTrue(is.numeric(test.elist$Eb[1,1]))
    checkTrue(is.data.frame(test.elist$targets))
    checkTrue(nrow(test.elist$targets) >= 1)
    checkTrue(is.data.frame(test.elist$genes))
    checkTrue(nrow(test.elist$genes) >= 1)
    checkTrue(is.character(test.elist$source))
    checkTrue(is.list(test.elist$printer))
}