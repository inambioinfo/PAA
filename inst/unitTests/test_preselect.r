test_preselect <- function() {
    cwd <- system.file("extdata", package="PAA")
    load(paste(cwd, "/Alzheimer.RData", sep=""))

    elist <- elist[1:100,]
    mMs.matrix1 <- mMs.matrix2 <- mMsMatrix(x=20, y=20)
    c1 <- paste(rep("AD",20), 1:20, sep="")
    c2 <- paste(rep("NDC",20), 1:20, sep="")
    test.results <- preselect(elist=elist, columns1=c1, columns2=c2,
     label1="AD", label2="NDC", discard.threshold=0.5, fold.thresh=1.5,
     discard.features=TRUE, mMs.above=1500, mMs.between=400,
     mMs.matrix1=mMs.matrix1, mMs.matrix2=mMs.matrix2,
     method="mMs")
    test.elist <- elist[-test.results$discard,]



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