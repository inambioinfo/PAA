test_batchAdjust <- function() {
    cwd <- system.file("extdata", package="PAA")
    load(paste(cwd, "/Alzheimer.RData", sep=""))

    elist <- elist[1:100,]
    test.elist <- batchAdjust(elist=elist, log=TRUE)



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