test_loadGPR <- function() {
    gpr <- system.file("extdata", package="PAA")
    targets <- list.files(system.file("extdata", package="PAA"),
     pattern = "dummy_targets", full.names=TRUE)
    test.elist <- loadGPR(gpr.path=gpr, targets.path=targets)



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