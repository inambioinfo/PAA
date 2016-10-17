test_preselect <- function() {
    cwd <- system.file("extdata", package="PAA")
    load(paste(cwd, "/Alzheimer.RData", sep=""))
    load(paste(cwd, "/selectFeaturesResultsFreq.RData", sep=""))

    test.results <- printFeatures(features=selectFeatures.results$features, elist=elist)



    checkTrue(is.data.frame(test.results))
    checkTrue(nrow(test.results) >= 1)
}