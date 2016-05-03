test_diffAnalysis <- function() {
    cwd <- system.file("extdata", package="PAA")
    load(paste(cwd, "/Alzheimer.RData", sep=""))

    elist <- elist[1:100,]
    output <- elist$E
    rownames(output) <- paste(elist$genes[,1], " ", elist$genes[,3], " ",
     elist$genes[,2], sep="")
    mMs.matrix1 <- mMs.matrix2 <- mMsMatrix(x=20, y=20)
    c1 <- paste(rep("AD",20), 1:20, sep="")
    c2 <- paste(rep("NDC",20), 1:20, sep="")
    test.results <- diffAnalysis(input=output, label1=c1, label2=c2,
     class1="AD", class2="NDC", output.path=NULL,
     mMs.matrix1=mMs.matrix1, mMs.matrix2=mMs.matrix2, above=1500,
     between=400)



    checkTrue(is.data.frame(test.results))
    checkTrue(nrow(test.results) >= 1)
}