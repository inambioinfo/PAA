##########################################################################
#                      ProteinArrayAnalyzer (PAA)                        #
##########################################################################

##########################################################################
#                            FUNCTIONS                                   #
##########################################################################

#+++++++++++++++++++++++++++ substitute +++++++++++++++++++++++++++++++++++++++
substitute <- function(brc=NULL) {
    if(is.null(brc)) {
        stop("ERROR: Not all mandatory arguments have been defined!")
    }
    tmp <- gsub('(.*)', 'x\\1x', brc)
    tmp <- gsub(' ', '_', tmp)
    return(tmp)
}
#+++++++++++++++++++++++++++ substitute +++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++ resubstitute ++++++++++++++++++++++++++++++++++++++
resubstitute <- function(xb_r_cx=NULL) {
    if(is.null(xb_r_cx)) {
        stop("ERROR: Not all mandatory arguments have been defined!")
    }
    tmp <- gsub('_', ' ', xb_r_cx)
    tmp <- gsub('x', '', tmp)
    return(tmp)
}
#+++++++++++++++++++++++++++ resubstitute ++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++ svm.rfe +++++++++++++++++++++++++++++++++++++++
svm.rfe <- function(iteration=NULL, data.table=NULL, output.path=NULL,
  label1="A", label2="B", n1=NULL, n2=NULL, panel.selection.criterion="accuracy",
  verbose=FALSE){
    if(is.null(iteration) || is.null(data.table) || is.null(n1) ||
    is.null(n2)) {
        stop("ERROR: Not all mandatory arguments have been defined!")
    }
    
    if(min(n1,n2) < 10){
        k <- min(n1,n2) #the 'k' for the classifier k-fold cross validation
    }else{
        k <- 10   #the 'k' for the classifier k-fold cross validation
    }
    brc <- substitute(rownames(data.table))
    data.table <- t(data.table)
    colnames(data.table) <- brc
    trainsample_urn <- rownames(data.table)
    testsample_urn1 <- rownames(data.table[1:n1,])
    testsample_urn2 <- rownames(data.table[{n1+1}:{n1+n2},])
    testsamplesize1 <- floor(n1/k)
    testsamplesize2 <- floor(n2/k)
    testsamplesize <- testsamplesize1+testsamplesize2
    trainsamplesize <- n1+n2-testsamplesize1-testsamplesize2
    classes_test <- factor(c(rep(label1, testsamplesize1),
      rep(label2, testsamplesize2)))
    classes_train <- factor(c(rep(label1, {n1-testsamplesize1}),
      rep(label2, {n2-testsamplesize2})))
    
    p <- ncol(data.table)
    survivingFeaturesIndexes <- seq(1:p)
    
    acc.best.set.accuracy <- 0
    acc.best.set <- c()
    sens.best.set.sensitivity <- 0
    sens.best.set <- c()
    spec.best.set.specificity <- 0
    spec.best.set <- c()
    if(!is.null(output.path)){
        cat(x="set\tAccuracy\tSensitivity\tSpecificity\tTP\tFP\tTN\tFN\n",
          file=paste(output.path, "/rfe_", iteration, ".txt", sep=""),
          append=FALSE)
    }
    while({p <- length(survivingFeaturesIndexes)} > 0){
        rankingCriteria <- matrix()
        tmp_trainsample_urn <- trainsample_urn
        tmp_testsample_urn1 <- testsample_urn1
        tmp_testsample_urn2 <- testsample_urn2
        crossValTest <- matrix(nrow=k,ncol=testsamplesize)
        crossValTrain <- matrix(nrow=k,ncol=trainsamplesize)
        for (i in 1:k){
            testsample1 <- sample(tmp_testsample_urn1, testsamplesize1,
              replace=FALSE)
            testsample2 <- sample(tmp_testsample_urn2, testsamplesize2,
              replace=FALSE)
            tmp_testsample_urn1 <-
              tmp_testsample_urn1[!(tmp_testsample_urn1 %in% testsample1)]
            tmp_testsample_urn2 <-
              tmp_testsample_urn2[!(tmp_testsample_urn2 %in% testsample2)]
            crossValTest[i,] <- testsample <- c(testsample1, testsample2)
            crossValTrain[i,] <-
              tmp_trainsample_urn[!(tmp_trainsample_urn %in% testsample)]
        }
        TP <- 0
        FP <- 0
        TN <- 0
        FN <- 0
        for (i in 1:k){
            test.dat <- data.table[crossValTest[i,], survivingFeaturesIndexes,
              drop=FALSE]
            train.dat <- data.table[crossValTrain[i,], survivingFeaturesIndexes,
              drop=FALSE]

            #~~~~~~~~~~~~~CLASSIFIER BEGIN~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #svm.model <- svm(train.dat, classes_train, type="C-classification",
            #  kernel="radial")
            #svm.model <- svm(train.dat, classes_train, cost=10,
            #  type="C-classification", kernel="linear")
            svm.model <- svm(train.dat, classes_train, type="C-classification",
              kernel="linear")
            w <- t(svm.model$coefs)%*%svm.model$SV
            if(i > 1){
                rankingCriteria <- rankingCriteria + {w*w}
            }else{
                rankingCriteria <- {w*w}
            }
            svm.pred <- predict(object=svm.model, newdata=test.dat)
            confusion_matrix <- table(observed=factor(classes_test,
              levels=sort(c(label1, label2))), predicted=svm.pred)
            #~~~~~~~~~~~~~CLASSIFIER END~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            TP <- TP + confusion_matrix[1,1]
            FP <- FP + confusion_matrix[2,1]
            TN <- TN + confusion_matrix[2,2]
            FN <- FN + confusion_matrix[1,2]
        }
        ACCURACY <- {TP+TN}/{TP+FP+TN+FN}
        SENSITIVITY <- TP/{TP+FN}   #=TPR
        SPECIFICITY <- TN/{TN+FP}   #=(1-FPR)
        #FPR_total <- FP/{FP+TN}   #=(1-SPECIFITY)
        if(!is.null(output.path)){
            cat(x=paste(p, ACCURACY, SENSITIVITY, SPECIFICITY, TP, FP, TN, FN,
              sep="\t"), file=paste(output.path, "/rfe_", iteration, ".txt",
              sep=""), append=TRUE)
            cat(x="\n", file=paste(output.path, "/rfe_", iteration, ".txt",
              sep=""), append=TRUE)
        }
        if(verbose){
            message(paste("svm.rfe - features: ", p, ", ACCURACY: ", ACCURACY,
              ", SENSI: ", SENSITIVITY, ", SPECI: ", SPECIFICITY, sep=""), "\n")
        }

        if(acc.best.set.accuracy <= ACCURACY){
            acc.best.set.accuracy <- ACCURACY
            acc.best.set <- colnames(data.table[,survivingFeaturesIndexes,
              drop=FALSE])
        }
        if(sens.best.set.sensitivity <= SENSITIVITY){
            sens.best.set.sensitivity <- SENSITIVITY
            sens.best.set <- colnames(data.table[,survivingFeaturesIndexes,
              drop=FALSE])
        }
        if(spec.best.set.specificity <= SPECIFICITY){
            spec.best.set.specificity <- SPECIFICITY
            spec.best.set <- colnames(data.table[,survivingFeaturesIndexes,
              drop=FALSE])
        }
        
        if(p >= 10){
            tenth <- floor(p*0.1)
        }else if(p < 10 && p != 1){
            tenth <- ceiling(p*0.1)
        }else{
            break
        }
        ranking <- sort(rankingCriteria, index.return = TRUE)$ix
        if(tenth > 1){
            survivingFeaturesIndexes <-
              survivingFeaturesIndexes[-ranking[1:tenth]]
        }else if(tenth == 1){
            survivingFeaturesIndexes <- survivingFeaturesIndexes[-ranking[1]]
        }
        
    } 
    if(panel.selection.criterion == "accuracy"){
        if(verbose){
            message(paste("feature selection - optimal number of features: ",
              length(acc.best.set), sep=""), "\n")
            message(paste("feature selection - best accuracy: ",
              acc.best.set.accuracy, sep=""), "\n")
        }
        return(resubstitute(acc.best.set))
    }else if(panel.selection.criterion == "sensitivity"){
        if(verbose){
            message(paste("feature selection - optimal number of features: ",
              length(sens.best.set), sep=""), "\n")
            message(paste("feature selection - best sensitivity: ",
              sens.best.set.sensitivity, sep=""), "\n")
        }
        return(resubstitute(sens.best.set))
    }else if(panel.selection.criterion == "specificity"){
        if(verbose){
            message(paste("feature selection - optimal number of features: ",
              length(spec.best.set), sep=""), "\n")
            message(paste("feature selection - best specificity: ",
              spec.best.set.specificity, sep=""), "\n")
        }
        return(resubstitute(spec.best.set))
    }
}
#+++++++++++++++++++++++++++ svm.rfe +++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++ rj.rfe +++++++++++++++++++++++++++++++++++++++
rj.rfe <- function(iteration=NULL, data.table=NULL, output.path=NULL,
  label1="A", label2="B", n1=NULL,
  n2=NULL, panel.selection.criterion="accuracy", importance.measure=NULL,
  ntree=NULL, mtry=NULL, verbose=FALSE){
    
    if(is.null(iteration) || is.null(data.table) || is.null(n1) ||
      is.null(n2)) {
        stop("ERROR: Not all mandatory arguments have been defined!")
    }
    
    if(is.null(output.path)){
        tmp.path <- tempdir()
        if(Sys.info()["sysname"] == "Windows"){
            rjungleInFile <- file.path(paste("\"", tmp.path,
              "/rjungleInFile.dat", "\"", sep=""))
        }else{
            rjungleInFile <-
              file.path(paste(tmp.path, "/rjungleInFile.dat", sep=""))
        }
        rjungleInFileR <-
          file.path(paste(tmp.path, "/rjungleInFile.dat", sep=""))
        if(Sys.info()["sysname"] == "Windows"){
            rjungleOutFile <-
              file.path(paste("\"", tmp.path, "/rj", "\"", sep=""))
        }else{
            rjungleOutFile <- file.path(paste(tmp.path, "/rj", sep=""))
        }
        rjungleGiniImportanceFile <- paste(tmp.path, "/rj.importance", sep="")
        rjunglePermImportanceFile <- paste(tmp.path, "/rj.importance2", sep="")
        rjungleConfusionFile <- paste(tmp.path, "/rj.confusion", sep="")
        #rjungleClassConfusionFile <- paste(tmp.path, "/rj.confusion2", sep="")
    } else {
        if(Sys.info()["sysname"] == "Windows"){
            rjungleInFile <-
              file.path(paste("\"", output.path, "/rjungleInFile.dat", "\"",
              sep=""))
        }else{
            rjungleInFile <-
              file.path(paste(output.path, "/rjungleInFile.dat", sep=""))
        }
        rjungleInFileR <-
          file.path(paste(output.path, "/rjungleInFile.dat", sep=""))
        if(Sys.info()["sysname"] == "Windows"){
            rjungleOutFile <-
              file.path(paste("\"", output.path, "/rj", "\"", sep=""))
        }else{
            rjungleOutFile <- file.path(paste(output.path, "/rj", sep=""))
        }
        rjungleGiniImportanceFile <-
          paste(output.path, "/rj.importance", sep="")
        rjunglePermImportanceFile <-
          paste(output.path, "/rj.importance2", sep="")
        rjungleConfusionFile <-
        paste(output.path, "/rj.confusion", sep="")
        #rjungleClassConfusionFile <-
        #paste(output.path, "/rj.confusion2", sep="")
    }

    brc <- substitute(rownames(data.table))
    data_table <- t(data.table)
    colnames(data_table) <- brc
    group <- c(rep(1, n1), rep(2,n2))
    data_table <- cbind(data_table, group)
    colnames(data_table)[ncol(data_table)] <- "Group" 
    write.table(data_table, file = rjungleInFileR, row.names = FALSE,
      quote = FALSE)
    
    acc.best.set.accuracy <- 0
    acc.best.set <- c()
    sens.best.set.sensitivity <- 0
    sens.best.set <- c()
    spec.best.set.specificity <- 0
    spec.best.set <- c()
    if(!is.null(output.path)){
        cat(x="set\tAccuracy\tSensitivity\tSpecificity\tTP\tFP\tTN\tFN\n",
          file=paste(output.path, "/rfe_", iteration, ".txt", sep=""),
          append=FALSE)
    }
    while({p <- ncol(data_table)-1} > 0){
        rjungleCMD <- ""
        if(is.null(importance.measure)){
            importance.measure.tmp <- "-i3"
        }else if(importance.measure == "MDA"){
            importance.measure.tmp <- "-i3"
        }else if(importance.measure == "MDG"){
            importance.measure.tmp <- "-i1"    
        }else{
            importance.measure.tmp <- "-i3"
        }
        
        if(is.null(mtry)){
            mtry.tmp <- paste("-m", sqrt(p), sep="")
        }else if(is.numeric(mtry) && mtry < 10){
            mtry.tmp <- paste("-m", mtry*sqrt(p), sep="")   
        }else{
            mtry.tmp <- paste("-m", mtry, sep="") 
        }
        
        if(is.null(ntree)){
            ntree.tmp <- "-t500"
        }else if(is.numeric(ntree) && ntree < 10){
            ntree.tmp <- paste("-t", max(ntree*p, 500), sep="")
        }else{
            ntree.tmp <- paste("-t", ntree, sep="")
        }
        
        if(Sys.info()["sysname"] == "Windows"){
            rjungleInFile <- gsub('^\"([A-Za-z]):', '\"/cygdrive/\\L\\1',
              rjungleInFile, perl=TRUE)
            rjungleOutFile <- gsub('^\"([A-Za-z]):', '\"/cygdrive/\\L\\1',
              rjungleOutFile, perl=TRUE)
        }
        if(p > 1000){
            rjungleCMD <- paste("rjungle -f", rjungleInFile,
              importance.measure.tmp, mtry.tmp, ntree.tmp, "-U4", "-D Group",
              "-o", rjungleOutFile)
            #rjungleCMD <- paste("rjungle -f", rjungleInFile, "-v",
            #  importance.measure.tmp, mtry.tmp, ntree.tmp, "-U2", "-D Group",
            #  "-o", rjungleOutFile)
        }else{
            rjungleCMD <- paste("rjungle -f", rjungleInFile,
              importance.measure.tmp, mtry.tmp, ntree.tmp, "-U1", "-D Group",
              "-o", rjungleOutFile)
            #rjungleCMD <- paste("rjungle -f", rjungleInFile, "-v",
            #  importance.measure.tmp, mtry.tmp, ntree.tmp, "-U1", "-D Group",
            # "-o", rjungleOutFile)
        }
        try(system(rjungleCMD, wait=TRUE))
                
            
        importanceGini <- read.table(rjungleGiniImportanceFile)
        importancePerm <- read.table(rjunglePermImportanceFile)
        confusion <- read.table(rjungleConfusionFile, skip=4, sep="\t")
            
        #Test/OOB set: (real outcome == rows / predicted outcome == columns )
        TP <- confusion[2,2]
        FP <- confusion[3,2]
        TN <- confusion[3,3]
        FN <- confusion[2,3]
        ACCURACY <- {TP+TN}/{TP+FP+TN+FN}
        SENSITIVITY <- TP/{TP+FN}   #=TPR
        SPECIFICITY <- TN/{TN+FP}   #=(1-FPR)
        
        if(!is.null(output.path)){    
            cat(x=paste(p, ACCURACY, SENSITIVITY, SPECIFICITY, TP, FP, TN, FN,
              sep="\t"), file=paste(output.path, "/rfe_", iteration, ".txt",
              sep=""), append=TRUE)
            cat(x="\n", file=paste(output.path, "/rfe_", iteration, ".txt",
              sep=""), append=TRUE)
        }
        if(verbose){
            message(paste("rj.rfe - features: ", p, ", ACCURACY: ", ACCURACY,
              ", SENSI: ", SENSITIVITY, ", SPECI: ", SPECIFICITY, sep=""), "\n")
        }
        
        #rf.importance <- data.frame(importancePerm[-1,3,drop=FALSE],
        #  importancePerm[-1,4,drop=FALSE],importanceGini[-1,4,drop=FALSE])
        #names(rf.importance) <- c('BRC', 'MDA', 'MDG')
        #write.table(x=rf.importance, file=paste(output.path, "/subrun",
        #  iteration, "/feature.elimination.importance_", p, ".txt", sep=""),
        #  sep="\t", eol="\n", row.names=FALSE)

        if(acc.best.set.accuracy <= ACCURACY){
            acc.best.set.accuracy <- ACCURACY
            acc.best.set <-
              colnames(data_table[,1:ncol(data_table)-1,drop=FALSE])
        }
        if(sens.best.set.sensitivity <= SENSITIVITY){
            sens.best.set.sensitivity <- SENSITIVITY
            sens.best.set <-
              colnames(data_table[,1:ncol(data_table)-1,drop=FALSE])
        }
        if(spec.best.set.specificity <= SPECIFICITY){
            spec.best.set.specificity <- SPECIFICITY
            spec.best.set <-
              colnames(data_table[,1:ncol(data_table)-1,drop=FALSE])
        }
        
        if(p >= 10){
            tenth <- floor(p*0.1)
        }else if(p < 10 && p != 1){
            tenth <- ceiling(p*0.1)
        }else{
            break
        }
        if(importance.measure == "MDG"){
            data_table <- data_table[,c(as.character(
              importanceGini[{tenth+2}:{ncol(data_table)},3]), "Group")]
        }else{
            data_table <- data_table[,c(as.character(
              importancePerm[{tenth+2}:{ncol(data_table)},3]), "Group")]
        }
        write.table(data_table, file = rjungleInFileR, row.names = FALSE,
          quote = FALSE)
    }
    if(panel.selection.criterion == "accuracy"){
        if(verbose){
            message(paste("feature selection - optimal number of features: ",
              length(acc.best.set), sep=""), "\n")
            message(paste("feature selection - best accuracy: ",
              acc.best.set.accuracy, sep=""), "\n")
        }
        return(resubstitute(acc.best.set))
    }else if(panel.selection.criterion == "sensitivity"){
        if(verbose){
            message(paste("feature selection - optimal number of features: ",
              length(sens.best.set), sep=""), "\n")
            message(paste("feature selection - best sensitivity: ",
              sens.best.set.sensitivity, sep=""), "\n")
        }
        return(resubstitute(sens.best.set))
    }else if(panel.selection.criterion == "specificity"){
        if(verbose){
            message(paste("feature selection - optimal number of features: ",
              length(spec.best.set), sep=""), "\n")
            message(paste("feature selection - best specificity: ",
              spec.best.set.specificity, sep=""), "\n")
        }
        return(resubstitute(spec.best.set))
    }
}
#+++++++++++++++++++++++++++ rj.rfe +++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++ rf.rfe +++++++++++++++++++++++++++++++++++++++
rf.rfe <- function(iteration=NULL, datamatrix=NULL, output.path=NULL,
  label1="A", label2="B", n1=NULL, n2=NULL,
  panel.selection.criterion="accuracy", importance.measure="MDA", ntree=500,
  mtry=NULL, verbose=FALSE){
    if(is.null(iteration) || is.null(datamatrix) || is.null(n1) ||
      is.null(n2)) {
        stop("ERROR: Not all mandatory arguments have been defined!")
    }
    repeats <- 10
    classes <- factor(c(rep(label1, n1), rep(label2, n2)))
    acc.best.set.accuracy <- 0
    acc.best.set <- c()
    sens.best.set.sensitivity <- 0
    sens.best.set <- c()
    spec.best.set.specificity <- 0
    spec.best.set <- c()
    if(!is.null(output.path)){ 
        cat(x="set\tAccuracy\tSensitivity\tSpecificity\tTP\tFP\tTN\tFN\n",
          file=paste(output.path, "/rfe_", iteration, ".txt", sep=""),
          append=FALSE)
    }
    
    if(is.null(ntree)){
        ntree <- 500
    }
    
    while({p <- nrow(datamatrix)} > 0){
        if(is.null(mtry)){
            mtry.tmp <- sqrt(p)
        }else{
            mtry.tmp <- mtry
        }
        TP <- 0
        FP <- 0
        TN <- 0
        FN <- 0
        rf.importance <- data.frame(resubstitute(rownames(datamatrix)),
          rep(0,nrow(datamatrix)),rep(0,nrow(datamatrix)))
        names(rf.importance) <- c('BRC', 'MDA', 'MDG')
        dat <- t(datamatrix)
        for (j in 1:repeats){
            #~~~~~~~~~~~~~~~CLASSIFIER BEGIN~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            model_randomForest <- randomForest(x=dat, y=classes,
              importance=TRUE, keep.forest=FALSE, ntree=ntree, mtry=mtry.tmp)

            importance <- data.frame(model_randomForest$importance)
            rf.importance.tmp <- data.frame(rownames(importance),
              importance[,3], importance[,4])
            names(rf.importance.tmp) <- c('BRC', 'MDA', 'MDG')
            rf.importance.tmp$BRC <- resubstitute(rf.importance.tmp$BRC)
            rf.importance$MDA <- rf.importance$MDA + rf.importance.tmp$MDA
            rf.importance$MDG <- rf.importance$MDG + rf.importance.tmp$MDG
            #~~~~~~~~~~~~~~~CLASSIFIER END~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            TP <- TP + model_randomForest$confusion[1,1]
            FP <- FP + model_randomForest$confusion[2,1]
            TN <- TN + model_randomForest$confusion[2,2]
            FN <- FN + model_randomForest$confusion[1,2]
        }
        ACCURACY <- {TP+TN}/{TP+FP+TN+FN}
        SENSITIVITY <- TP/{TP+FN}   #=TPR
        SPECIFICITY <- TN/{TN+FP}   #=(1-FPR)
        if(!is.null(output.path)){ 
            cat(x=paste(nrow(datamatrix), ACCURACY, SENSITIVITY, SPECIFICITY,
              TP, FP, TN, FN, sep="\t"), file=paste(output.path, "/rfe_",
              iteration, ".txt", sep=""), append=TRUE)
            cat(x="\n", file=paste(output.path, "/rfe_", iteration, ".txt",
              sep=""), append=TRUE)
        }
        if(verbose){
            message(paste("rf.rfe - features: ", p, ", ACCURACY: ", ACCURACY,
              ", SENSI: ", SENSITIVITY, ", SPECI: ", SPECIFICITY, sep=""), "\n")
        }

        if(acc.best.set.accuracy <= ACCURACY){
            acc.best.set.accuracy <- ACCURACY
            acc.best.set <- rownames(datamatrix)
        }
        if(sens.best.set.sensitivity <= SENSITIVITY){
            sens.best.set.sensitivity <- SENSITIVITY
            sens.best.set <- rownames(datamatrix)
        }
        if(spec.best.set.specificity <= SPECIFICITY){
            spec.best.set.specificity <- SPECIFICITY
            spec.best.set <- rownames(datamatrix)
        }

        rf.importance$MDA <- rf.importance$MDA / repeats
        rf.importance$MDG <- rf.importance$MDG / repeats
        if(importance.measure == "MDA"){
            rf.importance <- rf.importance[order(-rf.importance$MDA),]
        }else if(importance.measure == "MDG"){
            rf.importance <- rf.importance[order(-rf.importance$MDG),]
        }
        #write.table(x=rf.importance, file=paste(output.path, "/subrun",
        #  iteration, "/feature.elimination.importance_", nrow(datamatrix),
        #  ".txt", sep=""), sep="\t", eol="\n", row.names=FALSE)
        if(nrow(rf.importance) >= 10){
            next.set <- rf.importance$BRC[1:ceiling(nrow(importance)*0.9)]
        }else if(nrow(rf.importance) < 10 && nrow(rf.importance) != 1){
            next.set <- rf.importance$BRC[1:floor(nrow(importance)*0.9)]
        }else{
            break
        }
        next.set <-
          rownames(datamatrix)[resubstitute(rownames(datamatrix)) %in% next.set]
        datamatrix <- datamatrix[next.set,,drop=FALSE]
    }
    if(panel.selection.criterion == "accuracy"){
        if(verbose){
            message(paste("feature selection - optimal number of features: ",
              length(acc.best.set), sep=""), "\n")
            message(paste("feature selection - best accuracy: ",
              acc.best.set.accuracy, sep=""), "\n")
        }
        return(resubstitute(acc.best.set))
    }else if(panel.selection.criterion == "sensitivity"){
        if(verbose){
            message(paste("feature selection - optimal number of features: ",
              length(sens.best.set), sep=""), "\n")
            message(paste("feature selection - best sensitivity: ",
              sens.best.set.sensitivity, sep=""), "\n")
        }
        return(resubstitute(sens.best.set))
    }else if(panel.selection.criterion == "specificity"){
        if(verbose){
            message(paste("feature selection - optimal number of features: ",
              length(spec.best.set), sep=""), "\n")
            message(paste("feature selection - best specificity: ",
              spec.best.set.specificity, sep=""), "\n")
        }
        return(resubstitute(spec.best.set))
    }
}
#+++++++++++++++++++++++++++ rf.rfe +++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++ final.classify.rf +++++++++++++++++++++++++++++++++
final.classify.rf <- function(trainset=NULL, testset=NULL, classes_train=NULL,
  classes_test=NULL, label1="A", label2="B", cwd, iteration, plot=FALSE,
  performance=FALSE, verbose=FALSE){
    
    if(is.null(trainset) || is.null(testset) || is.null(classes_train) ||
      is.null(classes_test) || is.null(iteration)) {
        stop("ERROR: Not all mandatory arguments have been defined!")
    }
    
    TP_total <- 0
    FP_total <- 0
    TN_total <- 0
    FN_total <- 0

    train_dat1 <- t(trainset[,-1])
    test_dat1 <- t(testset[,-1]) 
    train_dat2 <- data.frame(train_dat1, classes_train)
    test_dat2 <- data.frame(test_dat1, classes_test)
    names(train_dat2) <- substitute(testset$BRC)
    names(test_dat2) <- names(train_dat2)
    names(train_dat2)[ncol(train_dat2)] <- "classes"
    names(test_dat2)[ncol(test_dat2)] <- "classes"

    if(performance){
        preds_total <- c()
        classes_total <- c()
    }
    for (i in 1:1){
        model_randomForest <- randomForest(classes ~ ., data=train_dat2,
          importance=TRUE, keep.forest=TRUE)
        pred_rf <- predict(object=model_randomForest, newdata=test_dat2,
          type="response", norm.votes=TRUE)
        confusion_matrix <- table(observed = factor(classes_test,
          levels=sort(c(label1, label2))), predicted = pred_rf)

        TP <- confusion_matrix[1,1]
        FP <- confusion_matrix[2,1]
        TN <- confusion_matrix[2,2]
        FN <- confusion_matrix[1,2]
        TP_total <- TP_total+TP
        FP_total <- FP_total+FP
        TN_total <- TN_total+TN
        FN_total <- FN_total+FN
        
        if(performance){
            pred_rf <- predict(object=model_randomForest, newdata=test_dat2,
              type="prob", norm.votes=TRUE)
            preds <- pred_rf[,2]
            perf <- performance(prediction(preds, classes_test), 'tpr', 'fpr')
            preds_total <- c(preds_total, preds)
            classes_total <- c(classes_total, classes_test)
        }
    }
    if(plot && !is.null(cwd)){
        perf <-
          performance(prediction(preds_total, classes_total), 'tpr', 'fpr')
        auc <- as.numeric(performance(prediction(preds_total, classes_total),
          'auc')@y.values)
        tiff(paste(cwd,"/roc_", iteration, ".tiff",sep=""), width=2000,
          height=2000, pointsize=15, compression="lzw", res=300)
            plot(perf, col="red", lwd=2, main=paste(label1, " vs. ", label2,
              " - average AUC: ", round(auc,4), sep=""))
        dev.off()
    }
    
    if(performance && !is.null(cwd)){
        cat(x=preds_total, file=paste(cwd, "/preds_total.txt", sep=""),
          append=TRUE, sep="\t")
        cat(x="\n", file=paste(cwd, "/preds_total.txt", sep=""), append=TRUE)
        cat(x=classes_total, file=paste(cwd, "/classes_total.txt", sep=""),
          append=TRUE, sep="\t")
        cat(x="\n", file=paste(cwd, "/classes_total.txt", sep=""), append=TRUE)
    }
    
    ACCURACY_total <- {TP_total+TN_total}/{TP_total+FP_total+TN_total+FN_total}
    SENSITIVITY_total <- TP_total/{TP_total+FN_total}   #=TPR
    SPECIFITY_total <- TN_total/{TN_total+FP_total}   #=(1-FPR)
    if(verbose){
        message(paste("test set classification - ACCURACY:", ACCURACY_total,
          ", SENSI: ", SENSITIVITY_total, ", SPECI: ", SPECIFITY_total, sep=""),
          "\n")
    }
    results <- list(accuracy=ACCURACY_total, sensitivity=SENSITIVITY_total,
      specificity=SPECIFITY_total, tp=TP_total, fp=FP_total, tn=TN_total,
      fn=FN_total)
    return(results)
}
#+++++++++++++++++++++++++++ final.classify.rf +++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++ final.classify.svm ++++++++++++++++++++++++++++++++
final.classify.svm <- function(trainset=NULL, testset=NULL, classes_train=NULL,
  classes_test=NULL, label1="A", label2="B", cwd, iteration, plot=FALSE,
  performance=FALSE, verbose=FALSE){
    
    if(is.null(trainset) || is.null(testset) || is.null(classes_train) ||
      is.null(classes_test) || is.null(iteration)) {
        stop("ERROR: Not all mandatory arguments have been defined!")
    }
    
    TP_total <- 0
    FP_total <- 0
    TN_total <- 0
    FN_total <- 0

    brc <- substitute(trainset$BRC)
    train.dat <- t(as.matrix(trainset[,-1]))
    test.dat <- t(as.matrix(testset[,-1])) 
    colnames(train.dat) <- brc
    colnames(test.dat) <- brc
    
    if(performance){
        preds_total <- c()
        classes_total <- c()
    }
    
    #tuneobj <-
    #  tune.svm(x=train.dat, y=classes_train, gamma=10^(-9:-1), cost=10^(0:4))
    tuneobj <- tune.svm(x=train.dat, y=classes_train, cost=10^(0:4))
    #bestGamma <- tuneobj$best.parameters[[1]]
    #bestCost <- tuneobj$best.parameters[[2]]
    bestCost <- tuneobj$best.parameters[[1]]
    #model.svm <- svm(train.dat, classes_train, probability=TRUE, cost=10,
    #  type="C-classification", kernel="linear")
    model.svm <- svm(train.dat, classes_train, probability=TRUE, cost=bestCost,
      type="C-classification", kernel="linear")
    #model.svm <- svm(train.dat, classes_train, probability=TRUE, cost=bestCost,
    #  gamma=bestGamma, type="C-classification", kernel="radial")
    pred.svm <- predict(object=model.svm, newdata=test.dat,
      decision.values = TRUE, probability = TRUE)
    confusion_matrix <- table(observed = factor(classes_test, levels=
      sort(c(label1, label2))), predicted = pred.svm)

    TP <- confusion_matrix[1,1]
    FP <- confusion_matrix[2,1]
    TN <- confusion_matrix[2,2]
    FN <- confusion_matrix[1,2]
    TP_total <- TP_total+TP
    FP_total <- FP_total+FP
    TN_total <- TN_total+TN
    FN_total <- FN_total+FN
        
    if(performance){
        perf <- performance(prediction(attributes(pred.svm)$decision.values,
          classes_test), 'tpr', 'fpr')
        preds_total <- c(preds_total, attributes(pred.svm)$decision.values)
        classes_total <- c(classes_total, classes_test)
    }
    
    if(plot && !is.null(cwd)){
        perf <-
          performance(prediction(preds_total, classes_total), 'tpr', 'fpr')
        auc <- as.numeric(performance(prediction(preds_total, classes_total),
          'auc')@y.values)
        tiff(paste(cwd,"/roc_", iteration, ".tiff",sep=""), width=2000,
          height=2000, pointsize=15, compression="lzw", res=300)
            plot(perf, col="red", lwd=2, main=paste(label1, " vs. ", label2,
            " - average AUC: ", round(auc,4), sep=""))
        dev.off()
    }
    
    if(performance && !is.null(cwd)){
        cat(x=preds_total, file=paste(cwd, "/preds_total.txt", sep=""),
          append=TRUE, sep="\t")
        cat(x="\n", file=paste(cwd, "/preds_total.txt", sep=""), append=TRUE)
        cat(x=classes_total, file=paste(cwd, "/classes_total.txt", sep=""),
          append=TRUE, sep="\t")
        cat(x="\n", file=paste(cwd, "/classes_total.txt", sep=""), append=TRUE)
    }
    
    ACCURACY_total <- {TP_total+TN_total}/{TP_total+FP_total+TN_total+FN_total}
    SENSITIVITY_total <- TP_total/{TP_total+FN_total}   #=TPR
    SPECIFITY_total <- TN_total/{TN_total+FP_total}   #=(1-FPR)
    if(verbose){
        message(paste("test set classification - ACCURACY:", ACCURACY_total,
          ", SENSI: ", SENSITIVITY_total, ", SPECI: ",
          SPECIFITY_total, sep=""), "\n")
    }
    results <- list(accuracy=ACCURACY_total, sensitivity=SENSITIVITY_total,
      specificity=SPECIFITY_total, tp=TP_total, fp=FP_total, tn=TN_total,
      fn=FN_total)
    return(results)
}
#+++++++++++++++++++++++++++ final.classify.svm ++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++ biomarkers.update +++++++++++++++++++++++++++++++++
biomarkers.update <- function(biomarker.list=NULL, update.vector=NULL){
    if(is.null(biomarker.list) || is.null(update.vector)) {
        stop("ERROR: Not all mandatory arguments have been defined!")
    }
    update.vector <- as.vector(update.vector)
    for(i in 1:length(update.vector)){
        if(update.vector[i] %in% names(biomarker.list)){
            biomarker.list[[update.vector[i]]] <-
              biomarker.list[[update.vector[i]]] + 1
        }else{
            biomarker.list[[update.vector[i]]] <- 1
        }
    }
    return(biomarker.list)
}
#+++++++++++++++++++++++++++ biomarkers.update +++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++ noFS +++++++++++++++++++++++++++++++++++++++
noFS <- function(elist=NULL, columns1=NULL, columns2=NULL,
  label1="A", label2="B"){
    
    if(is.null(elist) || is.null(columns1) || is.null(columns2)) {
        stop("ERROR: Not all mandatory arguments have been defined!")
    }
    
    row.number <- nrow(elist$E)
    n1 <- length(columns1)
    n2 <- length(columns2)

    testcols1 <- columns1[columns1 %in% sample(columns1,n1/3,replace=FALSE)]
    testcols2 <- columns2[columns2 %in% sample(columns2,n2/3,replace=FALSE)]
    testcols <- c(testcols1,testcols2)
    testcols1.length <- length(testcols1)
    testcols2.length <- length(testcols2)
    testcols.length <- length(testcols)
    testsample <- elist$E[,testcols]

    traincols1 <- columns1[!(columns1 %in% testcols1)]
    traincols2 <- columns2[!(columns2 %in% testcols2)]
    traincols <- c(traincols1,traincols2)
    traincols1.length <- length(traincols1)
    traincols2.length <- length(traincols2)
    traincols.length <- length(traincols)
    trainsample <- elist$E[,traincols]
    trainsample1 <- elist$E[,traincols1]
    trainsample2 <- elist$E[,traincols2]
    
    idx <- paste(elist$genes$Block,elist$genes$Row,elist$genes$Column)
    p <- rep(0, row.number)

    trainsample <- cbind(p, trainsample)
    colnames(trainsample) <- c("Dummy Score",traincols)
    rownames(trainsample) <- idx
    trainsample <- trainsample[order(trainsample[,1]),]

    testsample <- cbind(p, testsample)
    colnames(testsample) <- c("Dummy Score",testcols)
    rownames(testsample) <- idx
    testsample <- testsample[order(testsample[,1]),]

    return(list(trainsample=trainsample, testsample=testsample,
      traincols1=traincols1, traincols2=traincols2, testcols=testcols))
}
#+++++++++++++++++++++++++++ noFS +++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++ tTestFS +++++++++++++++++++++++++++++++++++++++
tTestFS <- function(elist=NULL, columns1=NULL, columns2=NULL, label1="A",
  label2="B"){
    
    if(is.null(elist) || is.null(columns1) || is.null(columns2)) {
        stop("ERROR: Not all mandatory arguments have been defined!")
    }
    
    row.number <- nrow(elist$E)
    n1 <- length(columns1)
    n2 <- length(columns2)

    testcols1 <- columns1[columns1 %in% sample(columns1,n1/3,replace=FALSE)]
    testcols2 <- columns2[columns2 %in% sample(columns2,n2/3,replace=FALSE)]
    testcols <- c(testcols1,testcols2)
    testcols1.length <- length(testcols1)
    testcols2.length <- length(testcols2)
    testcols.length <- length(testcols)
    testsample <- elist$E[,testcols]

    traincols1 <- columns1[!(columns1 %in% testcols1)]
    traincols2 <- columns2[!(columns2 %in% testcols2)]
    traincols <- c(traincols1,traincols2)
    traincols1.length <- length(traincols1)
    traincols2.length <- length(traincols2)
    traincols.length <- length(traincols)
    trainsample <- elist$E[,traincols]
    trainsample1 <- elist$E[,traincols1]
    trainsample2 <- elist$E[,traincols2]

    idx <- paste(elist$genes$Block,elist$genes$Row,elist$genes$Column)
    p <- vector(mode="numeric", length=row.number)
    for(zeile in 1:row.number) {
        p.tmp <- try(t.test(x=trainsample1[zeile,] ,
          y=trainsample2[zeile,])$p.value, TRUE)
        if(!is.numeric(p.tmp)){p.tmp <- 1}
        p[zeile] <- p.tmp
    }

    trainsample <- cbind(p, trainsample)
    colnames(trainsample) <- c("p value",traincols)
    rownames(trainsample) <- idx
    trainsample <- trainsample[order(trainsample[,1]),]
    
    testsample <- cbind(p, testsample)
    colnames(testsample) <- c("p value",testcols)
    rownames(testsample) <- idx
    testsample <- testsample[order(testsample[,1]),]

    return(list(trainsample=trainsample, testsample=testsample,
      traincols1=traincols1, traincols2=traincols2, testcols=testcols))
}
#+++++++++++++++++++++++++++ tTestFS +++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++ tTest +++++++++++++++++++++++++++++++++++++++
tTest <- function(elist=NULL, columns1=NULL, columns2=NULL, label1="A",
  label2="B", discard.threshold=0.5, fold.thresh=1.5, discard.features=FALSE){
    
    if(is.null(elist) || is.null(columns1) || is.null(columns2)) {
        stop("ERROR: Not all mandatory arguments have been defined!")
    }
    row.number <- nrow(elist$E)
    col.number <- ncol(elist$E)
    n1 <- length(columns1)
    n2 <- length(columns2)
    discard <- c()
    count <- 0

    output.matrix <- matrix(nrow=row.number,ncol=8+n1+n2)
    for(zeile in 1:row.number) {
        p <- try(t.test(x=elist$E[zeile,1:n1] ,
          y=elist$E[zeile,(n1+1):col.number])$p.value, TRUE)
        if(!is.numeric(p)){p <- 1}

        mean1 <- mean(elist$E[zeile,1:n1])
        mean2 <- mean(elist$E[zeile,(n1+1):col.number])
        fold <- 2^(mean1 - mean2) # --> data must be in log2 scale!
        #if(round(p,1) < round(discard.threshold,1)) {   
        #---> ISSUE: there  may be no exact representation for floats!!!
        if(p < discard.threshold && 
          (fold > fold.thresh || fold < 1/fold.thresh)) {
            
            idx.tmp <- paste(elist$genes$Block[zeile],elist$genes$Row[zeile],
              elist$genes$Column[zeile])
            count <- count+1
            output.matrix[count,] <- c(idx.tmp, elist$genes$Description[zeile],
              elist$genes$Name[zeile], elist$genes$ID[zeile], p, mean1, mean2,
              (mean1 / mean2), elist$E[zeile,])
        }else if(discard.features) {                                                             
            discard <- c(discard,zeile)
        }
        #message('tTest() - computing p values: ',
        #  round(zeile*100/{row.number-1}), '%\r', sep="")
    }
    #message('tTest() - computing p values: ', '100', '%\n', sep="")

    output.matrix <- data.frame(output.matrix[1:count,])
    colnames(output.matrix) <- c("BRC", "Description", "Name", "ID", "p Value",
      paste(label1, " Mean", sep=""), paste(label2, " Mean", sep=""),
      "Fold Change", colnames(elist$E))

    if(discard.features){
        return(list(results=output.matrix, discard=discard))
    }else{
        return(output.matrix)
    }
}
#+++++++++++++++++++++++++++ tTest +++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++ mrmr +++++++++++++++++++++++++++++++++++++++
mrmr <- function(elist=NULL, columns1=NULL, columns2=NULL, label1="A",
  label2="B", discard.threshold=300, discard.features=FALSE){
    
    
    if(is.null(elist) || is.null(columns1) || is.null(columns2)) {
        stop("ERROR: Not all mandatory arguments have been defined!")
    }
    row.number <- nrow(elist$E)
    col.number <- ncol(elist$E)
    n1 <- length(columns1)
    n2 <- length(columns2)
    discard <- 1:row.number

    output.matrix <- matrix(nrow=row.number,ncol=8+n1+n2)
    mrmr.dat <-
      mRMR.data(data=as.data.frame(t(rbind(c(rep(1,n1),rep(2,n2)),elist$E))))
    mrmr.results <-
      mRMR.classic(data = mrmr.dat, target_indices = c(1),
      feature_count = discard.threshold)
    mrmr.idx <- solutions(mrmr.results)$`1`[,1]-1
    mrmr.scores <- mrmr.results@scores$`1`[,1]
    discard <- discard[!(discard %in% mrmr.idx)]
    
    elist$E <- elist$E[mrmr.idx,]
    elist$genes <- elist$genes[mrmr.idx,]
    idx.tmp <-
      paste(elist$genes$Block,elist$genes$Row,elist$genes$Column) 
    mean1 <- apply(elist$E[,columns1], 1, mean)
    mean2 <- apply(elist$E[,columns2], 1, mean)
    fold <- 2^(mean1 - mean2) # --> data must be in log2 scale!
    output.matrix <- cbind(idx.tmp, elist$genes$Description, elist$genes$Name,
      elist$genes$ID, mrmr.scores, mean1, mean2, fold, elist$E)

    output.matrix <- data.frame(output.matrix)
    colnames(output.matrix) <- c("BRC", "Description", "Name", "ID", "Score",
      paste(label1, " Mean", sep=""), paste(label2, " Mean", sep=""),
      "Fold Change", colnames(elist$E))

    if(discard.features){
        return(list(results=output.matrix, discard=discard))
    }else{
        return(output.matrix)
    }
}
#+++++++++++++++++++++++++++ mrmr +++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++ mrmrFS +++++++++++++++++++++++++++++++++++++++
mrmrFS <- function(elist=NULL, columns1=NULL, columns2=NULL, label1="A",
  label2="B", candidate.number=300){
    
    
    if(is.null(elist) || is.null(columns1) || is.null(columns2)) {
        stop("ERROR: Not all mandatory arguments have been defined!")
    }
    
    n1 <- length(columns1)
    n2 <- length(columns2)

    testcols1 <- columns1[columns1 %in% sample(columns1,n1/3,replace=FALSE)]
    testcols2 <- columns2[columns2 %in% sample(columns2,n2/3,replace=FALSE)]
    testcols <- c(testcols1,testcols2)
    testcols1.length <- length(testcols1)
    testcols2.length <- length(testcols2)
    testcols.length <- length(testcols)
    testsample <- elist$E[,testcols]

    traincols1 <- columns1[!(columns1 %in% testcols1)]
    traincols2 <- columns2[!(columns2 %in% testcols2)]
    traincols <- c(traincols1,traincols2)
    traincols1.length <- length(traincols1)
    traincols2.length <- length(traincols2)
    traincols.length <- length(traincols)
    trainsample <- elist$E[,traincols]
    trainsample1 <- elist$E[,traincols1]
    trainsample2 <- elist$E[,traincols2]
    
    mrmr.dat <- mRMR.data(data=as.data.frame(t(rbind(c(rep(1,n1),rep(2,n2)),
      trainsample))))
    mrmr.results <- mRMR.classic(data = mrmr.dat, target_indices = c(1),
      feature_count = candidate.number)
    mrmr.idx <- solutions(mrmr.results)$`1`[,1]-1
    mrmr.scores <- mrmr.results@scores$`1`[,1]      
    
    trainsample <- trainsample[mrmr.idx,]
    trainsample1 <- trainsample1[mrmr.idx,]
    trainsample2 <- trainsample2[mrmr.idx,]
    testsample <- testsample[mrmr.idx,]
    genes <- elist$genes[mrmr.idx,]
    idx <- paste(genes$Block,genes$Row,genes$Column)
    
    trainsample <- cbind(mrmr.scores, trainsample)
    colnames(trainsample) <- c("mrmr score",traincols)
    rownames(trainsample) <- idx
    trainsample <- trainsample[order(trainsample[,1]),]
    
    testsample <- cbind(mrmr.scores, testsample)
    colnames(testsample) <- c("mrmr score",testcols)
    rownames(testsample) <- idx
    testsample <- testsample[order(testsample[,1]),]

    return(list(trainsample=trainsample, testsample=testsample,
      traincols1=traincols1, traincols2=traincols2, testcols=testcols))
}
#+++++++++++++++++++++++++++ mrmrFS +++++++++++++++++++++++++++++++++++++++

##+++++++++++++++++++++++++++ mMsFS +++++++++++++++++++++++++++++++++++++++
#mMsFS <- function(elist=NULL, columns1=NULL, columns2=NULL, label1="A",
#  label2="B", above=1500, between=400, mMs.matrix1=NULL, mMs.matrix2=NULL){
#    
#    if(is.null(elist) || is.null(columns1) || is.null(columns2) ||
#      is.null(mMs.matrix1) || is.null(mMs.matrix2)) {
#        stop("ERROR: Not all mandatory arguments have been defined!")
#    }
#    
#    row.number <- nrow(elist$E)
#    n1 <- length(columns1)
#    n2 <- length(columns2)
#
#    testcols1 <- columns1[columns1 %in% sample(columns1,n1/3,replace=FALSE)]
#    testcols2 <- columns2[columns2 %in% sample(columns2,n2/3,replace=FALSE)]
#    testcols <- c(testcols1,testcols2)
#    testcols1.length <- length(testcols1)
#    testcols2.length <- length(testcols2)
#    testcols.length <- length(testcols)
#    testsample <- elist$E[,testcols]
#
#    traincols1 <- columns1[!(columns1 %in% testcols1)]
#    traincols2 <- columns2[!(columns2 %in% testcols2)]
#    traincols <- c(traincols1,traincols2)
#    traincols1.length <- length(traincols1)
#    traincols2.length <- length(traincols2)
#    traincols.length <- length(traincols)
#    trainsample <- elist$E[,traincols]
#    trainsample1 <- elist$E[,traincols1]
#    trainsample2 <- elist$E[,traincols2]
#
#    idx <- paste(elist$genes$Block,elist$genes$Row,elist$genes$Column)
#    p <- vector(mode="numeric", length=row.number)
#    x <- c(1:traincols2.length)
#    y <- c(1:traincols1.length)
#    for(zeile in 1:row.number) {
#        mCount.results <- 1 + matrix(unlist(mapply("mCount", idx=x,
#          MoreArgs=list(vector1=trainsample1[zeile,],
#          vector2=trainsample2[zeile,], above=above, between=between))),2,
#          traincols2.length)
#        indexes <- array(c(x, mCount.results[1,]), dim=c(traincols2.length,2))
#        pij.results <- mMs.matrix1[indexes]
#        min.order1 <- which.min(pij.results)
#        min.p1 <- pij.results[min.order1]
#        #min.count1 <- mCount.results[1,min.order1]
#        #min.cutoff1 <- mCount.results[2,min.order1]
#
#        mCount.results <- 1 + matrix(unlist(mapply("mCount", idx=y,
#          MoreArgs=list(vector1=trainsample2[zeile,],
#          vector2=trainsample1[zeile,], above=above, between=between))),2,
#          traincols1.length)
#        indexes <- array(c(y, mCount.results[1,]), dim=c(traincols1.length,2))
#        pij.results <- mMs.matrix2[indexes]
#        min.order2 <- which.min(pij.results)
#        min.p2 <- pij.results[min.order2]
#        #min.count2 <- mCount.results[1,min.order2]
#        #min.cutoff2 <- mCount.results[2,min.order2]
#        p[zeile] <- min(min.p1,min.p2)
#    }
#
#    trainsample <- cbind(p, trainsample)
#    colnames(trainsample) <- c("p value",traincols)
#    rownames(trainsample) <- idx
#    trainsample <- trainsample[order(trainsample[,1]),]
#    
#    testsample <- cbind(p, testsample)
#    colnames(testsample) <- c("p value",testcols)
#    rownames(testsample) <- idx
#    testsample <- testsample[order(testsample[,1]),]
#
#    return(list(trainsample=trainsample, testsample=testsample,
#      traincols1=traincols1, traincols2=traincols2, testcols=testcols))
#}
##+++++++++++++++++++++++++++ mMsFS +++++++++++++++++++++++++++++++++++++++
#
##+++++++++++++++++++++++++++ mMs +++++++++++++++++++++++++++++++++++++++
##column 01: BRC
##column 02: Description
##column 03: Name
##column 04: ID
##column 05: Group1 Order
##column 06: Group1 Count
##column 07: Group1 Cutoff
##column 08: Group2 Order
##column 09: Group2 Count
##column 10: Group2 Cutoff
##column 11: minimum M statistic
##column 12: Candidate
##column 13: Group1 Mean
##column 14: Group2 Mean
##column 15: Fold Change
##column 16: Group1 Count Prospector
##column 17: Group2 Count Prospector
##column 18: Cutoff Prospector
##columns 19 - (18+n1+n2): colnames(elist$E)
#mMs <- function(elist=NULL, columns1=NULL, columns2=NULL, label1="A",
#  label2="B", above=1500, between=400, discard.threshold=0.5, fold.thresh=1.5,
#  discard.features=FALSE, mMs.matrix1=NULL, mMs.matrix2=NULL){
#    
#    if(is.null(elist) || is.null(columns1) || is.null(columns2) ||
#      is.null(mMs.matrix1) || is.null(mMs.matrix2)) {
#        stop("ERROR: Not all mandatory arguments have been defined!")
#    }
#    
#    row.number <- nrow(elist$E)
#    col.number <- ncol(elist$E)
#    n1 <- length(columns1)
#    n2 <- length(columns2)
#    discard <- c()
#    count <- 0
#
#    output.matrix <- matrix(nrow=row.number,ncol=18+n1+n2)
#    x <- 1:n2
#    y <- 1:n1
#    for(zeile in 1:row.number) {
#        mCount.results <- 1 + matrix(unlist(mapply("mCount", idx=x,
#          MoreArgs=list(vector1=elist$E[zeile,1:n1],
#          vector2=elist$E[zeile,(n1+1):col.number], above=above,
#          between=between))),2,n2)
#        indexes <- array(c(x, mCount.results[1,]), dim=c(n2,2))
#        pij.results <- mMs.matrix1[indexes]
#        min.order1 <- which.min(pij.results)
#        min.p1 <- pij.results[min.order1]
#        min.count1 <- mCount.results[1,min.order1]
#        min.cutoff1 <- mCount.results[2,min.order1]
#
#        mCount.results <- 1 + matrix(unlist(mapply("mCount", idx=y,
#          MoreArgs=list(vector1=elist$E[zeile,(n1+1):col.number],
#          vector2=elist$E[zeile,1:n1], above=above, between=between))),2,n1)
#        indexes <- array(c(y, mCount.results[1,]), dim=c(n1,2))
#        pij.results <- mMs.matrix2[indexes]
#        min.order2 <- which.min(pij.results)
#        min.p2 <- pij.results[min.order2]
#        min.count2 <- mCount.results[1,min.order2]
#        min.cutoff2 <- mCount.results[2,min.order2]
#        
#        if(min.p1 < min.p2) {
#            prosp.count1 <- min.count1 - 1
#            prosp.count2 <- min.order1 - 1 
#            prosp.cutoff <- min.cutoff1 - 1 
#        }else if(min.p1 > min.p2) {
#            prosp.count2 <- min.count2 - 1
#            prosp.count1 <- min.order2 - 1
#            prosp.cutoff <- min.cutoff2 - 1
#        }else{
#            prosp.count1 <- min.count1 - 1
#            prosp.count2 <- min.count2 - 1
#            prosp.cutoff <- min.cutoff2 - 1   
#            #Prospector seems to use group 2 cutoff when winner group is
#            #indeterminate (checked for some examples)
#        }
#
#        mean1 <- mean(elist$E[zeile,1:n1])
#        mean2 <- mean(elist$E[zeile,(n1+1):col.number])
#        fold <- mean1 / mean2
#        if(min(min.p1,min.p2) < discard.threshold && (fold > fold.thresh ||
#          fold < 1/fold.thresh)) {
#            
#            idx.tmp <- paste(elist$genes$Block[zeile],elist$genes$Row[zeile],
#              elist$genes$Column[zeile])
#            count <- count+1
#            
#            output.matrix[count,] <- c(idx.tmp, elist$genes$Description[zeile],
#              elist$genes$Name[zeile], elist$genes$ID[zeile], min.order1,
#              min.count1, min.cutoff1, min.order2, min.count2, min.cutoff2,
#              min(min.p1,min.p2),
#              if(min.p2 < min.p1) label2 else if (min.p2 > min.p1) label1
#              else "indeterminate", mean1, mean2, fold, prosp.count1,
#              prosp.count2, prosp.cutoff, elist$E[zeile,])
#            
#            #output.matrix[count,] <- c(idx.tmp, elist$genes[zeile,4],
#            # elist$genes[zeile,5], elist$genes[zeile,6], min.order1,
#            # min.count1, min.cutoff1, min.order2, min.count2, min.cutoff2,
#            # round(min(min.p1,min.p2),9),
#            # if(min.p2 < min.p1) label2 else if (min.p2 > min.p1) label1 else
#            # "indeterminate", mean1, mean2, fold, prosp.count1, prosp.count2,
#            # prosp.cutoff, elist$E[zeile,])
#        }else if(discard.features) {
#            discard <- c(discard,zeile)
#        }
#        #message('computing minimum M statistics: ',
#        # round(zeile*100/{row.number-1}), '%\r', sep="")
#    }
#    #message('computing minimum M statistics: ', '100', '%\n', sep="")
#
#    output.matrix <- data.frame(output.matrix[1:count,])
#    
#    colnames(output.matrix) <- c("BRC","Description","Name","ID",
#     paste(label1," Order",sep=""),paste(label1," Count",sep=""),
#     paste(label1," Cutoff",sep=""),paste(label2," Order",sep=""),
#     paste(label2," Count",sep=""),paste(label2," Cutoff",sep=""),
#     "minimum M statistic","Candidate",paste(label1," Mean",sep=""),
#     paste(label2," Mean",sep=""),"Fold Change",
#     paste(label1,".Count Prosp.",sep=""),
#     paste(label2,".Count Prosp.",sep=""),"Cutoff Prosp.",colnames(elist$E))
#
#    if(discard.features){
#        return(list(results=output.matrix, discard=discard))
#    }else{
#        return(output.matrix)
#    }
#}
##+++++++++++++++++++++++++++ mMs +++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++ mMsFS2 +++++++++++++++++++++++++++++++++++++++
mMsFS2 <- function(elist=NULL, columns1=NULL, columns2=NULL, label1="A",
  label2="B", above=1500, between=400, mMs.matrix1=NULL, mMs.matrix2=NULL){
    
    if(is.null(elist) || is.null(columns1) || is.null(columns2) ||
      is.null(mMs.matrix1) || is.null(mMs.matrix2)) {
        stop("ERROR: Not all mandatory arguments have been defined!")
    }
    
    row.number <- nrow(elist$E)
    n1 <- length(columns1)
    n2 <- length(columns2)
    
    tc1 <- sort(sample(0:{n1-1},n1/3,replace=FALSE))
    tc2 <- sort(sample({n1}:{n1+n2-1},n2/3,replace=FALSE))
    sampling.results <- sampling(x=elist$E, c=c(columns1,columns2), c1=columns1,
      c2=columns2, tc1=tc1, tc2=tc2)
    
    testcols1 <- sampling.results$testcols1 
    testcols2 <- sampling.results$testcols2 
    testcols <- sampling.results$testcols
    testsample <- sampling.results$testmatrix
    
    traincols1 <- sampling.results$traincols1 
    traincols2 <- sampling.results$traincols2 
    traincols <- sampling.results$traincols
    trainsample <- sampling.results$trainmatrix
    trainsample1 <- sampling.results$trainmatrix1
    trainsample2 <- sampling.results$trainmatrix2 

    mCount.results1 <-
      mCount(x=trainsample1, y=trainsample2, z=mMs.matrix1, a=above, b=between)
    mCount.results2 <-
      mCount(x=trainsample2, y=trainsample1, z=mMs.matrix2, a=above, b=between)
    joined.mCount.results <-
      joinMCountResults(p1=mCount.results1$mMs, p2=mCount.results2$mMs,
      m1=mCount.results2$means2, m2=mCount.results1$means2, l1=label1,
      l2=label2)
    idx <- paste(elist$genes$Block,elist$genes$Row,elist$genes$Column)
    
    trainsample <- cbind(joined.mCount.results$p, trainsample)
    colnames(trainsample) <- c("M Statistic", traincols)
    rownames(trainsample) <- idx
    trainsample <- trainsample[order(trainsample[,1]),]
    
    testsample <- cbind(joined.mCount.results$p, testsample)
    colnames(testsample) <- c("M Statistic", testcols)
    rownames(testsample) <- idx
    testsample <- testsample[order(testsample[,1]),]
    
    return(list(trainsample=trainsample, testsample=testsample,
      traincols1=traincols1, traincols2=traincols2, testcols=testcols))
}
#+++++++++++++++++++++++++++ mMsFS2 +++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++ mMs2 +++++++++++++++++++++++++++++++++++++++
#column 01: BRC
#column 02: Description
#column 03: Name
#column 04: ID
#column 05: Group1 Order
#column 06: Group1 Count               
#column 07: Group1 Cutoff
#column 08: Group2 Order
#column 09: Group2 Count
#column 10: Group2 Cutoff
#column 11: minimum M statistic
#column 12: Candidate
#column 13: Group1 Mean
#column 14: Group2 Mean
#column 15: Fold Change
#column 16: Group1 Count Prospector
#column 17: Group2 Count Prospector
#column 18: Cutoff Prospector
#columns 19 - (18+n1+n2): colnames(elist$E)
mMs2 <- function(elist=NULL, columns1=NULL, columns2=NULL, label1="A",
  label2="B", above=1500, between=400, discard.threshold=0.5, fold.thresh=1.5,
  discard.features=FALSE, mMs.matrix1=NULL, mMs.matrix2=NULL){
    
    if(is.null(elist) || is.null(columns1) || is.null(columns2) ||
      is.null(mMs.matrix1) || is.null(mMs.matrix2)) {
        stop("ERROR: Not all mandatory arguments have been defined!")
    }
    
    row.number <- nrow(elist$E)
    col.number <- ncol(elist$E)
    n1 <- length(columns1)
    n2 <- length(columns2)
    discard <- c()
    count <- 0
    message(paste0("row.numer: ", row.number))
    message(paste0("col.numer: ", col.number))

    output.matrix <- matrix(nrow=row.number,ncol=18+n1+n2)
    mCount.results1 <- mCount(x=elist$E[,1:n1], y=elist$E[,{n1+1}:col.number],
      z=mMs.matrix1, a=above, b=between)
    mCount.results2 <- mCount(x=elist$E[,{n1+1}:col.number], y=elist$E[,1:n1],
      z=mMs.matrix2, a=above, b=between)
    
    for(zeile in 1:row.number) {
        m.order1 <- mCount.results1$order[zeile]
        m.count1 <- mCount.results1$count[zeile]
        m.cutoff1 <- mCount.results1$cutoff[zeile]
        mean2 <- mCount.results1$means2[zeile]
        #mean1 <- mean(elist$E[zeile,1:n1])
        
        m.order2 <- mCount.results2$order[zeile]
        m.count2 <- mCount.results2$count[zeile]
        m.cutoff2 <- mCount.results2$cutoff[zeile]
        mean1 <- mCount.results2$means2[zeile]
        #mean2 <- mean(elist$E[zeile,(n1+1):col.number])
        
        fold <- mean1 / mean2
        if(mCount.results1$mMs[zeile] < mCount.results2$mMs[zeile]){
            p <- mCount.results1$mMs[zeile] 
            label <- label1
            prosp.count1 <- m.count1 - 1
            prosp.count2 <- m.order1 - 1
            prosp.cutoff <- m.cutoff1 
        }else if(mCount.results1$mMs[zeile] > mCount.results2$mMs[zeile]){
            p <- mCount.results2$mMs[zeile]
            label <- label2
            prosp.count2 <- m.count2 - 1
            prosp.count1 <- m.order2 - 1
            prosp.cutoff <- m.cutoff2
        }else{
            p <- mCount.results1$mMs[zeile]
            label <- "indeterminate"
            prosp.count1 <- m.count1 - 1
            prosp.count2 <- m.count2 - 1
            prosp.cutoff <- m.cutoff2   
            #Prospector seems to use group 2 cutoff when winner group is
            #indeterminate (checked for some examples)
        }
        
        if(p < discard.threshold && (fold > fold.thresh ||
          fold < 1/fold.thresh)) {
            
            idx.tmp <- paste(elist$genes$Block[zeile],elist$genes$Row[zeile],
              elist$genes$Column[zeile])
            count <- count+1
            output.matrix[count,] <- c(idx.tmp, elist$genes$Description[zeile],
              elist$genes$Name[zeile], elist$genes$ID[zeile], m.order1,
              m.count1,m.cutoff1, m.order2, m.count2, m.cutoff2, round(p,9),
              label, mean1, mean2, fold, prosp.count1, prosp.count2,
              prosp.cutoff, elist$E[zeile,])
              
        }else if(discard.features) {
            discard <- c(discard,zeile)
        }
    }
    
    if(count==0){
        warning("No features fulfill the prefiltering thresholds!")
    }

    output.matrix <- data.frame(output.matrix[1:count,])
    colnames(output.matrix) <- c("BRC","Description","Name","ID",
      paste(label1," Order",sep=""),paste(label1," Count",sep=""),
      paste(label1," Cutoff",sep=""),paste(label2," Order",sep=""),
      paste(label2," Count",sep=""),paste(label2," Cutoff",sep=""),
      "minimum M statistic","Candidate",paste(label1," Mean",sep=""),
      paste(label2," Mean",sep=""),"Fold Change",
      paste(label1,".Count Prosp.",sep=""),paste(label2,".Count Prosp.",sep=""),
      "Cutoff Prosp.",colnames(elist$E))

    if(discard.features){
        return(list(results=output.matrix, discard=discard))
    }else{
        return(output.matrix)
    }
}
#+++++++++++++++++++++++++++ mMs2 +++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++ loadGPR +++++++++++++++++++++++++++++++++++++++
loadGPR <- function(gpr.path=NULL, targets.path=NULL, array.type=NULL,
  aggregation="none", array.columns=list(E="F635 Median", Eb="B635 Median"),
  array.annotation=c("Block", "Column", "Row", "Description", "Name", "ID"),
  description=NULL, description.features=NULL, description.discard=NULL){

    if(is.null(gpr.path) || is.null(targets.path) || is.null(array.type)) {
        stop("ERROR: Not all mandatory arguments have been defined!")
    }
    if(!is.null(description) && is.null(description.features)) {
        stop(paste0("ERROR: If 'description' is defined 'description.features'",
        " is mandatory!"))
    }
    if(!is.null(description) && is.null(description.discard)) {
        stop(paste0("ERROR: If 'description' is defined 'description.discard' ",
        "is mandatory!"))
    }
    targets <- readTargets(targets.path)
    elist.tmp <- read.maimages(files=targets, path=gpr.path,
      source="genepix.median", columns=array.columns,
      annotation=array.annotation)
    if(length(setdiff(array.annotation, colnames(elist.tmp$genes))) > 0){
        warning(paste0("The following columns were not imported: ",
        paste(setdiff(array.annotation, colnames(elist.tmp$genes)),
        collapse=", ")),
        ". If one of these columns is mandatory (see the loadGPR() ",
        "documentation to check whether a column is mandatory) this may cause ",
        "serious errors.")
    }
     
    # For ProtoArrays: construct missing column 'Description' from the column
    # 'Name'. Alternatively, another column in elist$genes and other
    # regular expressions can be specified.  
    if(is.null(elist.tmp$genes$Description) && array.type=="ProtoArray"){
        elist.tmp$genes <- cbind(elist.tmp$genes,
          rep("NA",nrow(elist.tmp$genes)), stringsAsFactors=FALSE)
        colnames(elist.tmp$genes)[ncol(elist.tmp$genes)] <- "Description"
        features <- grep("^Hs~", elist.tmp$genes$Name)
        controls <- grep("^Hs~", elist.tmp$genes$Name, invert=TRUE)
        empty <- grep("(^Empty|^EMPTY)", elist.tmp$genes$Name) 
        #controls <-
        #  grep("(^HumanIg|^Anti-HumanIg|^V5|Anti-human-Ig|V5Control)",
        #  elist.tmp$genes$Name)
        elist.tmp$genes$Description[features] <- "Feature"
        elist.tmp$genes$Description[controls] <- "Control"
        elist.tmp$genes$Description[empty] <- "Empty"  
    }else if(is.null(elist.tmp$genes$Description) && !is.null(description)
      && !is.null(description.features) && !is.null(description.discard)){
      
        elist.tmp$genes <- cbind(elist.tmp$genes,
          rep("NA",nrow(elist.tmp$genes)), stringsAsFactors=FALSE)
        colnames(elist.tmp$genes)[ncol(elist.tmp$genes)] <- "Description"
        
        features <- grep(description.features, elist.tmp$genes[,description])
        controls <- grep(description.features, elist.tmp$genes[,description],
          invert=TRUE)
        empty <- grep(description.discard, elist.tmp$genes[,description])
        
        elist.tmp$genes$Description[features] <- "Feature"
        elist.tmp$genes$Description[controls] <- "Control"
        elist.tmp$genes$Description[empty] <- "Empty"    
    } 

    # If no match, grep will return 'integer(0)' resulting in an empty elist$E
    # in the following lines. Hence the following if-statements are necessary:
    if(any(grep("(^Empty$|^EMPTY$|^NA$)", elist.tmp$genes$Description))){
      elist.tmp <- elist.tmp[-grep("(^Empty$|^EMPTY$|^NA$)",
       elist.tmp$genes$Description),]
    }
    if(any(grep("(^Control$|^CONTROL$)", elist.tmp$genes$Description))){
      elist <- elist.tmp[-grep("(^Control$|^CONTROL$)", elist.tmp$genes$Description),]
    }

    colnames(elist$E) <- targets$ArrayID
    colnames(elist$Eb) <- colnames(elist$E)

    if(aggregation=="none"){
        cat("No aggregation performed.\n")
    }else if(array.type=="ProtoArray" && aggregation=="min"){
        row.len <- nrow(elist$E)
        col.len <- ncol(elist$E)
        tmp.col.len <- (row.len*col.len)/2
        elist$E[row(elist$E)[,1]%%2==1] <-
          matrix(apply(matrix(elist$E,2,tmp.col.len),2,min),row.len/2,col.len)
        elist <- elist[-row(elist)[,1]%%2==1,]
    }else if(array.type=="ProtoArray" && aggregation=="mean"){
        elist$E[row(elist$E)[,1]%%2==1,] <-
          (elist$E[row(elist$E)[,1]%%2==1,]+elist$E[row(elist$E)[,1]%%2==0,])/2
        elist <- elist[-row(elist)[,1]%%2==1,]
    }else{
        stop(paste0("ERROR: This aggregation approach is for the specified ",
        "array type not supported.\nPlease contact the PAA maintainer for a ",
        "possible implementation of the requested agrregation approach for ",
        "this array type."))
    }
    
    # Adding controls data to the EListRaw object for ProtoArrays in order to
    # use it for rlm normalization.
    if(array.type == "ProtoArray"){
        elist.tmp <-
          elist.tmp[grep("^Control$", elist.tmp$genes$Description),]
        elist.tmp$genes$Name <-
          gsub('^([0-9A-Za-z_-]*)~(.*)', '\\1', elist.tmp$genes$Name)
        elist$C <- elist.tmp$E
        elist$Cb <- elist.tmp$Eb
        elist$cgenes <- elist.tmp$genes
        colnames(elist$C) <- colnames(elist$E)
        colnames(elist$Cb) <- colnames(elist$E)
    }else if(array.type != "ProtoArray" && !is.null(description)){
        elist.tmp <-
          elist.tmp[grep("^Control$", elist.tmp$genes$Description),]
        elist.tmp$genes$Name <-
          gsub('^([0-9A-Za-z_-]*)~(.*)', '\\1', elist.tmp$genes$Name)
        elist$C <- elist.tmp$E
        elist$Cb <- elist.tmp$Eb
        elist$cgenes <- elist.tmp$genes
        colnames(elist$C) <- colnames(elist$E)
        colnames(elist$Cb) <- colnames(elist$E)
    }
    
    # Adding this custom component in order to check the microarray type in
    # type-specific functions
    elist$array.type <- array.type
    
    return(elist)
}
#+++++++++++++++++++++++++++ loadGPR +++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++ backgroundCorrect +++++++++++++++++++++++++++++++++
#-->limma: backgroundCorrect(RG, method="auto", offset=0, printer=RG$printer,
#  normexp.method="saddle", verbose=TRUE)
#+++++++++++++++++++++++++++ backgroundCorrect +++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++ batchAdjust +++++++++++++++++++++++++++++++++++++++
batchAdjust <- function(elist=NULL, log=NULL){
    if(is.null(elist) || is.null(log)) {
        stop("ERROR: Not all mandatory arguments have been defined!")
    }
    if(log){
        elist$E <-
          ComBat(dat=elist$E, mod=factor(elist$targets$Group),
          batch=factor(elist$targets$Batch))
    }else if(!log){
        elist$E <-
          2^ComBat(dat=log2(elist$E), mod=factor(elist$targets$Group),
          batch=factor(elist$targets$Batch))
    }else{
        stop("ERROR in batchAdjust: log not 'TRUE' or 'FALSE'.")
    }
    return(elist)
}
#+++++++++++++++++++++++++++ batchAdjust +++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++ batchFilter +++++++++++++++++++++++++++++++++++++++
batchFilter <- function(elist=NULL, lot1=NULL, lot2=NULL, log=NULL,
  p.thresh=0.05, fold.thresh=1.5, output.path=NULL){
    
    if(is.null(elist) || is.null(lot1) || is.null(lot2) || is.null(log)) {
        stop("batchFilter - ERROR: Not all mandatory arguments defined!")
    }
    
    if(log==FALSE){
        E <- log2(elist$E)
    }else{
        E <- elist$E
    }
    
    row.number <- nrow(E)
    col.number <- ncol(E)
    discard <- c()
    output <- matrix(nrow=row.number, ncol=4)
    colnames(output) <- c("BRC", "P-value", "Fold change", "Description")

    for(zeile in 1:row.number) {
        p <- try(t.test(x=E[zeile,lot1],
          y=E[zeile,lot2])$p.value, FALSE)
        if(!is.numeric(p)){p <- 1}
        
        fold <- 2^(mean(E[zeile, lot1]) - mean(E[zeile, lot2]))
                # --> data must be in log2 scale!
        
        output[zeile,] <- c(paste(elist$genes$Block[zeile],
          elist$genes$Row[zeile], elist$genes$Column[zeile]), p, fold,
          elist$genes$Description[zeile])
        if(p < p.thresh && (fold > fold.thresh || fold < 1/fold.thresh)) {
            #---> ISSUE: there  may be no exact representation for floats!!!
            discard <- c(discard,zeile)
        }
    }
    message(paste("batchFilter - number of features to discard: ",
      length(discard), sep=""), "\n")
    if(is.null(output.path)) {
        plot(x=log2(as.numeric(output[,3])), y=-log10(as.numeric(output[,2])),
          xlab="log2(fold change)", ylab="-log10(p-value)",
          main="Batch Filter Volcano")
        if(length(discard)>0){
            points(x=log2(as.numeric(output[discard,3])),
              y=-log10(as.numeric(output[discard,2])), col="red")
        }
        abline(h=-log10(p.thresh), lty=2, col="red")
        abline(v=log2(fold.thresh), lty=2, col="red")
        abline(v=log2(1/fold.thresh), lty=2, col="red")
        if(length(discard)>0){
            elist <- elist[-discard,]
        }
    }else{
        write.table(x=output, file=paste(output.path,
          "/batch_filter.txt",sep=""), sep="\t", eol="\n", row.names=FALSE,
          quote=FALSE)
        write.table(x=output[discard,], file=paste(output.path,
          "/batch_filter_discarded.txt",sep=""), sep="\t", eol="\n",
          row.names=FALSE, quote=FALSE)
        tiff(paste(output.path,"/batch_filter.tiff",sep=""), width=2000,
          height=2000, pointsize=15, compression="lzw", res=300)
            plot(x=log2(as.numeric(output[,3])),
              y=-log10(as.numeric(output[,2])), xlab="log2(fold change)",
              ylab="-log10(p-value)", main="Batch Filter Volcano")
            if(length(discard)>0){
                points(x=log2(as.numeric(output[discard,3])),
                  y=-log10(as.numeric(output[discard,2])), col="red")
            }
            abline(h=-log10(p.thresh), lty=2, col="red")
            abline(v=log2(fold.thresh), lty=2, col="red")
            abline(v=log2(1/fold.thresh), lty=2, col="red")
        dev.off()
        if(length(discard)>0){
            tiff(paste(output.path,"/batch_filter_discarded.tiff",sep=""),
              width=2000, height=2000, pointsize=15, compression="lzw", res=300)
                plot(x=log2(as.numeric(output[-discard,3])),
                  y=-log10(as.numeric(output[-discard,2])),
                  xlab="log2(fold change)", ylab="-log10(p-value)",
                  main="Batch Filter Volcano")
            dev.off()
            elist <- elist[-discard,]
        }
    }
    return(elist)
}
#+++++++++++++++++++++++++++ batchFilter +++++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++++ normalizeRLM() +++++++++++++++++++++++++++++++++++++++
normalizeRLM <- function(elist=NULL, controls="internal", output.path=NULL){
    
    if(elist$array.type != "ProtoArray"){
        stop(paste0("ERROR: Wrong microarray type! RLM normalization is", 
        " currently implemented only for 'array.type=ProtoArray'.\n",
        "If you need this normalization approach for other types of protein",
        " microarrays please contact the PAA maintainer."))
    }
    
    controls.elist <- list(E = elist$C, genes = elist$cgenes)  
        
    if(controls == "internal"){
        controls.elist$E <-
          controls.elist$E[grep("(^HumanIg|^Anti-HumanIg)",
          controls.elist$genes$Name),]
        controls.elist$genes <-
          controls.elist$genes[grep("(^HumanIg|^Anti-HumanIg)",
          controls.elist$genes$Name),]
          
    }else if(controls == "external"){
        controls.elist$E <-
          controls.elist$E[grep("^V5", controls.elist$genes$Name),]
        controls.elist$genes <-
          controls.elist$genes[grep("^V5", controls.elist$genes$Name),]
    }else if(controls == "both"){
        controls.elist$E <-
          controls.elist$E[grep("(^HumanIg|^Anti-HumanIg|^V5)", 
          controls.elist$genes$Name),]
        controls.elist$genes <-
          controls.elist$genes[grep("(^HumanIg|^Anti-HumanIg|^V5)", 
          controls.elist$genes$Name),]
    }else{
        controls.elist$E <-
          try(controls.elist$E[grep(controls, controls.elist$genes$Name),],
          TRUE)
        controls.elist$genes <-
          try(controls.elist$genes[grep(controls, controls.elist$genes$Name),],
          TRUE)
        if(nrow(controls.elist$E) == 0){
            stop("Error! Please check your regular expression!\n")
        }
    }
    
    controls.elist$E <- log2(controls.elist$E)
    rownames(controls.elist$E) <- controls.elist$genes$Name

    contr.names <- unique(rownames(controls.elist$E))
    contr.names.len <- length(contr.names)
    contr.mapping <- matrix(nrow=contr.names.len, ncol=1)
    rownames(contr.mapping) <- contr.names
    contr.mapping[,1] <- 1:contr.names.len
    
    p.features <- nrow(elist$E)
    p.controls <- nrow(controls.elist$E)
    n.arrays <- ncol(elist$E)
    n.blocks <- max(controls.elist$genes$Block)
    
    y <- c(controls.elist$E)
    dummies <-
      matrix(0, ncol={n.arrays+n.blocks+contr.names.len-3}, nrow=length(y))
    rnames <- vector(length=length(y))
    a.cols <- paste("a",1:{n.arrays-1},sep="")
    b.cols <- paste("b",1:{n.blocks-1},sep="")
    t.cols <- paste("t",1:{contr.names.len-1},sep="")
    colnames(dummies) <- c(a.cols,b.cols,t.cols)
    
    a.params <- matrix(nrow=n.arrays, ncol=2)
    rownames(a.params) <- colnames(controls.elist)
    b.params <- matrix(nrow=n.blocks, ncol=2)


    idx <- 1
    for(i in 1:n.arrays){
        a.tmp <- paste("a",i,sep="")
        a.params[i,1] <- a.tmp
        for(j in 1:p.controls){
            rnames[idx] <- rownames(controls.elist$E)[j]
            b.idx <- controls.elist$genes$Block[j]
            b.tmp <- paste("b",b.idx,sep="")
            b.params[b.idx,1] <- b.tmp
            t.idx <- contr.mapping[rownames(controls.elist$E)[j], 1]
            t.tmp <- paste("t",t.idx,sep="")
            if(i == n.arrays){
                dummies[idx,a.cols] <- -1
            }else{
                dummies[idx,a.tmp] <- 1
            }
            
            if(b.idx == n.blocks){
                dummies[idx,b.cols] <- -1
            }else{
                dummies[idx,b.tmp] <- 1
            }
            
            if(t.idx == contr.names.len){
                dummies[idx,t.cols] <- -1
            }else{
                dummies[idx,t.tmp] <- 1
            }
            
            idx <- idx+1
        }
    }
    dummies <- cbind(y,dummies)
    rownames(dummies) <- rnames



    e <- rnorm(length(y),0,1)
    rlm.result <- rlm(y~.+e, data=data.frame(dummies))


    
    a.params[-n.arrays,2] <- rlm.result$coefficients[a.params[-n.arrays,1]]
    a.params.sum <- sum(as.numeric(a.params[-n.arrays,2]))
    a.params[n.arrays,2] <- -a.params.sum
    
    b.params[-n.blocks,2] <- rlm.result$coefficients[b.params[-n.blocks,1]]
    b.params.sum <- sum(as.numeric(b.params[-n.blocks,2]))
    b.params[n.blocks,2] <- -b.params.sum



    a <- t(matrix(rep(as.numeric(a.params[,2]),p.features), ncol=p.features))
    b <- matrix(mapply(function(i){
      as.numeric(b.params[elist$genes$Block[i],2])},
      rep(1:p.features, n.arrays)),
      ncol=n.arrays)
    normalized <- matrix(mapply(function(x,y,z){x-y-z}, log2(elist$E), a, b),
      ncol=n.arrays)
    colnames(normalized) <- colnames(elist$E) 
    result <- elist
    result$E <- normalized

    a <- t(matrix(rep(as.numeric(a.params[,2]),p.controls), ncol=p.controls))
    b <- matrix(mapply(function(i){
      as.numeric(b.params[controls.elist$genes$Block[i],2])},
      rep(1:p.controls, n.arrays)),
      ncol=n.arrays)
    contr.normalized <- matrix(mapply(function(x,y,z){x-y-z},
      controls.elist$E, a, b),
      ncol=n.arrays)
    colnames(contr.normalized) <- colnames(elist$E)



    if(!is.null(output.path)){
        dir.create(paste(output.path, "/rlm", sep=""))

        tiff(paste(output.path,"/rlm/rlm_raw.tiff",sep=""), width=4000,
          height=2000, pointsize=10, compression="lzw", res=300)
            boxplot(log2(elist$E), col="gray", cex=0.6, cex.axis=0.6, las=3)
        dev.off()
        
        tiff(paste(output.path,"/rlm/rlm_normalized.tiff",sep=""), width=4000,
          height=2000, pointsize=10, compression="lzw", res=300)
            boxplot(normalized, col="gray", cex=0.6, cex.axis=0.6, las=3)
        dev.off()

        tiff(paste(output.path,"/rlm/rlm_controls.tiff",sep=""), width=4000,
          height=2000, pointsize=10, compression="lzw", res=300)
            boxplot(controls.elist$E, col="gray", cex=0.6, cex.axis=0.6, las=3)
        dev.off()
        
        tiff(paste(output.path,"/rlm/rlm_controls_normalized.tiff",sep=""),
          width=4000, height=2000, pointsize=10, compression="lzw", res=300)
            boxplot(contr.normalized, col="gray", cex=0.6, cex.axis=0.6, las=3)
        dev.off()
    }

    return(result)
}
#+++++++++++++++++++++++++++ normalizeRLM() ++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++ normalizeArrays +++++++++++++++++++++++++++++++++++
normalizeArrays <- function(elist=NULL, method="quantile",
  cyclicloess.method="pairs", controls="internal", group1=NULL, group2=NULL,
  output.path=NULL){
    
    if(is.null(elist)) {
        stop("ERROR: Not all mandatory arguments have been defined!")
    }
    if(!is.null(group1) && !is.null(group2)) {
        if(method=="quantile"){
            tmp1 <-
              normalizeBetweenArrays(elist[,group1], method=method, ties=FALSE)
            tmp2 <-
              normalizeBetweenArrays(elist[,group2], method=method, ties=FALSE)
            elist.normalized <- cbind(tmp1, tmp2)
        }else if(method=="cyclicloess"){
            tmp1 <- normalizeBetweenArrays(elist[,group1], method=method,
              cyclicloess.method)  #cyclicloess.method: ("pairs"|"fast"|"affy")
            tmp2 <- normalizeBetweenArrays(elist[,group2], method=method,
              cyclicloess.method)
            elist.normalized <- cbind(tmp1, tmp2)
        }else if(method=="vsn"){
            tmp1 <- normalizeVSN(elist[,group1])
            tmp2 <- normalizeVSN(elist[,group2])
            elist.normalized <- cbind(tmp1, tmp2)
        }else if(method=="rlm"){
            tmp1 <- normalizeRLM(elist=elist[,group1],controls=controls,
              output.path=output.path)
            tmp2 <- normalizeRLM(elist=elist[,group2],controls=controls,
              output.path=output.path)
            elist.normalized <- cbind(tmp1, tmp2)
        }
    }else{
        if(method=="quantile"){
            elist.normalized <-
              normalizeBetweenArrays(elist, method=method, ties=FALSE)
        }else if(method=="cyclicloess"){
            elist.normalized <-
              normalizeBetweenArrays(elist, method=method, cyclicloess.method)
              #cyclicloess.method: ("pairs"|"fast"|"affy")
        }else if(method=="vsn"){
            elist.normalized <- normalizeVSN(elist)
        }else if(method=="rlm"){
            elist.normalized <-
              normalizeRLM(elist=elist,controls=controls,
              output.path=output.path)
        }
    }
    return(elist.normalized)
}
#+++++++++++++++++++++++ normalizeArrays +++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++ plotNormMethods +++++++++++++++++++++++++++++++++++++++
plotNormMethods <- function(elist=NULL, include.rlm=FALSE, controls="internal",
  output.path=NULL){
    
    if(is.null(elist)) {
        stop("ERROR: 'elist' has not been defined!")
    }
    
    if(!is.null(output.path)){
        dir.create(paste(output.path, "/norm_methods_plots", sep=""))
    }
    elist.quantile <- normalizeArrays(elist=elist, method="quantile")
    elist.cyclicloess.pairs <-
      normalizeArrays(elist=elist, method="cyclicloess",
      cyclicloess.method="pairs")   #("pairs"|"fast"|"affy")
    elist.vsn <- normalizeArrays(elist=elist, method="vsn")
    if(include.rlm){
        elist.rlm <-
          normalizeArrays(elist=elist, method="rlm", controls=controls,
          output.path=NULL)
    }

    if(is.null(output.path)){
        par(ask=TRUE)
            boxplot(log2(elist$E), col="gray", cex=0.6, cex.axis=0.6, las=3,
              main="raw values")
            boxplot(elist.quantile$E, col="gray", cex=0.6, cex.axis=0.6, las=3,
              main="quantile normalization")
            boxplot(elist.cyclicloess.pairs$E, col="gray", cex=0.6,
              cex.axis=0.6, las=3, main="cyclic loess normalization")
            boxplot(elist.vsn$E, col="gray", cex=0.6, cex.axis=0.6, las=3,
              main="vsn normalization")
            if(include.rlm){
                boxplot(elist.rlm$E, col="gray", cex=0.6, cex.axis=0.6, las=3,
                  main="rlm normalization")
                dev.off()
            }else{
                dev.off()
            }
    }else{
        tiff(paste(output.path,"/norm_methods_plots/boxplots_raw.tiff",sep=""),
          width=4000, height=2000, pointsize=10, compression="lzw", res=300)
            boxplot(log2(elist$E),col="gray", cex=0.6, cex.axis=0.6, las=3)
        dev.off()
        tiff(paste(output.path,"/norm_methods_plots/boxplots_quantile.tiff",
          sep=""), width=4000, height=2000, pointsize=10, compression="lzw",
          res=300)
            boxplot(elist.quantile$E,col="gray", cex=0.6, cex.axis=0.6, las=3)
        dev.off()
        tiff(paste(output.path,
          "/norm_methods_plots/boxplots_cyclicloess_pairs.tiff",sep=""),
          width=4000, height=2000, pointsize=10, compression="lzw", res=300)
            boxplot(elist.cyclicloess.pairs$E,col="gray", cex=0.6, cex.axis=0.6,
              las=3)
        dev.off()
        tiff(paste(output.path,"/norm_methods_plots/boxplots_vsn.tiff",sep=""),
          width=4000, height=2000, pointsize=10, compression="lzw", res=300)
            boxplot(elist.vsn$E,col="gray", cex=0.6, cex.axis=0.6, las=3)
        dev.off()
        if(include.rlm){
            tiff(paste(output.path,"/norm_methods_plots/boxplots_rlm.tiff",
              sep=""), width=4000, height=2000, pointsize=10, compression="lzw",
              res=300)
                boxplot(elist.rlm$E,col="gray", cex=0.6, cex.axis=0.6, las=3)
            dev.off()
        }
    }
}
#+++++++++++++++++++++++ plotNormMethods +++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++ preselect +++++++++++++++++++++++++++++++++++++++

preselect <- function(elist=NULL, columns1=NULL, columns2=NULL, label1="A",
  label2="B", log=NULL, discard.threshold=0.5, fold.thresh=1.5,
  discard.features=TRUE, mMs.above=1500, mMs.between=400, mMs.matrix1=NULL,
  mMs.matrix2=NULL, method=NULL){
    
    if(is.null(log)) {
        stop("preselect - ERROR: Not all mandatory arguments defined!")
    }
    
    if(method=="mMs" && log==TRUE){
        elist$E <- 2^(elist$E)
    }else if(method=="tTest" && log==FALSE){
        elist$E <- log2(elist$E)
    }else if(method=="mrmr" && log==FALSE){
        elist$E <- log2(elist$E)
    }
    
    if(method=="mMs"){
        mMs2(elist=elist, columns1=columns1, columns2=columns2, label1=label1,
          label2=label2, above=mMs.above, between=mMs.between,
          discard.threshold=discard.threshold, fold.thresh=fold.thresh,
          discard.features=discard.features, mMs.matrix1=mMs.matrix1,
          mMs.matrix2=mMs.matrix2)
    }else if(method=="tTest"){
        tTest(elist=elist, columns1=columns1, columns2=columns2, label1=label1,
          label2=label2, discard.threshold=discard.threshold,
          fold.thresh=fold.thresh, discard.features=discard.features)
    }else if(method=="mrmr"){
        if(discard.threshold == 0.5){discard.threshold <- 300}
        mrmr(elist=elist, columns1=columns1, columns2=columns2, label1=label1,
          label2=label2, discard.threshold=discard.threshold,
          discard.features=discard.features)
    }else{
        stop("preselect(): Error - Unknown preselection method!!!")
    }
}

#+++++++++++++++++++++++++++ preselect +++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++ selectFeatures.frequency.cv +++++++++++++++++++++++
selectFeatures.frequency.cv <- function(elist=NULL, columns1=NULL,
  columns2=NULL, n1=NULL, n2=NULL, label1="A", label2="B", cutoff=10,
  selection.method="rf.rfe", preselection.method="mMs", subruns=100,
  candidate.number=300, above=1500, between=400, mMs.matrix1=NULL,
  mMs.matrix2=NULL, panel.selection.criterion="accuracy",
  importance.measure="MDA", ntree=500, mtry=NULL, plot=TRUE, output.path=NULL,
  verbose=FALSE){
    
    if(is.null(elist) || is.null(columns1) || is.null(columns2) || is.null(n1)
      || is.null(n2)) {
        stop("ERROR: Not all mandatory arguments have been defined!")
    }
    classes_test <- factor(c(rep(label1, floor(n1/3)),
      rep(label2, floor(n2/3))))
    classes_train <- factor(c(rep(label1, {n1-floor(n1/3)}),
      rep(label2, {n2-floor(n2/3)})))
    if(!is.null(output.path)){
        dir.create(output.path)
        cat(x="Iteration\tAccuracy\tSensitivity\tSpecificity\tProteins\n",
          file=paste(output.path, "/selectFeaturesSubruns.txt", sep=""),
          append=FALSE)
    }
    biomarker.list <- list()
    row.number <- nrow(elist$E)
    for(it in 1:subruns) {
        if(preselection.method == "mMs"){
            mMsFS.results <- mMsFS2(elist, columns1, columns2, label1, label2,
              above=above, between=between, mMs.matrix1=mMs.matrix1,
              mMs.matrix2=mMs.matrix2)
        }else if(preselection.method == "tTest"){
            mMsFS.results <- tTestFS(elist, columns1, columns2, label1, label2)
        }else if(preselection.method == "mrmr"){
            mMsFS.results <- mrmrFS(elist=elist, columns1=columns1,
              columns2=columns2, label1=label1, label2=label2,
              candidate.number=candidate.number)
        }else if(preselection.method == "none"){
            mMsFS.results <- noFS(elist, columns1, columns2, label1, label2)
            candidate.number <- row.number
        }else{
            stop(paste("selectFeatures.frequency.cv()-Error:", 
              " Unknown pre-selection method!", sep=""))
        }

        traincols1 <- mMsFS.results$traincols1
        traincols2 <- mMsFS.results$traincols2
        traincols <- c(traincols1,traincols2)
        testcols <- mMsFS.results$testcols
        traincols1.length <- length(traincols1)
        traincols2.length <- length(traincols2)
        prospector <- mMsFS.results$prospector
        trainsample <- mMsFS.results$trainsample
        testsample <- mMsFS.results$testsample

        #if(!is.null(output.path)){
        #   if(preselection.method == "mMs"){
        #       write.table(x=trainsample[1:candidate.number,-1],
        #         file=paste(output.path, "/subrun", it, "/trainsample_", it,
        #         ".txt", sep=""), sep="\t", eol="\n", row.names=FALSE,
        #         quote=FALSE)
        #       write.table(x=testsample[1:candidate.number,-1],
        #         file=paste(output.path, "/subrun", it, "/testsample_", it,
        #         ".txt", sep=""), sep="\t", eol="\n", row.names=FALSE,
        #         quote=FALSE)
        #   }else if(preselection.method == "tTest"){
        #       write.table(x=trainsample[1:candidate.number,-1],
        #         file=paste(output.path, "/subrun", it, "/trainsample_", it,
        #         ".txt", sep=""), sep="\t", eol="\n", row.names=FALSE,
        #         quote=FALSE)
        #       write.table(x=testsample[1:candidate.number,-1],
        #         file=paste(output.path, "/subrun", it, "/testsample_", it,
        #         ".txt", sep=""), sep="\t", eol="\n", row.names=FALSE,
        #         quote=FALSE)
        #   }else if(preselection.method == "mrmr"){
        #       write.table(x=trainsample[1:candidate.number,-1],
        #         file=paste(output.path, "/subrun", it, "/trainsample_", it,
        #         ".txt", sep=""), sep="\t", eol="\n", row.names=FALSE,
        #         quote=FALSE)
        #       write.table(x=testsample[1:candidate.number,-1],
        #         file=paste(output.path, "/subrun", it, "/testsample_", it,
        #         ".txt", sep=""), sep="\t", eol="\n", row.names=FALSE,
        #         quote=FALSE)
        #   }else if(preselection.method == "none"){
        #       write.table(x=trainsample[,-1], file=paste(output.path,
        #         "/subrun", it, "/trainsample_", it, ".txt", sep=""),
        #         sep="\t", eol="\n", row.names=FALSE, quote=FALSE)
        #       write.table(x=testsample[,-1], file=paste(output.path,
        #         "/subrun", it, "/testsample_", it, ".txt", sep=""), sep="\t",
        #         eol="\n", row.names=FALSE, quote=FALSE)
        #   }else{
        #      stop(paste("selectFeatures.frequency.cv()-Error:", 
        #        " Unknown pre-selection method!", sep=""))
        #   }
        #}

        if(selection.method == "rf.rfe"){
            best.list <- rf.rfe(it,trainsample[1:candidate.number,-1],
              output.path,label1,label2,traincols1.length,traincols2.length,
              panel.selection.criterion=panel.selection.criterion,
              importance.measure=importance.measure,ntree=ntree,mtry=mtry,
              verbose=verbose)
        }else if(selection.method == "svm.rfe"){
            best.list <- svm.rfe(it, trainsample[1:candidate.number,-1],
              output.path=output.path,label1=label1,label2=label2,
              n1=traincols1.length,n2=traincols2.length,
              panel.selection.criterion=panel.selection.criterion,
              verbose=verbose)
        }else if(selection.method == "rj.rfe"){
            best.list <- rj.rfe(iteration=it,trainsample[1:candidate.number,-1],
                output.path=output.path,label1=label1,label2=label2,
                n1=traincols1.length,n2=traincols2.length,
                panel.selection.criterion=panel.selection.criterion,
                importance.measure=importance.measure,ntree=ntree,
                mtry=mtry,verbose=verbose) 
        }else{
            stop(paste("selectFeatures.frequency.cv()-Error:", 
              " Unknown selection method!", sep=""))
        }
        best.set <- which(rownames(testsample) %in% best.list)
        traindata <- data.frame(best.list,rbind(elist$E[best.set,traincols]))
        testdata <- data.frame(best.list,rbind(elist$E[best.set,testcols]))
        names(testdata)[1] <- names(traindata)[1] <- "BRC"
        if(selection.method == "rf.rfe"){
            fin.classi.results <- final.classify.rf(trainset=traindata, 
              testset=testdata, classes_train=classes_train, 
              classes_test=classes_test, label1=label1, label2=label2, 
              cwd=output.path, iteration=it, plot=FALSE, performance=plot,
              verbose=verbose)
        }else if(selection.method == "rj.rfe"){
            fin.classi.results <- final.classify.rf(trainset=traindata,
              testset=testdata, classes_train=classes_train,
              classes_test=classes_test, label1=label1, label2=label2,
              cwd=output.path, iteration=it, plot=FALSE, performance=plot,
              verbose=verbose)
        }else if(selection.method == "svm.rfe"){
            fin.classi.results <- final.classify.svm(trainset=traindata,
              testset=testdata, classes_train=classes_train,
              classes_test=classes_test, label1=label1, label2=label2,
              cwd=output.path, iteration=it, plot=FALSE, performance=plot,
              verbose=verbose)
        }else{
            stop(paste("selectFeatures.frequency.cv()-Error:", 
              " Unknown selection method!", sep=""))
        }
        if(!is.null(output.path)){
            cat(it, fin.classi.results$accuracy, fin.classi.results$sensitivity,
              fin.classi.results$specificity, best.list, file=paste(output.path,
              "/selectFeaturesSubruns.txt", sep=""), sep="\t", append=TRUE)
            cat("\n", file=paste(output.path, "/selectFeaturesSubruns.txt",
              sep=""), sep="\t", append=TRUE)
        }
        biomarker.list <-
          biomarkers.update(biomarker.list=biomarker.list,
          update.vector=rownames(testsample)[best.set])

        if(!is.null(output.path)){
            biomarker.stat <- data.frame(names(biomarker.list),
              unlist(biomarker.list),unlist(biomarker.list)/subruns)
            colnames(biomarker.stat) <- c("BRC","Count","Percentage")
            write.table(x=biomarker.stat[order(-biomarker.stat$"Count"),],
              file=paste(output.path, "/biomarkerStat.txt", sep=""), sep="\t",
              eol="\n", row.names=FALSE, quote=FALSE)
        }
    }
    if(plot && !is.null(output.path)){
        preds_total <- as.matrix(read.table(paste(output.path,
          "/preds_total.txt", sep=""), header=FALSE, sep="\t", na.strings="NA",
          dec=".", strip.white=TRUE, quote=""))
        classes_total <- as.matrix(read.table(paste(output.path,
          "/classes_total.txt", sep=""), header=FALSE, sep="\t",
          na.strings="NA", dec=".", strip.white=TRUE, quote=""))
        perf <- performance(prediction(as.vector(preds_total[1,]),
          as.vector(classes_total[1,])), 'tpr', 'fpr')
        tiff(paste(output.path,"/roc_all.tif",sep=""), width=2000, height=2000,
          pointsize=15, compression="lzw", res=300)
            plot(perf, plot=FALSE)
            auc <- 0
            for(i in 1:nrow(preds_total)) {
                perf <- performance(prediction(as.vector(preds_total[i,]),
                  as.vector(classes_total[i,])), 'tpr', 'fpr')
                auc <- auc + as.numeric(performance(prediction(as.vector(
                  preds_total[i,]),
                  as.vector(classes_total[i,])), 'auc')@y.values)
                plot(perf, col="red", add=TRUE, lwd=2)
            }
            title(main=paste(label1, " vs. ", label2, " - average AUC: ",
              round(auc / nrow(preds_total),4), sep=""))
        dev.off()
    }
    
    biomarker.stat <- data.frame(names(biomarker.list),unlist(biomarker.list),
      unlist(biomarker.list)/subruns, stringsAsFactors=FALSE)
    colnames(biomarker.stat) <- c("BRC","Count","Percentage")
    output <- biomarker.stat[order(-biomarker.stat$"Count"),]
    if(!is.null(output.path)){
        write.table(x=output, file=paste(output.path, "/biomarkerStat.txt",
          sep=""), sep="\t", eol="\n", row.names=FALSE, quote=FALSE)
    }
    return(na.omit(output[1:cutoff,]))
}
#+++++++++++++++++ selectFeatures.frequency.cv +++++++++++++++++++++++++++++++++

#++++++++++++++++++++ selectFeatures.frequency +++++++++++++++++++++++++++++++++
selectFeatures.frequency <- function(elist=NULL, n1=NULL, n2=NULL, label1="A",
  label2="B", cutoff=10, selection.method="rf.rfe", preselection.method="mMs",
  subruns=100, k=10, candidate.number=300, above=1500, between=400,
  panel.selection.criterion="accuracy", importance.measure="MDA", ntree=500,
  mtry=NULL, plot=FALSE, output.path=NULL, verbose=FALSE){
    
    if(is.null(elist) || is.null(n1) || is.null(n2)) {
        stop("ERROR: Not all mandatory arguments have been defined!")
    }
    
    trainsample_urn <- colnames(elist$E)
    testsample_urn1 <- colnames(elist$E[,1:n1])
    testsample_urn2 <- colnames(elist$E[,{n1+1}:{n1+n2}])
    testsamplesize1 <- floor(n1/k)
    testsamplesize2 <- floor(n2/k)
    testsamplesize <- testsamplesize1+testsamplesize2
    trainsamplesize <- n1+n2-testsamplesize1-testsamplesize2
    classes_test <- 
      factor(c(rep(label1, testsamplesize1), rep(label2, testsamplesize2)))
    classes_train <- 
      factor(c(rep(label1, {n1-testsamplesize1}),
      rep(label2, {n2-testsamplesize2})))

    crossValTest <- matrix(nrow=k,ncol=testsamplesize)
    crossValTrain <- matrix(nrow=k,ncol=trainsamplesize)
    selectFeatures.results <- list()
    for (i in 1:k){
        testsample1 <- sample(testsample_urn1, testsamplesize1, replace=FALSE)
        testsample2 <- sample(testsample_urn2, testsamplesize2, replace=FALSE)
        testsample_urn1 <- testsample_urn1[!(testsample_urn1 %in% testsample1)]
        testsample_urn2 <- testsample_urn2[!(testsample_urn2 %in% testsample2)]
        crossValTest[i,] <- testsample <- c(testsample1, testsample2)
        crossValTrain[i,] <- trainsample_urn[!(trainsample_urn %in% testsample)]
    }

    if(preselection.method == "mMs"){
        mMs.matrix1 <- mMsMatrix({n1-testsamplesize1}, {n2-testsamplesize2})
        mMs.matrix2 <- mMsMatrix({n2-testsamplesize2}, {n1-testsamplesize1})
    }else{
        mMs.matrix1 <- NULL
        mMs.matrix2 <- NULL
    }
    tmp.brc <- matrix(nrow=nrow(elist),ncol=2)
    tmp.brc[,1] <- 
      paste(elist$genes$Block,elist$genes$Row,elist$genes$Column)
    tmp.brc[,2] <- 1:nrow(elist)
    selected.features <- c()
    avg.accuracy <- 0
    avg.sensitivity <- 0
    avg.specificity <- 0
    if(!is.null(output.path)){
        cat(x="Accuracy\tSensitivity\tSpecificity\tTP\tFP\tTN\tFN\tProteins\n",
          file=paste(output.path, "/selectFeaturesCV.txt", sep=""),
          append=FALSE)
    }
    for (i in 1:k){
        if(!is.null(output.path)){
            tmp.output.path <- paste(output.path, "/cv_iteration", i,  sep="")
        } else {
            tmp.output.path <- NULL
        }
        selectFeatures.results[[i]] <-
          selectFeatures.frequency.cv(elist[,crossValTrain[i,]],
          crossValTrain[i,1:{n1-testsamplesize1}],
          crossValTrain[i,{n1-testsamplesize1+1}:trainsamplesize],
          n1=n1-testsamplesize1, n2=n2-testsamplesize2, label1=label1,
          label2=label2, cutoff=cutoff, selection.method=selection.method,
          preselection.method=preselection.method, subruns=subruns,
          candidate.number=candidate.number, above=above, between=between,
          mMs.matrix1=mMs.matrix1, mMs.matrix2=mMs.matrix2,
          panel.selection.criterion=panel.selection.criterion,
          importance.measure=importance.measure, ntree=ntree, mtry=mtry,
          plot=plot, output.path=tmp.output.path, verbose=verbose)
          
        selected.features <- 
          c(selected.features, selectFeatures.results[[i]]$BRC)
        indexes <-
          as.numeric(tmp.brc[tmp.brc[,1] %in% selectFeatures.results[[i]]$BRC,,
          drop=FALSE][,2])
        traindata <- data.frame(paste(elist$genes$Block[indexes],
          elist$genes$Row[indexes],elist$genes$Column[indexes]),
          rbind(elist$E[indexes,crossValTrain[i,]]))
        testdata <- data.frame(paste(elist$genes$Block[indexes],
          elist$genes$Row[indexes],elist$genes$Column[indexes]),
          rbind(elist$E[indexes,crossValTest[i,]]))
        names(testdata)[1] <- names(traindata)[1] <- "BRC"
        
        if(selection.method == "rf.rfe"){
            fin.classi.results <- final.classify.rf(trainset=traindata,
              testset=testdata, classes_train=classes_train,
              classes_test=classes_test, label1=label1, label2=label2,
              cwd=output.path, iteration=i, plot=plot, performance=plot)
        }else if(selection.method == "rj.rfe"){
            fin.classi.results <- final.classify.rf(trainset=traindata,
              testset=testdata, classes_train=classes_train,
              classes_test=classes_test, label1=label1, label2=label2,
              cwd=output.path, iteration=i, plot=plot, performance=plot)
        }else if(selection.method == "svm.rfe"){
            fin.classi.results <- final.classify.svm(trainset=traindata,
              testset=testdata, classes_train=classes_train,
              classes_test=classes_test, label1=label1, label2=label2,
              cwd=output.path, iteration=i, plot=plot, performance=plot)
        }else{
            stop(paste("selectFeatures.frequency()-Error:", 
              " Unknown selection method!", sep=""))
        }
        avg.accuracy <- avg.accuracy + fin.classi.results$accuracy
        avg.sensitivity <- avg.sensitivity + fin.classi.results$sensitivity
        avg.specificity <- avg.specificity + fin.classi.results$specificity
        if(!is.null(output.path)){
            cat(fin.classi.results$accuracy, fin.classi.results$sensitivity,
              fin.classi.results$specificity, fin.classi.results$tp,
              fin.classi.results$fp, fin.classi.results$tn,
              fin.classi.results$fn, selectFeatures.results[[i]]$BRC,
              file=paste(output.path, "/selectFeaturesCV.txt", sep=""),
              sep="\t", append=TRUE)
            cat("\n", file=paste(output.path, "/selectFeaturesCV.txt", sep=""),
              sep="\t", append=TRUE)
        }
        message("selectFeatures - current cross validation set:", i,
          " of ", k, "\n", sep="  ")
        message("selectFeatures - selected features: ")
        message(selectFeatures.results[[i]]$BRC, sep=", ")
        message("\n")
        message("selectFeatures - classification results - accuracy: ",
          fin.classi.results$accuracy, ", sensitivity: ",
          fin.classi.results$sensitivity, ", specificity: ",
          fin.classi.results$specificity, "\n", sep="")
        message("selectFeatures - current set of non-redundant features: ")
        message(unique(selected.features), sep=", ")
        message("\n\n")
    }
    avg.accuracy <- avg.accuracy / k
    avg.sensitivity <- avg.sensitivity / k
    avg.specificity <- avg.specificity / k
    message(paste("selectFeatures",
      " - average classification resultds - accuracy: ", avg.accuracy,
      ", sensitivity: ", avg.sensitivity, ", specificity: ", avg.specificity,
      sep=""), "\n\n\n")
    return(list(accuracy=avg.accuracy, sensitivity=avg.sensitivity,
      specificity=avg.specificity, features=unique(selected.features),
      all.results=selectFeatures.results))
}
#++++++++++++++++++++++ selectFeatures.frequency +++++++++++++++++++++++++++++++

#++++++++++++++++++++++++ subsample.data +++++++++++++++++++++++++++++++++++++++
subsample.data <- function(columns=NULL, n1=NULL, n2=NULL, k=10,
  subsamples=10) {
    
    if(is.null(columns) || is.null(n1) || is.null(n2)) {
        stop("ERROR: Not all mandatory arguments have been defined!")
    }
    n <- n1+n2
    test.col.size1 <- ceiling(n1/k)
    test.col.size2 <- ceiling(n2/k)
    test.col.size <- test.col.size1 + test.col.size2
    train.col.size <- n - test.col.size
    columns1 <- columns[1:n1]
    columns2 <- columns[{n1+1}:n]
    test.columns <- matrix(nrow=test.col.size, ncol=subsamples)
    train.columns <- matrix(nrow=train.col.size, ncol=subsamples)

    for(i in 1:subsamples){
        testcols1 <-
          columns1[columns1 %in% sample(columns1,test.col.size1,replace=FALSE)]
        testcols2 <-
          columns2[columns2 %in% sample(columns2,test.col.size2,replace=FALSE)]
        testcols <- c(testcols1, testcols2)
        
        traincols1 <- columns1[!(columns1 %in% testcols1)]
        traincols2 <- columns2[!(columns2 %in% testcols2)]
        traincols <- c(traincols1,traincols2)
        
        test.columns[,i]  <- testcols
        train.columns[,i] <- traincols
    }
    
    return(list(train.columns=train.columns, test.columns=test.columns))
}
#+++++++++++++++++++++++++ subsample.data ++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++ bootstrap.data ++++++++++++++++++++++++++++++++++++++
bootstrap.data <- function(columns=NULL, n1=NULL, n2=NULL, bootstraps=100) {
    if(is.null(columns) || is.null(n1) || is.null(n2)) {
        stop("ERROR: Not all mandatory arguments have been defined!")
    }
    n <- n1+n2
    test.columns <- list()
    train.columns <- matrix(nrow=n, ncol=bootstraps)

    for(i in 1:bootstraps){
        traincols <- sample(columns,n,replace=TRUE)
        testcols <- columns[!(columns %in% traincols)]
        train.columns[,i] <- traincols
        test.columns[[i]]  <- testcols
    }
    return(list(train.columns=train.columns, test.columns=test.columns))
}
#+++++++++++++++++++++++++++ bootstrap.data ++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++ svm.rfe.ensemble ++++++++++++++++++++++++++++++++++++
svm.rfe.ensemble <- function(x=NULL, label1="A", label2="B", n1=NULL, n2=NULL){
    if(is.null(x) || is.null(n1) || is.null(n2)) {
        stop("ERROR: Not all mandatory arguments have been defined!")
    }

    data.table <- t(x)
    classes <- factor(gsub('\\d+', '', colnames(x)))
    p <- ncol(data.table)
    survivingFeaturesIndexes <- seq(1:p)
    final.ranking <- c()

    while({p <- length(survivingFeaturesIndexes)} > 0){
        train.dat <- data.table[, survivingFeaturesIndexes, drop=FALSE]

        #~~~~~~~~~~~~~~~CLASSIFIER BEGIN~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        svm.model <- svm(train.dat, classes, cost=10, type="C-classification",
          kernel="linear")
        w <- t(svm.model$coefs)%*%svm.model$SV
        rankingCriteria <- w*w
        #~~~~~~~~~~~~~~~CLASSIFIER END~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if(p >= 10){
            tenth <- floor(p*0.1)
        }else{
            tenth <- ceiling(p*0.1)
        }
        ranking <- sort(rankingCriteria, index.return = TRUE)$ix
        if(tenth > 1){
            final.ranking <- c(survivingFeaturesIndexes[ranking[1:tenth]],
              final.ranking)
            survivingFeaturesIndexes <-
              survivingFeaturesIndexes[-ranking[1:tenth]]
        }else if(tenth == 1){
            final.ranking <-
              c(survivingFeaturesIndexes[ranking[1]], final.ranking)
            survivingFeaturesIndexes <- survivingFeaturesIndexes[-ranking[1]]
        }
    }
    result <- rep(0, length(final.ranking))
    for(i in 1:length(final.ranking)){
        result[final.ranking[i]] <- i
    }
    return(result)
}
#++++++++++++++++++++++++ svm.rfe.ensemble +++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++ classify.svm.ensemble +++++++++++++++++++++++++++++++++++
classify.svm.ensemble <- function(iteration=NULL, trainset=NULL, testset=NULL,
  label1="A", label2="B", plot=FALSE, output.path=NULL, verbose=verbose){
    
    if(is.null(iteration) || is.null(trainset) || is.null(testset)) {
        stop("ERROR: Not all mandatory arguments have been defined!")
    }

    classes_train <- factor(gsub('\\d+', '', colnames(trainset)))
    classes_test <- factor(gsub('\\d+', '', colnames(testset)))
    train.dat <- t(as.matrix(trainset))
    test.dat <- t(as.matrix(testset))

    #tuneobj <-
    #  tune.svm(x=train.dat, y=classes_train, gamma=10^(-9:-1), cost=10^(0:4))
    tuneobj <- tune.svm(x=train.dat, y=classes_train, cost=10^(0:4))
    #bestGamma <- tuneobj$best.parameters[[1]]
    #bestCost <- tuneobj$best.parameters[[2]]
    bestCost <- tuneobj$best.parameters[[1]]
    #model.svm <-
    #  svm(train.dat, classes_train, probability=TRUE, cost=10,
    #  type="C-classification", kernel="linear")
    model.svm <- svm(train.dat, classes_train, probability=TRUE, cost=bestCost,
      type="C-classification", kernel="linear")
    #model.svm <- svm(train.dat, classes_train, probability=TRUE, cost=bestCost,
    #  gamma=bestGamma, type="C-classification", kernel="radial")
    pred.svm <- predict(object=model.svm, newdata=test.dat,
      decision.values = TRUE, probability = TRUE)
    confusion_matrix <- table(observed = factor(classes_test,
      levels=sort(c(label1, label2))), predicted = pred.svm)

    TP <- confusion_matrix[1,1]
    FP <- confusion_matrix[2,1]
    TN <- confusion_matrix[2,2]
    FN <- confusion_matrix[1,2]
    if(verbose){
        print(confusion_matrix)
    }

    pred.to.roc <- attr(pred.svm, "probabilities")[,1]
    labs <- colnames(attr(pred.svm, "probabilities"))
    pred.rocr <- prediction(pred.to.roc, classes_test,
      label.ordering=c(labs[2], labs[1]))
    #pred.rocr <- prediction(pred.to.roc, classes_test,
    #  label.ordering=c(label2, label1))
    perf.rocr.tpr <- performance(pred.rocr, measure="tpr", x.measure="fpr")
    perf.rocr.auc <- performance(pred.rocr, measure="auc")
    auc <- as.numeric(perf.rocr.auc@y.values)
    if(!is.null(output.path) && plot){
        tiff(paste(output.path,"/subsample", iteration, "/roc.tiff",sep=""),
          width=2000, height=2000, pointsize=40, compression="lzw")
            plot(perf.rocr.tpr, col="red", lwd=2, main=
              paste(label1, " vs. ", label2, " - average AUC: ", round(auc,4),
              sep=""))
        dev.off()
    }

    ACCURACY <- {TP+TN}/{TP+FP+TN+FN}
    SENSITIVITY <- TP/{TP+FN}   #=TPR
    SPECIFICITY <- TN/{TN+FP}   #=(1-FPR)
    message(paste("selectFeatures - classification results - accuracy:",
      ACCURACY, ", sensitivity: ", SENSITIVITY, ", specificity: ", SPECIFICITY,
      sep=""), "\n\n")
    results <- list(auc=auc, accuracy=ACCURACY, sensitivity=SENSITIVITY,
      specificity=SPECIFICITY, tp=TP, fp=FP, tn=TN, fn=FN)
    return(results)
}
#+++++++++++++++++++++++ classify.svm.ensemble +++++++++++++++++++++++++++++++++

#++++++++++++++++++++++++ assess.stability +++++++++++++++++++++++++++++++++++++
assess.stability <- function(x=NULL, p=NULL, subsamples=NULL) {
    if(is.null(x) || is.null(p) || is.null(subsamples)) {
        stop("ERROR: Not all mandatory arguments have been defined!")
    }
    stability <- 0
    for(i in 1:{subsamples-1}){
        for(j in {i+1}:subsamples){
            s.i <- length(x[[i]])
            s.j <- length(x[[j]])
            s.min <- min(s.i, s.j)
            r <- length(intersect(x[[i]], x[[j]]))
            kuncheva.index <- {r - {s.min^2/p}} / {s.min - {s.min^2/p}}
            stability <- stability + kuncheva.index
        }
    }
    stability <- 2*stability / {subsamples*{subsamples - 1}}
    return(stability)
}
#++++++++++++++++++++++++ assess.stability +++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++ selectFeatures.ensemble +++++++++++++++++++++++++++++++++
selectFeatures.ensemble <- function(elist=NULL, n1=NULL, n2=NULL, label1="A",
  label2="B", cutoff=10, k=10, subsamples=10, bootstraps=100, plot=FALSE,
  output.path=NULL, verbose=FALSE) {
    
    if(is.null(elist) || is.null(n1) || is.null(n2)) {
        stop("ERROR: Not all mandatory arguments have been defined!")
    }
    x <- elist$E
    rownames(x) <- paste(elist$genes$Block,elist$genes$Row,elist$genes$Column) 
    subsample.results <- subsample.data(columns=colnames(x), n1=n1, n2=n2, k=k,
      subsamples=subsamples)
    train.columns <- subsample.results$train.columns
    test.columns <- subsample.results$test.columns
    n.test1 <- ceiling(n1/k)
    n.test2 <- ceiling(n2/k)
    n.train1 <- n1 - n.test1
    n.train2 <- n2 - n.test2
    n.train <- nrow(train.columns)
    n.test <- nrow(test.columns)
    p <- nrow(x)
    subsample.features <- list()
    subsample.rankings <- matrix(nrow=p, ncol=subsamples)
    rownames(subsample.rankings) <- rownames(x)

    avg.accuracy <- 0
    avg.sensitivity <- 0
    avg.specificity <- 0
    if(!is.null(output.path)){
        cat(x=paste("iteration\tfeature number\tAUC\tAccuracy\tSensitivity\t",
          "Specificity\tTP\tFP\tTN\tFN\tfeatures\n", sep=""),
          file=paste(output.path, "/selection.results.txt", sep=""),
          append=FALSE)
    }
    for(i in 1:subsamples){
        message(paste("selectFeatures - current subsample: ", i, sep=""),
          " of ", subsamples, "\n")
        bootstrap.results <- bootstrap.data(columns=
          colnames(x[,train.columns[,i]]), n1=n.train1, n2=n.train2,
          bootstraps=100)
        bootstrap.columns <- bootstrap.results$train.columns
        #oob.columns <- bootstrap.results$test.columns
        rankings <- matrix(nrow=p, ncol=bootstraps)
        rownames(rankings) <- rownames(x)
        for(j in 1:bootstraps){
            if(verbose){message(paste("bootstrap sample: ", j, sep=""), "\n")}
            rankings[,j] <- svm.rfe.ensemble(x[,bootstrap.columns[,j]],
              label1=label1, label2=label2, n1=n.train1, n2=n.train2)
        }
        rankings <- cbind(rankings, apply(rankings, 1, sum))
        rankings[,bootstraps+1] <- rankings[,bootstraps+1] / bootstraps
        rankings.ordered <- rankings[order(rankings[,bootstraps+1]),]
        if(!is.null(output.path)){
            dir.create(paste(output.path, "/subsample", i, sep=""))
            write.table(x=rankings.ordered, file=paste(output.path,
              "/subsample", i, "/bootstrap.rankings.txt", sep=""), sep="\t",
              eol="\n", row.names=TRUE, col.names=FALSE)
        }
        subsample.rankings[,i] <- rankings[,bootstraps+1]
        
        features <- rownames(rankings.ordered)[1:cutoff]
        subsample.features[[i]] <- features
        message("selectFeatures - best features for current subsample: ")
        message(features, sep=", ")
        message("\n")
        
        classi <- classify.svm.ensemble(iteration=i, trainset=x[features,
          train.columns[,i], drop=FALSE], testset=x[features, test.columns[,i],
          drop=FALSE], label1=label1, label2=label2, plot=plot,
          output.path=output.path, verbose=verbose)
        
        avg.accuracy <- avg.accuracy + classi$accuracy
        avg.sensitivity <- avg.sensitivity + classi$sensitivity
        avg.specificity <- avg.specificity + classi$specificity
        
        if(!is.null(output.path)){
            cat(x=paste(i, cutoff, classi$auc, classi$accuracy,
              classi$sensitivity, classi$specificity, classi$tp, classi$fp,
              classi$tn, classi$fn, sep="\t"), file=paste(output.path,
              "/selection.results.txt", sep=""), append=TRUE)
            cat(x="\t", file=paste(output.path, "/selection.results.txt",
              sep=""), append=TRUE)
            cat(x=features, sep="\t", file=paste(output.path,
              "/selection.results.txt", sep=""), append=TRUE)
            cat(x="\n", file=paste(output.path, "/selection.results.txt",
              sep=""), append=TRUE)
        }
    }
    avg.accuracy <- avg.accuracy / subsamples
    avg.sensitivity <- avg.sensitivity / subsamples
    avg.specificity <- avg.specificity / subsamples
    
    stability <-
      assess.stability(x=subsample.features, p=p, subsamples=subsamples)
    if(verbose){message(paste("total stability: ", stability, sep=""), "\n")}
    subsample.rankings <-
      cbind(subsample.rankings, apply(subsample.rankings,1,mean))
    colnames(subsample.rankings) <- c(1:subsamples, "Mean")
    subsample.rankings <-
      subsample.rankings[order(subsample.rankings[,"Mean"]),]
    
    if(!is.null(output.path)){
        cat(x=paste("Total stability: ", stability, sep=""),
          file=paste(output.path, "/stability.txt", sep=""), append=FALSE)
        write.table(x=cbind(c("ID",rownames(subsample.rankings)),
          subsample.rankings), file=paste(output.path,
          "/subsample.rankings.txt", sep=""), sep="\t", eol="\n",
          row.names=FALSE, col.names=TRUE)
    }
    message("selectFeatures - finally selected features: ")
    message(rownames(subsample.rankings)[1:cutoff], sep=", ")
    message("\n")
    message(paste(
      "selectFeatures - average classification resultds - accuracy: ",
      avg.accuracy, ", sensitivity: ", avg.sensitivity, ", specificity: ",
      avg.specificity, sep=""), "\n\n\n")
   
    return(list(accuracy=avg.accuracy, sensitivity=avg.sensitivity,
      specificity=avg.specificity,
      features=rownames(subsample.rankings)[1:cutoff],
      all.results=subsample.rankings, stability=stability))
}
#++++++++++++++++++++ selectFeatures.ensemble ++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++ selectFeatures ++++++++++++++++++++++++++++++++++++++

selectFeatures <- function(elist=NULL, n1=NULL, n2=NULL, label1="A", label2="B",
  log=NULL, cutoff=10, selection.method="rf.rfe", preselection.method="mMs",
  subruns=100, k=10, subsamples=10, bootstraps=10, candidate.number=300,
  above=1500, between=400, panel.selection.criterion="accuracy",
  importance.measure="MDA", ntree=500, mtry=NULL, plot=FALSE, output.path=NULL,
  verbose=FALSE, method="frequency"){
    
    if(k < 2){
        stop("ERROR: k must be > 1!")
    }else if(nrow(elist$E) < candidate.number){
        stop(paste("ERROR: number of features (here: ", nrow(elist$E),
          ") must be >= candidate.number (here: ", candidate.number, ")!",
          sep=""))
    }
    
    
    
    if(is.null(log)) {
        stop("selectFeatures - ERROR: Not all mandatory arguments defined!")
    }
    
    if(preselection.method=="mMs" && log==TRUE){
        elist$E <- 2^(elist$E)
    }else if(preselection.method=="tTest" && log==FALSE){
        elist$E <- log2(elist$E)
    }else if(preselection.method=="mrmr" && log==FALSE){
        elist$E <- log2(elist$E)
    }
    
    
    
    if(method=="ensemble"){
        selectFeatures.ensemble(elist=elist, n1=n1, n2=n2, label1=label1,
          label2=label2, cutoff=cutoff, k=k, subsamples=subsamples,
          bootstraps=bootstraps, plot=plot, output.path=output.path,
          verbose=verbose)    
    }else if(method=="frequency"){
        selectFeatures.frequency(elist=elist, n1=n1, n2=n2, label1=label1,
          label2=label2, cutoff=cutoff, selection.method=selection.method,
          preselection.method=preselection.method, subruns=subruns, k=k,
          candidate.number=candidate.number, above=above, between=between,
          panel.selection.criterion=panel.selection.criterion,
          importance.measure=importance.measure, ntree=ntree, mtry=mtry,
          plot=plot, output.path=output.path, verbose=verbose) 
    }else{
        stop("selectFeatures(): Error - Unknown feature selection method!!!")
    }
}

#+++++++++++++++++++++++++ selectFeatures ++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++ plotFeatures ++++++++++++++++++++++++++++++++++++++
plotFeatures <- function(features=NULL, elist=NULL, n1=NULL, n2=NULL,
  group1="group1", group2="group2", output.path=NULL){
    
    if(is.null(features) || is.null(elist) || is.null(n1) || is.null(n2)) {
        stop("ERROR: Not all mandatory arguments have been defined!")
    }
    datamatrix <- cbind(paste(elist$genes$Block,elist$genes$Row,
      elist$genes$Column),elist$E)
    datamatrix <- datamatrix[datamatrix[,1] %in% features,]
    len <- length(datamatrix[,1])
    if(!is.null(output.path)){
        tiff(paste(output.path,"/biomarker_plots.tiff",sep=""), width=2500,
          height=500*ceiling({len+1}/5), pointsize=6, compression="lzw",
          res=300) 
        par(mfrow=c(ceiling({len+1}/5),5), cex.main=3, cex.lab=1, cex.axis=2.5,
          font.axis=2)
            for (i in 1:{ceiling((len+1)/5)*5}){
                if(i <= len){
                    plot(as.numeric(datamatrix[i,2:{n1+n2+1}]), col="red",
                      xlim=c(1,max(n1,n2)), xlab="", ylab="", pch=16,
                      main=datamatrix[i,1], axes=FALSE, type="n")
                    axis(2)
                    box()
                    points(as.numeric(datamatrix[i,2:{n1+1}]), col="red",
                      pch=16, cex=2.5)
                    points(as.numeric(datamatrix[i,{n1+2}:{1+n1+n2}]),
                      col="blue", pch=17, cex=2.5)
                }else{
                    plot(x=1,xlab="",ylab="",axes=FALSE,type="n")
                }
            }
            legend("topleft", legend=c(group1, group2), col=c("red", "blue"),
              pch=c(16,17), cex=5)
        dev.off()
    } else {
        #south, west, north, east
        par(mfrow=c(ceiling({len+1}/5),5), cex.main=1.5, cex.lab=1, cex.axis=1,
          font.axis=2, mar=c(2,2,1.5,0.5), oma=c(1.5,2,1,1))
            for (i in 1:{ceiling((len+1)/5)*5}){
                if(i <= len){
                    plot(as.numeric(datamatrix[i,2:{n1+n2+1}]), col="red",
                      xlim=c(1,max(n1,n2)), xlab="", ylab="", pch=16,
                      main=datamatrix[i,1], axes=FALSE, type="n")
                    axis(2)
                    box()
                    points(as.numeric(datamatrix[i,2:{n1+1}]), col="red",
                      pch=16, cex=2)
                    points(as.numeric(datamatrix[i,{n1+2}:{1+n1+n2}]),
                      col="blue", pch=17, cex=2)
                }else{
                    plot(x=1,xlab="",ylab="",axes=FALSE,type="n")
                }
            }
        legend("topleft", legend=c(group1, group2), col=c("red", "blue"),
          pch=c(16,17), cex=1.5)
    }
}
#++++++++++++++++++++++++++ plotFeatures +++++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++ plotFeaturesHeatmap ++++++++++++++++++++++++++++++++++++
plotFeaturesHeatmap <- function(features=NULL, elist=NULL, n1=NULL, n2=NULL,
  output.path=NULL, description=FALSE){
    
    if(is.null(features) || is.null(elist) || is.null(n1) || is.null(n2)) {
        stop("ERROR: Not all mandatory arguments have been defined!")
    }
    features <- c(na.exclude(features))
    datamatrix <- elist$E
    brc <- paste(elist$genes$Block,elist$genes$Row,elist$genes$Column)
    #datamatrix <- datamatrix[rownames(datamatrix) %in% features,]
    rows <- brc %in% features
    datamatrix <- datamatrix[rows,]
    if(description){
        rownames(datamatrix) <- elist$genes$Description[rows]
        mar2 <- 15
        cxRow=0.5
    }else{
        rownames(datamatrix) <- brc[rows]
        mar2 <- 5
        cxRow=1
    }
    #my.dist <- function(x) as.dist(1 - abs(cor(t(x), method="pearson")))
    my.dist <- function(x) as.dist((1-cor(t(x)))/2)
      #"euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
    my.hclust <- function(d) hclust(d, method="complete")
      #"ward", "single", "complete", "average", "mcquitty",
      #"median" or "centroid"
    #if(!is.null(output.path)){
    #    tiff(paste(output.path,"/biomarker_heatmap.tiff",sep=""), width=2000,
    #      height=2000, pointsize=15, compression="lzw", res=300) 
    #        heatmap(datamatrix,Rowv=NA,Colv=NA,distfun=my.dist,
    #          hclustfun=my.hclust,col = cm.colors(300), cexCol=0.5,
    #          cexRow=cxRow, margins=c(2, mar2))
    #          #col:cm.colors,terrain.colors,heat.colors,topo.colors,rainbow,gray
    #    dev.off()
    #} else {
    #    heatmap(datamatrix,Rowv=NA,Colv=NA,distfun=my.dist,
    #      hclustfun=my.hclust,col = cm.colors(300), cexCol=0.5, cexRow=cxRow,
    #      margins=c(2, mar2))
    #      #col:cm.colors,terrain.colors,heat.colors,topo.colors,rainbow,gray
    #}
    
    if(!is.null(output.path)){
        tiff(paste(output.path,"/biomarker_heatmap.tiff",sep=""), width=2000,
          height=2000, pointsize=15, compression="lzw", res=300) 
            heatmap(datamatrix,Colv=NA,keep.dendro=FALSE,distfun=my.dist,
              hclustfun=my.hclust,scale="row",
              #col = bluered(300),
              col = blueyellow(300),
              cexCol=0.5,cexRow=cxRow, margins=c(2, mar2))
              #col:cm.colors,terrain.colors,heat.colors,topo.colors,rainbow,gray
        dev.off()
    } else {
        heatmap(datamatrix,Colv=NA,keep.dendro=FALSE,distfun=my.dist,
          hclustfun=my.hclust,scale="row",
          #col = bluered(300),
          col = blueyellow(300),
          cexCol=0.5, cexRow=cxRow, margins=c(2, mar2))
          #col:cm.colors,terrain.colors,heat.colors,topo.colors,rainbow,gray
    }
}
#++++++++++++++++++++++ plotFeaturesHeatmap ++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++ printFeatures +++++++++++++++++++++++++++++++++++++++
printFeatures <- function(features=NULL, elist=NULL, output.path=NULL){
    if(is.null(features) || is.null(elist)) {
        stop("ERROR: Not all mandatory arguments have been defined!")
    }
    datamatrix <- cbind(paste(elist$genes$Block,elist$genes$Row,
      elist$genes$Column), elist$genes$Description,elist$E)
    datamatrix <- datamatrix[datamatrix[,1] %in% features,]
    colnames(datamatrix) <- c("BRC","Description", colnames(elist$E))
    if(is.null(output.path)){
        return(data.frame(datamatrix))
    } else {
        write.table(x=datamatrix,
          file=paste(output.path, "/candidates.txt", sep=""),
          sep="\t", eol="\n", row.names=FALSE, quote=FALSE)
    }
}
#+++++++++++++++++++++++++ printFeatures +++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++ plotMAPlots +++++++++++++++++++++++++++++++++++++++
plotMAPlots <- function(elist=NULL, idx="all", include.rlm=FALSE,
  controls="internal", output.path=NULL){
    
    if(is.null(elist) || is.null(idx)) {
        stop("ERROR: Not all mandatory arguments have been defined!")
    }
    if(!is.null(output.path)){
        dir.create(paste(output.path, "/MA_plots", sep=""))
    }
    
    col.number <- ncol(elist$E)
    
    if(idx == "all"){
        if(is.null(output.path)){
            stop(paste("plotMAPlots - ERROR: an output.path has to be",
              "defined when argument 'idx' is set to 'all'!!!", sep=""))
        }
    }else if(is.numeric(idx)){
        idx <- round(abs(idx))
        if(idx < 1 || idx > col.number){
            stop(paste("plotMAPlots - ERROR: 'idx' has to be a positive",
              "integer that is not larger than the number of columns",
              "(or the string 'all')!", sep=""))
        }
    }else{
        stop(paste("plotMAPlots - ERROR: argument 'idx' has to be an",
          "integer or the string 'all'!", sep=""))
    }
    elist.quantile <- normalizeArrays(elist=elist, method="quantile")
    elist.cyclicloess.pairs <- normalizeArrays(elist=elist,
      method="cyclicloess", cyclicloess.method="pairs")
      #--> cyclicloess.method = ('pairs'|'fast'|'affy')
    elist.vsn <- normalizeArrays(elist=elist, method="vsn")
    if(include.rlm){
        elist.rlm <- normalizeArrays(elist=elist, method="rlm",
          controls=controls, output.path=NULL)
    }
    
    if(idx == "all"){
        for(i in 1:col.number) {
            A_raw <- apply(log2(elist$E),MARGIN=1,FUN=median)
            M_raw <- log2(elist$E[,i]) - A_raw      
            A_vsn <- apply(elist.vsn$E,MARGIN=1,FUN=median)
            M_vsn <- elist.vsn$E[,i] - A_vsn       
            A_quant <- apply(elist.quantile$E,MARGIN=1,FUN=median)
            M_quant <- elist.quantile$E[,i] - A_quant       
            A_cyclic_pairs <-
              apply(elist.cyclicloess.pairs$E,MARGIN=1,FUN=median)
            M_cyclic_pairs <- elist.cyclicloess.pairs$E[,i] - A_cyclic_pairs        
            if(include.rlm){
                A_rlm <- apply(elist.rlm$E,MARGIN=1,FUN=median)
                M_rlm <- elist.rlm$E[,i] - A_rlm
              
                tiff(paste(output.path,"/MA_plots/maplots_", colnames(elist)[i],
                  ".tiff",sep=""), width=6000, height=4000, pointsize=20,
                  compression="lzw", res=300)
                par(mfrow=c(2,3))   # --> c(nr, nc)
                    plot(A_raw,M_raw, main=paste(colnames(elist)[i], " raw",
                      sep=""), xlab="A", ylab="M", cex=0.5, pch=1)
                    abline(0,0,col="blue",lwd=5)
                    lines(lowess(M_raw~A_raw),col="red",lwd=5)
                
                    plot(A_vsn,M_vsn, main=paste(colnames(elist)[i], " vsn",
                      sep=""), xlab="A", ylab="M", cex=0.5, pch=1)
                    abline(0,0,col="blue",lwd=5)
                    lines(lowess(M_vsn~A_vsn),col="red",lwd=5)
                
                    plot(A_quant,M_quant, main=paste(colnames(elist)[i],
                      " quantile", sep=""), xlab="A", ylab="M", cex=0.5, pch=1)
                    abline(0,0,col="blue",lwd=5)
                    lines(lowess(M_quant~A_quant),col="red",lwd=5)
                
                    plot(A_cyclic_pairs,M_cyclic_pairs,
                      main=paste(colnames(elist)[i], " cyclic loess", sep=""),
                      xlab="A", ylab="M", cex=0.5, pch=1)
                    abline(0,0,col="blue",lwd=5)
                    lines(lowess(M_cyclic_pairs~A_cyclic_pairs),col="red",lwd=5)
                
                    plot(A_rlm,M_rlm, main=paste(colnames(elist)[i], " rlm",
                      sep=""), xlab="A", ylab="M", cex=0.5, pch=1)
                    abline(0,0,col="blue",lwd=5)
                    lines(lowess(M_rlm~A_rlm),col="red",lwd=5)
                dev.off()
            }else{
                tiff(paste(output.path,"/MA_plots/maplots_", colnames(elist)[i],
                  ".tiff",sep=""), width=4000, height=4000, pointsize=20,
                  compression="lzw", res=300)
                par(mfrow=c(2,2))   # --> c(nr, nc)
                    plot(A_raw,M_raw, main=paste(colnames(elist)[i],
                      " raw", sep=""), xlab="A", ylab="M", cex=0.5, pch=1)
                    abline(0,0,col="blue",lwd=5)
                    lines(lowess(M_raw~A_raw),col="red",lwd=5)
                
                    plot(A_vsn,M_vsn, main=paste(colnames(elist)[i],
                      " vsn", sep=""), xlab="A", ylab="M", cex=0.5, pch=1)
                    abline(0,0,col="blue",lwd=5)
                    lines(lowess(M_vsn~A_vsn),col="red",lwd=5)
                
                    plot(A_quant,M_quant, main=paste(colnames(elist)[i],
                      " quantile", sep=""), xlab="A", ylab="M", cex=0.5, pch=1)
                    abline(0,0,col="blue",lwd=5)
                    lines(lowess(M_quant~A_quant),col="red",lwd=5)
                
                    plot(A_cyclic_pairs,M_cyclic_pairs,
                      main=paste(colnames(elist)[i], " cyclic loess", sep=""),
                      xlab="A", ylab="M", cex=0.5, pch=1)
                    abline(0,0,col="blue",lwd=5)
                    lines(lowess(M_cyclic_pairs~A_cyclic_pairs),col="red",lwd=5)
                dev.off()
            }
            message('drawing MA plots: ', round(i*100/{col.number-1}),
              '%\r', sep="")
        }
        message('drawing MA plots: ', '100', '%\n', sep="")
    }else if(is.numeric(idx)){
        A_raw <- apply(log2(elist$E),MARGIN=1,FUN=median)
        M_raw <- log2(elist$E[,idx]) - A_raw      
        A_vsn <- apply(elist.vsn$E,MARGIN=1,FUN=median)
        M_vsn <- elist.vsn$E[,idx] - A_vsn       
        A_quant <- apply(elist.quantile$E,MARGIN=1,FUN=median)
        M_quant <- elist.quantile$E[,idx] - A_quant       
        A_cyclic_pairs <- apply(elist.cyclicloess.pairs$E,MARGIN=1,FUN=median)
        M_cyclic_pairs <- elist.cyclicloess.pairs$E[,idx] - A_cyclic_pairs
        if(include.rlm){
            A_rlm <- apply(elist.rlm$E,MARGIN=1,FUN=median)
            M_rlm <- elist.rlm$E[,idx] - A_rlm
            if(is.null(output.path)){
                par(mfrow=c(2,3))   # --> c(nr, nc)
                    plot(A_raw,M_raw, main=paste(colnames(elist)[idx], " raw",
                      sep=""), xlab="A", ylab="M")
                    abline(0,0,col="blue",lwd=5)
                    lines(lowess(M_raw~A_raw),col="red",lwd=5)
                
                    plot(A_vsn,M_vsn, main=paste(colnames(elist)[idx], " vsn",
                      sep=""), xlab="A", ylab="M")
                    abline(0,0,col="blue",lwd=5)
                    lines(lowess(M_vsn~A_vsn),col="red",lwd=5)
                
                    plot(A_quant,M_quant, main=paste(colnames(elist)[idx],
                      " quantile", sep=""), xlab="A", ylab="M")
                    abline(0,0,col="blue",lwd=5)
                    lines(lowess(M_quant~A_quant),col="red",lwd=5)
                
                    plot(A_cyclic_pairs,M_cyclic_pairs,
                      main=paste(colnames(elist)[idx], " cyclic loess", sep=""),
                      xlab="A", ylab="M")
                    abline(0,0,col="blue",lwd=5)
                    lines(lowess(M_cyclic_pairs~A_cyclic_pairs),col="red",lwd=5)
                
                    plot(A_rlm,M_rlm, main=paste(colnames(elist)[idx], " rlm",
                      sep=""), xlab="A", ylab="M")
                    abline(0,0,col="blue",lwd=5)
                    lines(lowess(M_rlm~A_rlm),col="red",lwd=5)
            }else{
                tiff(paste(output.path,"/MA_plots/maplots_",
                  colnames(elist)[idx], ".tiff",sep=""), width=6000,
                  height=4000, pointsize=20, compression="lzw", res=300)
                par(mfrow=c(2,3))   # --> c(nr, nc)
                    plot(A_raw,M_raw, main=paste(colnames(elist)[idx], " raw",
                      sep=""), xlab="A", ylab="M", cex=0.5, pch=1)
                    abline(0,0,col="blue",lwd=5)
                    lines(lowess(M_raw~A_raw),col="red",lwd=5)
                
                    plot(A_vsn,M_vsn, main=paste(colnames(elist)[idx],
                      " vsn", sep=""), xlab="A", ylab="M", cex=0.5, pch=1)
                    abline(0,0,col="blue",lwd=5)
                    lines(lowess(M_vsn~A_vsn),col="red",lwd=5)
                
                    plot(A_quant,M_quant, main=paste(colnames(elist)[idx],
                      " quantile", sep=""), xlab="A", ylab="M", cex=0.5, pch=1)
                    abline(0,0,col="blue",lwd=5)
                    lines(lowess(M_quant~A_quant),col="red",lwd=5)
                
                    plot(A_cyclic_pairs,M_cyclic_pairs,
                      main=paste(colnames(elist)[idx], " cyclic loess", sep=""),
                      xlab="A", ylab="M", cex=0.5, pch=1)
                    abline(0,0,col="blue",lwd=5)
                    lines(lowess(M_cyclic_pairs~A_cyclic_pairs),col="red",lwd=5)
                
                    plot(A_rlm,M_rlm, main=paste(colnames(elist)[idx], " rlm",
                      sep=""), xlab="A", ylab="M", cex=0.5, pch=1)
                    abline(0,0,col="blue",lwd=5)
                    lines(lowess(M_rlm~A_rlm),col="red",lwd=5)
                dev.off()
            }
        }else{
            if(is.null(output.path)){
                par(mfrow=c(2,2))   # --> c(nr, nc)
                    plot(A_raw,M_raw, main=paste(colnames(elist)[idx], " raw",
                      sep=""), xlab="A", ylab="M")
                    abline(0,0,col="blue",lwd=5)
                    lines(lowess(M_raw~A_raw),col="red",lwd=5)
                
                    plot(A_vsn,M_vsn, main=paste(colnames(elist)[idx], " vsn",
                      sep=""), xlab="A", ylab="M")
                    abline(0,0,col="blue",lwd=5)
                    lines(lowess(M_vsn~A_vsn),col="red",lwd=5)
                
                    plot(A_quant,M_quant, main=paste(colnames(elist)[idx],
                      " quantile", sep=""), xlab="A", ylab="M")
                    abline(0,0,col="blue",lwd=5)
                    lines(lowess(M_quant~A_quant),col="red",lwd=5)
                
                    plot(A_cyclic_pairs,M_cyclic_pairs,
                      main=paste(colnames(elist)[idx], " cyclic loess",
                      sep=""), xlab="A", ylab="M")
                    abline(0,0,col="blue",lwd=5)
                    lines(lowess(M_cyclic_pairs~A_cyclic_pairs),col="red",lwd=5)
            }else{
                tiff(paste(output.path,"/MA_plots/maplots_",
                  colnames(elist)[idx], ".tiff",sep=""), width=4000,
                  height=4000, pointsize=20, compression="lzw", res=300)
                par(mfrow=c(2,2))   # --> c(nr, nc)
                    plot(A_raw,M_raw, main=paste(colnames(elist)[idx], " raw",
                      sep=""), xlab="A", ylab="M", cex=0.5, pch=1)
                    abline(0,0,col="blue",lwd=5)
                    lines(lowess(M_raw~A_raw),col="red",lwd=5)
                
                    plot(A_vsn,M_vsn, main=paste(colnames(elist)[idx], " vsn",
                      sep=""), xlab="A", ylab="M", cex=0.5, pch=1)
                    abline(0,0,col="blue",lwd=5)
                    lines(lowess(M_vsn~A_vsn),col="red",lwd=5)
                
                    plot(A_quant,M_quant, main=paste(colnames(elist)[idx],
                      " quantile", sep=""), xlab="A", ylab="M", cex=0.5, pch=1)
                    abline(0,0,col="blue",lwd=5)
                    lines(lowess(M_quant~A_quant),col="red",lwd=5)
                
                    plot(A_cyclic_pairs,M_cyclic_pairs,
                      main=paste(colnames(elist)[idx], " cyclic loess", sep=""),
                      xlab="A", ylab="M", cex=0.5, pch=1)
                    abline(0,0,col="blue",lwd=5)
                    lines(lowess(M_cyclic_pairs~A_cyclic_pairs),col="red",lwd=5)
                dev.off()
            }
        }
    }
}
#+++++++++++++++++++++++++++ plotMAPlots +++++++++++++++++++++++++++++++++++++++  

#+++++++++++++++++++++++++++ shuffleData +++++++++++++++++++++++++++++++++++++++

shuffleData <- function(elist=NULL, n1=NULL, n2=NULL, label1="A", label2="B"){
    if(is.null(elist) || is.null(n1) || is.null(n2)) {
        stop("ERROR: Not all mandatory arguments have been defined!")
    }
    array.id <-
      c(paste(rep(label1, n1), 1:n1, sep=""), paste(rep(label2, n2), 1:n2,
      sep=""))
    array.id <- sample(array.id)
    colnames(elist$E) <- array.id
    elist$targets$SerumID <- array.id
    elist$targets$ArrayID <- array.id
    elist$targets$Group <- gsub("\\d+", "", array.id)
    elist$E <-
      elist$E[ ,order(elist$targets$Group, as.numeric(sub("\\D*", "",
      elist$targets$ArrayID)))] #Spalten wegen 'order' zuerst sortieren!
    elist$targets <- elist$targets[order(elist$targets$Group,
      as.numeric(sub("\\D*", "", elist$targets$ArrayID))), ]  #Zeilen
    return(elist)
}

#+++++++++++++++++++++++++++ shuffleData +++++++++++++++++++++++++++++++++++++++        

#+++++++++++++++++++++++++++ pvaluePlot +++++++++++++++++++++++++++++++++++++++
pvaluePlot <- function(elist=NULL, group1=NULL, group2=NULL, log=NULL,
  method="tTest", output.path=NULL, tag="", mMs.matrix1=NULL, mMs.matrix2=NULL,
  above=1500, between=400, adjust=FALSE){
    
    if(is.null(elist) || is.null(group1) || is.null(group2) || is.null(log)) {
        stop("ERROR: Not all mandatory arguments have been defined!")
    }
    
    if(method=="mMs" && log==TRUE){
        elist$E <- 2^(elist$E)
    }
    
    row.number <- nrow(elist$E)
    col.number <- ncol(elist$E)
    pvalues <- vector(mode="numeric", length=row.number)

    if(method=="mMs"){
        mCount.results1 <-
          mCount(x=elist$E[,group1], y=elist$E[,group2], z=mMs.matrix1,
          a=above, b=between)
        mCount.results2 <-
          mCount(x=elist$E[,group2], y=elist$E[,group1], z=mMs.matrix2,
          a=above, b=between)
    }
    for(zeile in 1:row.number) {
        #fold <- mean(elist$E[zeile, group1])/mean(elist$E[zeile, group2])
        if(method=="tTest"){
            p <- try(t.test(x=elist$E[zeile,group1],
              y=elist$E[zeile,group2])$p.value, TRUE)
            if(!is.numeric(p)){p <- 1}
        }else if(method=="mMs"){
            min.p1 <- mCount.results1$mMs[zeile]
            min.p2 <- mCount.results2$mMs[zeile]
            p <- min(min.p1,min.p2)
        }
        pvalues[zeile] <- p
    }
    if(adjust){
        pvalues <- p.adjust(p=pvalues, method="fdr")
    }

    number1 <- length(pvalues[pvalues < 0.05])
    number2 <- length(pvalues[pvalues < 0.01])
    number3 <- length(pvalues[pvalues < 0.05/row.number])
    ylimit <- min(log2(0.04/row.number), log2(min(pvalues))-1)

    if(is.null(output.path)) {
        if(adjust){
            plot(sort(pvalues), ylab="fdr", main=paste("FDRs\n(", number1,
              " < 0.05, ", number2, " < 0.01)", sep=""))
            abline(h=0.01, lty=2, col="red")
            abline(h=0.05, lty=2, col="red")
        }else{
            plot(log2(sort(pvalues)), ylim=c(ylimit,0), ylab="log2(p-value)",
              main=paste("p-values\n(", number1, " < 0.05, ", number2,
              " < 0.01, ", number3, " < 0.05 after Bonferroni)", sep=""))
            abline(h=log2(0.05), lty=2, col="red")
            abline(h=log2(0.01), lty=2, col="green")
            abline(h=log2(0.05/row.number), lty=2, col="blue")
            legend("bottomright", legend=c("p-value = 0.05", "p-value = 0.01",
              "Bonferroni"), col=c("red", "green", "blue"), lty=c(2,2,2))
        }
    } else {
        if(adjust){
            tiff(paste(output.path,"/pvalues_fdr", tag, ".tiff",sep=""),
              width=2000, height=2000, pointsize=15, compression="lzw", res=300)
                plot(sort(pvalues), ylab="fdr", main=paste("FDRs\n(", number1,
                  " < 0.05, ", number2, " < 0.01)", sep=""))
                abline(h=0.01, lty=2, col="red")
                abline(h=0.05, lty=2, col="red")
            dev.off()
        }else{
            tiff(paste(output.path,"/pvalues", tag, ".tiff",sep=""), width=2000,
              height=2000, pointsize=15, compression="lzw", res=300)
                plot(log2(sort(pvalues)), ylim=c(ylimit,0),
                  ylab="log2(p-value)", main=paste("p-values\n(", number1,
                  " < 0.05, ", number2, " < 0.01, ", number3,
                  " < 0.05 after Bonferroni)", sep=""))
                abline(h=log2(0.05), lty=2, col="red")
                abline(h=log2(0.01), lty=2, col="green")
                abline(h=log2(0.05/row.number), lty=2, col="blue")
                legend("bottomright", legend=c("p-value = 0.05",
                  "p-value = 0.01", "Bonferroni"), col=c("red", "green",
                  "blue"), lty=c(2,2,2))
            dev.off()
        }
    }
}
#+++++++++++++++++++++++++++ pvaluePlot ++++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++ volcanoPlot +++++++++++++++++++++++++++++++++++++++
volcanoPlot <- function(elist=NULL, group1=NULL, group2=NULL, log=NULL,
  method="tTest", p.thresh=NULL, fold.thresh=NULL, output.path=NULL, tag="",
  mMs.matrix1=NULL, mMs.matrix2=NULL, above=1500, between=400){
    
    if(is.null(elist) || is.null(group1) || is.null(group2) || is.null(log)) {
        stop("ERROR: Not all mandatory arguments have been defined!")
    }
    
    if(method=="tTest" && log==FALSE){
        elist$E <- log2(elist$E)
    }else if(method=="mMs" && log==TRUE){
        elist$E <- 2^(elist$E)
    }
    
    row.number <- nrow(elist$E)
    col.number <- ncol(elist$E)
    mark <- c()
    output <- matrix(nrow=row.number, ncol={4+ncol(elist$E)})
    colnames(output) <- c("BRC", "P-value", "Fold change", "Description",
      colnames(elist$E))

    if(method=="mMs"){
        mCount.results1 <-
          mCount(x=elist$E[,group1], y=elist$E[,group2], z=mMs.matrix1,
          a=above, b=between)
        mCount.results2 <-
          mCount(x=elist$E[,group2], y=elist$E[,group1], z=mMs.matrix2,
          a=above, b=between)
    }
    for(zeile in 1:row.number) {
        if(log == TRUE){
            fold <- 
              2^(mean(elist$E[zeile, group1]) - mean(elist$E[zeile, group2]))
        }else{
            fold <- mean(elist$E[zeile, group1]) / mean(elist$E[zeile, group2])
        }
        if(method=="tTest"){
            p <- try(t.test(x=elist$E[zeile,group1],
              y=elist$E[zeile,group2])$p.value, TRUE)
            if(!is.numeric(p)){p <- 1}
        }else if(method=="mMs"){
            min.p1 <- mCount.results1$mMs[zeile]
            min.p2 <- mCount.results2$mMs[zeile]
            p <- min(min.p1,min.p2)
        }
        if(!is.null(p.thresh) && !is.null(fold.thresh)){   # --> A & B
            if(p < p.thresh || fold > fold.thresh || fold < 1/fold.thresh) {
            #---> ISSUE: there  may be no exact representation for floats!!!
            #if((round(p,3) < round(p.thresh,3)) || fold > fold.thresh ||
            #  fold < 1/fold.thresh) {
            #  ---> ISSUE: there  may be no exact representation for floats!!!
                mark <- c(mark,zeile)
            }
        }else if(!is.null(p.thresh) && is.null(fold.thresh)) {    # --> A & -B
            if(p < p.thresh) {
            #---> ISSUE: there  may be no exact representation for floats!!!
                mark <- c(mark,zeile)
            }
        }else if(is.null(p.thresh) && !is.null(fold.thresh)) {    # --> -A & B
            if(fold > fold.thresh || fold < 1/fold.thresh) {
            #---> ISSUE: there  may be no exact representation for floats!!!
                mark <- c(mark,zeile)
            }
        }
        output[zeile,] <- c(paste(elist$genes$Block[zeile],
          elist$genes$Row[zeile], elist$genes$Column[zeile]), p, fold,
          elist$genes$Description[zeile], log2(elist$E[zeile,]))
    }
    #write.table(x=output, file=paste(output.path, "/volcano.txt", sep=""),
    #  sep="\t", eol="\n", row.names=FALSE, quote=FALSE)
    if(is.null(output.path) && is.null(p.thresh) && is.null(fold.thresh)){
      # --> -A & -B & -C
        plot(x=log2(as.numeric(output[,3])), y=-log10(as.numeric(output[,2])),
          xlab="log2(fold change)", ylab="-log10(p-value)", main="Volcano Plot")
    }else if(!is.null(output.path) && is.null(p.thresh) &&
      is.null(fold.thresh)) {    # --> A & -B & -C
        tiff(paste(output.path,"/volcano", tag, ".tiff",sep=""), width=2000,
          height=2000, pointsize=15, compression="lzw", res=300)
            plot(x=log2(as.numeric(output[,3])),
              y=-log10(as.numeric(output[,2])), xlab="log2(fold change)",
              ylab="-log10(p-value)", main="Volcano Plot")
        dev.off()
    }else if(is.null(output.path) && !is.null(p.thresh) &&
      is.null(fold.thresh)) { # --> -A & B & -C
        plot(x=log2(as.numeric(output[,3])), y=-log10(as.numeric(output[,2])),
          xlab="log2(fold change)", ylab="-log10(p-value)", main="Volcano Plot")
        points(x=log2(as.numeric(output[mark,3])),
          y=-log10(as.numeric(output[mark,2])), col="red")
        abline(h=-log10(p.thresh), lty=2, col="red") 
    }else if(!is.null(output.path) && !is.null(p.thresh) &&
      is.null(fold.thresh)) {   # --> A & B & -C
        tiff(paste(output.path,"/volcano", tag, ".tiff",sep=""), width=2000,
          height=2000, pointsize=15, compression="lzw", res=300)
        plot(x=log2(as.numeric(output[,3])), y=-log10(as.numeric(output[,2])),
          xlab="log2(fold change)", ylab="-log10(p-value)", main="Volcano Plot")
            points(x=log2(as.numeric(output[mark,3])),
              y=-log10(as.numeric(output[mark,2])), col="red")
            abline(h=-log10(p.thresh), lty=2, col="red")
        dev.off()  
    }else if(is.null(output.path) && is.null(p.thresh) &&
      !is.null(fold.thresh)) {   # --> -A & -B & C
        plot(x=log2(as.numeric(output[,3])), y=-log10(as.numeric(output[,2])),
          xlab="log2(fold change)", ylab="-log10(p-value)", main="Volcano Plot")
        points(x=log2(as.numeric(output[mark,3])),
          y=-log10(as.numeric(output[mark,2])), col="red")
        abline(v=log2(fold.thresh), lty=2, col="red")
        abline(v=log2(1/fold.thresh), lty=2, col="red")
    }else if(!is.null(output.path) && is.null(p.thresh) &&
      !is.null(fold.thresh)) {    # --> A & -B & C
        tiff(paste(output.path,"/volcano", tag, ".tiff",sep=""), width=2000,
          height=2000, pointsize=15, compression="lzw", res=300)
        plot(x=log2(as.numeric(output[,3])), y=-log10(as.numeric(output[,2])),
          xlab="log2(fold change)", ylab="-log10(p-value)", main="Volcano Plot")
            points(x=log2(as.numeric(output[mark,3])),
              y=-log10(as.numeric(output[mark,2])), col="red")
            abline(v=log2(fold.thresh), lty=2, col="red")
            abline(v=log2(1/fold.thresh), lty=2, col="red")
        dev.off()
    }else if(is.null(output.path) && !is.null(p.thresh) &&
      !is.null(fold.thresh)) {    # --> -A & B & C
        plot(x=log2(as.numeric(output[,3])), y=-log10(as.numeric(output[,2])),
          xlab="log2(fold change)", ylab="-log10(p-value)", main="Volcano Plot")
        points(x=log2(as.numeric(output[mark,3])),
          y=-log10(as.numeric(output[mark,2])), col="red")
        abline(h=-log10(p.thresh), lty=2, col="red")
        abline(v=log2(fold.thresh), lty=2, col="red")
        abline(v=log2(1/fold.thresh), lty=2, col="red")
    }else if(!is.null(output.path) && !is.null(p.thresh) &&
      !is.null(fold.thresh)) {    # --> A & B & C
        tiff(paste(output.path,"/volcano", tag, ".tiff",sep=""), width=2000,
          height=2000, pointsize=15, compression="lzw", res=300)
            plot(x=log2(as.numeric(output[,3])),
              y=-log10(as.numeric(output[,2])), xlab="log2(fold change)",
              ylab="-log10(p-value)", main="Volcano Plot")
            points(x=log2(as.numeric(output[mark,3])),
              y=-log10(as.numeric(output[mark,2])), col="red")
            abline(h=-log10(p.thresh), lty=2, col="red")
            abline(v=log2(fold.thresh), lty=2, col="red")
            abline(v=log2(1/fold.thresh), lty=2, col="red")
        dev.off()
    }
}
#+++++++++++++++++++++++++++ volcanoPlot +++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++ diffAnalysis ++++++++++++++++++++++++++++++++++++++
diffAnalysis <- function(input=NULL, label1=NULL, label2=NULL, class1=NULL,
  class2=NULL, output.path=NULL, mMs.matrix1=NULL, mMs.matrix2=NULL,
  above=1500, between=400, features=NULL, feature.names=NULL){
    
    if(is.null(input) || is.null(label1) || is.null(label2) || is.null(class1)
      || is.null(class2) || is.null(mMs.matrix1) || is.null(mMs.matrix2)) {
        stop("ERROR: Not all mandatory arguments have been defined!")
    }
    
    row.number <- nrow(input)
    col.number <- ncol(input)
    output <- matrix(nrow=row.number, ncol=12)
    colnames(output) <- c("BRC", "t test", "FDR(t)", "min. M stat. (mMs)",
      "FDR(mMs)", "fold change", paste("mean ", class1, sep=""),
      paste("mean ", class2, sep=""), paste("median ", class1, sep=""),
      paste("median ", class2, sep=""), paste("sd ", class1, sep=""),
      paste("sd ", class2, sep=""))

    mCount.results1 <-
      mCount(x=input[,label1], y=input[,label2], z=mMs.matrix1, a=above,
      b=between)
    mCount.results2 <-
      mCount(x=input[,label2], y=input[,label1], z=mMs.matrix2, a=above,
      b=between)
    for(zeile in 1:row.number) {
        output[zeile,1] <- rownames(input)[zeile]

        p <- try(t.test(x=as.numeric(input[zeile,label1]),
          y=as.numeric(input[zeile,label2]))$p.value, TRUE)
        if(!is.numeric(p)){output[zeile,2] <- 1}else{output[zeile,2] <- p}
        
        min.p1 <- mCount.results1$mMs[zeile]
        m.order1 <- mCount.results1$order[zeile]
        m.count1 <- mCount.results1$count[zeile]
        m.cutoff1 <- mCount.results1$cutoff[zeile]
        mean2 <- mCount.results1$means2[zeile]
        min.p2 <- mCount.results2$mMs[zeile]
        m.order2 <- mCount.results2$order[zeile]
        m.count2 <- mCount.results2$count[zeile]
        m.cutoff2 <- mCount.results2$cutoff[zeile]
        mean1 <- mCount.results2$means2[zeile]
        output[zeile,4] <- min(min.p1,min.p2) 

        output[zeile,6] <- mean1/mean2
        output[zeile,7] <- mean1
        output[zeile,8] <- mean2
        output[zeile,9] <- median(as.numeric(input[zeile, label1]))
        output[zeile,10] <- median(as.numeric(input[zeile, label2]))
        output[zeile,11] <- sd(as.numeric(input[zeile, label1]))       
        output[zeile,12] <- sd(as.numeric(input[zeile, label2]))
    }
    output[,3] <- p.adjust(p=as.numeric(output[,2]), method="fdr")
    output[,5] <- p.adjust(p=as.numeric(output[,4]), method="fdr")
    
    if(!is.null(output.path)){
        if(is.null(features)){
            write.table(x=output, file=paste(output.path, "/diffAnalysis.txt",
              sep=""), sep="\t", eol="\n", row.names=FALSE, quote=FALSE)
        }else{
            if(is.null(feature.names)){
                write.table(x=output[features,], file=paste(output.path,
                  "/diffAnalysis.txt", sep=""), sep="\t", eol="\n",
                  row.names=FALSE, quote=FALSE)
            }else{
                write.table(x=cbind(output[features,],feature.names),
                  file=paste(output.path, "/diffAnalysis.txt", sep=""),
                  sep="\t", eol="\n", row.names=FALSE, quote=FALSE)
            }
        }
    }
    return(data.frame(output))
}
#++++++++++++++++++++++++++ diffAnalysis +++++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++++++++ plotArrayTiff ++++++++++++++++++++++++++++++++++++
plotArrayTiff <- function (plot.ma=NULL, tag=NULL, idx=NULL, array.label=NULL,
  block.dim=NULL, colpal="heat.colors", output.path=NULL) {
    
  if(!is.null(output.path)){
      #Create output directory if not present. Warning if already present will
      #be supressed
      dir.create(file.path(output.path, "array_plots"), showWarnings = FALSE)
      #Defining the dimensions and layout of the heatmap and its color key
      # 1. Heatmap
      # 2. Row dendrogram
      # 3. Column dendrogram
      # 4. Key
      #
      # Layout here:
      # 0   3
      # 2   1
      # 0   4
      #
      #
      # --> lwid and lhei are relative values: 1 = 100%, 2 = 200%, etc.
      lmat <-  rbind(c(0,3),c(2,1),c(0,4)) 
      lwid <-  c(1, 5)
      lhei <-  c(0.75,6,1.5)
    
      #Opening tiff device
      tiff(
        paste0(output.path,"/array_plots/",Sys.Date(),"_array_plot_",tag,"_",
        idx,".tif"),
        width = 2000, 
        height = 5000, 
        compression = "lzw", 
        pointsize = 10, 
        res = 300
      )
      heatmap.2(
        plot.ma,              # Source data
        dendrogram="none",    # Col and Row Dendrogram         
        trace="none",         # Trace line in the main histogram
        Rowv=FALSE,           # Row clustering
        Colv=FALSE,           # Column clustering
        na.color="black",     # Color of NA values
        col=colpal,           # Color palette
        breaks = 50,          # Number of breaks between lowest & highest color
        lmat = lmat,          # Defined above
        lhei = lhei,          # Defined above
        lwid = lwid,          # Defined above
        main = paste0(array.label," Array Plot"), # Heatmap title
        margins = c(5,10),    # Label margins
        colsep = c((1:4)*block.dim), # Separation of Columns
        rowsep = c((1:12)*22),  # Separation of Rows
        sepwidth = c(0.0025, 0.0025),  # Separation width
        cexCol = 0.6,         # Margins for Col names
        cexRow = 0.6,         # Margins for Row names
        denscol = "black",    # Density map Histogram color
        keysize = 0.01,          # Color key size
        symkey = FALSE,       # Symmetry of color key around 0
        symbreaks = FALSE,    # Symmetry of breaks around 0
        na.rm = TRUE          # Remove NA values
      )
      dev.off()
    }else{
      #Defining the dimensions and layout of the heatmap and its color key
      lmat <-  rbind(c(0,3,0),c(2,1,0),c(0,4,0)) 
      lwid <-  c(0.7, 1, 0.3)
      #lhei <-  c(0.15,1,0.3)
      lhei <-  c(0.17,1,0.3)
      
      heatmap.2(
        plot.ma,              # Source data
        dendrogram="none",    # Col and Row Dendrogram         
        trace="none",         # Trace line in the main histogram
        Rowv=FALSE,           # Row clustering
        Colv=FALSE,           # Column clustering
        na.color="black",     # Color of NA values
        col=colpal,           # Color palette
        breaks = 50,          # Number of breaks between lowest & highest color
        lmat = lmat,          # Defined above
        lhei = lhei,          # Defined above
        lwid = lwid,          # Defined above
        main = paste0(array.label," Array Plot"), # Heatmap title
        margins = c(5,10),    # Label margins
        colsep = c((1:4)*block.dim), # Separation of Columns
        rowsep = c((1:12)*22),  # Separation of Rows
        sepwidth = c(0.0025, 0.0025),  # Separation width
        cexCol = 0.6,         # Margins for Col names
        cexRow = 0.6,         # Margins for Row names            
        denscol = "black",    # Density map Histogram color
        keysize = 1,          # Color key size
        symkey = FALSE,       # Symmetry of color key around 0
        symbreaks = FALSE,    # Symmetry of breaks around 0
        na.rm = TRUE          # Remove NA values
      )
    }
} 
#++++++++++++++++++++++++++++ plotArrayTiff ++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++++ plotArrayPng ++++++++++++++++++++++++++++++++++++
plotArrayPng <- function (plot.ma=NULL, tag=NULL, idx=NULL, array.label=NULL,
  block.dim=NULL, colpal="heat.colors", output.path=NULL) {
    
  if(!is.null(output.path)){
      #Create output directory if not present. Warning if already present will
      #be supressed
      dir.create(file.path(output.path, "array_plots"), showWarnings = FALSE)
      #Defining the dimensions and layout of the heatmap and its color key
      # 1. Heatmap
      # 2. Row dendrogram
      # 3. Column dendrogram
      # 4. Key
      #
      # Layout here:
      # 0   3
      # 2   1
      # 0   4
      #
      #
      # --> lwid and lhei are relative values: 1 = 100%, 2 = 200%, etc.
      lmat <-  rbind(c(0,3),c(2,1),c(0,4)) 
      lwid <-  c(1, 5)
      lhei <-  c(0.75,6,1.5)
    
      #Opening png device
      png(
        paste0(output.path,"/array_plots/",Sys.Date(),"_array_plot_",tag,"_",
        idx,".png"),
        width = 2000, 
        height = 5000,  
        pointsize = 10, 
        res = 300
      )
      heatmap.2(
        plot.ma,              # Source data
        dendrogram="none",    # Col and Row Dendrogram         
        trace="none",         # Trace line in the main histogram
        Rowv=FALSE,           # Row clustering
        Colv=FALSE,           # Column clustering
        na.color="black",     # Color of NA values
        col=colpal,           # Color palette
        breaks = 50,          # Number of breaks between lowest & highest color
        lmat = lmat,          # Defined above
        lhei = lhei,          # Defined above
        lwid = lwid,          # Defined above
        main = paste0(array.label," Array Plot"), # Heatmap title
        margins = c(5,10),    # Label margins
        colsep = c((1:4)*block.dim), # Separation of Columns
        rowsep = c((1:12)*22),  # Separation of Rows
        sepwidth = c(0.0025, 0.0025),  # Separation width
        cexCol = 0.6,         # Margins for Col names
        cexRow = 0.6,         # Margins for Row names
        denscol = "black",    # Density map Histogram color
        keysize = 0.01,          # Color key size
        symkey = FALSE,       # Symmetry of color key around 0
        symbreaks = FALSE,    # Symmetry of breaks around 0
        na.rm = TRUE          # Remove NA values
      )
      dev.off()
    }else{
      #Defining the dimensions and layout of the heatmap and its color key
      lmat <-  rbind(c(0,3,0),c(2,1,0),c(0,4,0)) 
      lwid <-  c(0.7, 1, 0.3)
      #lhei <-  c(0.15,1,0.3)
      lhei <-  c(0.17,1,0.3)
      
      heatmap.2(
        plot.ma,              # Source data
        dendrogram="none",    # Col and Row Dendrogram         
        trace="none",         # Trace line in the main histogram
        Rowv=FALSE,           # Row clustering
        Colv=FALSE,           # Column clustering
        na.color="black",     # Color of NA values
        col=colpal,           # Color palette
        breaks = 50,          # Number of breaks between lowest & highest color
        lmat = lmat,          # Defined above
        lhei = lhei,          # Defined above
        lwid = lwid,          # Defined above
        main = paste0(array.label," Array Plot"), # Heatmap title
        margins = c(5,10),    # Label margins
        colsep = c((1:4)*block.dim), # Separation of Columns
        rowsep = c((1:12)*22),  # Separation of Rows
        sepwidth = c(0.0025, 0.0025),  # Separation width
        cexCol = 0.6,         # Margins for Col names
        cexRow = 0.6,         # Margins for Row names            
        denscol = "black",    # Density map Histogram color
        keysize = 1,          # Color key size
        symkey = FALSE,       # Symmetry of color key around 0
        symbreaks = FALSE,    # Symmetry of breaks around 0
        na.rm = TRUE          # Remove NA values
      )
    }
} 
#+++++++++++++++++++++++++++++ plotArrayPng ++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++++++ plotArray +++++++++++++++++++++++++++++++++++++
plotArray <- function(elist=NULL, idx=NULL, data.type="fg", log=NULL,
  normalized=NULL, aggregation=NULL, colpal="heat.colors",
  graphics.device="tiff", output.path=NULL) {  

  if(is.null(elist) || is.null(idx) || is.null(log) || is.null(normalized) ||
    is.null(aggregation)) {
        stop("ERROR: Not all mandatory arguments have been defined!")
  }
  
  n.blocks <- elist$printer$ngrid.c*elist$printer$ngrid.r
  n.block.c <- elist$printer$nspot.c
  n.block.r <- elist$printer$nspot.r
  
  if(idx == "all"){
      idx <- 1:ncol(elist$E)
  }
  
  if(aggregation=="min"){
      block.dim <- n.block.c/2
      row.len <- nrow(elist$C)
      col.len <- ncol(elist$C)
      tmp.col.len <- (row.len*col.len)/2
      elist$C[row(elist$C)[,1]%%2==1] <-
        matrix(apply(matrix(elist$C,2,tmp.col.len),2,min),row.len/2,col.len)
      elist$C <- elist$C[row(elist$C)[,1]%%2==1,]
      elist$Cb[row(elist$Cb)[,1]%%2==1] <-
        matrix(apply(matrix(elist$Cb,2,tmp.col.len),2,min),row.len/2,col.len)
      elist$Cb <- elist$Cb[row(elist$Cb)[,1]%%2==1,]
      elist$cgenes <- elist$cgenes[row(elist$cgenes)[,1]%%2==1,]
      
      ref.df <- data.frame("Block"=rep(1:n.blocks, each=n.block.c*n.block.r),
        "Column"=1:n.block.c, "Row"=rep(1:n.block.r, each=n.block.c))
      ref.df <- ref.df[row(ref.df)[,1]%%2==1,]
  }else if(aggregation=="mean"){
      block.dim <- n.block.c/2    # NEW!!!
      elist$C[row(elist$C)[,1]%%2==1,] <-
        (elist$C[row(elist$C)[,1]%%2==1,]+elist$C[row(elist$C)[,1]%%2==0,])/2
      elist$C <- elist$C[row(elist$C)[,1]%%2==1,]
      elist$Cb[row(elist$Cb)[,1]%%2==1,] <-
        (elist$Cb[row(elist$Cb)[,1]%%2==1,]+
        elist$Cb[row(elist$Cb)[,1]%%2==0,])/2
      elist$Cb <- elist$Cb[row(elist$Cb)[,1]%%2==1,]
      elist$cgenes <- elist$cgenes[row(elist$cgenes)[,1]%%2==1,]
      
      ref.df <- data.frame("Block"=rep(1:n.blocks, each=n.block.c*n.block.r),
        "Column"=1:n.block.c, "Row"=rep(1:n.block.r, each=n.block.c))
      ref.df <- ref.df[row(ref.df)[,1]%%2==1,]
  }else if(aggregation=="none"){
      block.dim <- n.block.c
      ref.df <- data.frame("Block"=rep(1:n.blocks, each=n.block.c*n.block.r),
        "Column"=1:n.block.c, "Row"=rep(1:n.block.r, each=n.block.c))
  }else{
      stop("ERROR: Unknown aggregation approach!")
  }
  
  if(is.null(output.path) && length(idx)>1){
      par(ask=TRUE)
  }
  
  #For any number of idx  
  for(i in 1:length(idx)) {    
    array.label <- colnames(elist$E)[idx[i]]   
    genes.merged <- rbind(elist$genes, elist$cgenes)                   
         
    if(normalized==TRUE && data.type=="fg") {
      #During normalization the data is logarithmized. 
      #For a proper plot, control data has to be logarithmized too. 
      data.merged <- c(elist$E[,idx[i]], log2(elist$C[,idx[i]]))
    } else if(normalized==FALSE && data.type=="fg"){
      data.merged <- c(elist$E[,idx[i]], elist$C[,idx[i]])       
    } else if(normalized==FALSE && data.type=="bg"){
      #Load background data
      data.merged <- c(elist$Eb[,idx[i]], elist$Cb[,idx[i]])
    }
    
    transfer.list <- cbind(genes.merged, data.merged)                    
    colnames(transfer.list)[ncol(transfer.list)] <- array.label            
    tl <- transfer.list  
  
    #Merge array data into reference data frame
    #alternative UTF-workaround: "sort=FALSE" added to keep rows orders
    full.df <- merge(ref.df, tl, by = c("Block", "Column", "Row"), all.x = TRUE)
  
    #Modulo used to define first and last block in a row
    # UTF-workaround: instead of x <- unique(full.df$Block)%%4 
    x <- sort(unique(full.df$Block))%%4
    # ++++++++++++++++++++++++++++++++++++++++++
    
    blocks <- 1   #Loop counter for blocks
    block.rows <- 1   #Loop counter for block.rows
    #Build blocks with 22 rows and 22 columns and append 4 of them for 
    #each block.row
    while(blocks <= length(unique(full.df$Block))) {
      tmp.block.ma <- matrix(full.df[which(full.df$Block == blocks), 
        array.label], byrow = FALSE, ncol = block.dim)
      
      #If first block in a row, create new row  
      if(x[blocks] == 1) {
        tmp.block.row.ma <- tmp.block.ma  
      } else { 
        tmp.block.row.ma <- cbind(tmp.block.row.ma, tmp.block.ma)
      }
      
      #If its the 4th block in a row, once create plot.ma, after that 
      #just append the last row to plot.ma
      if(x[blocks] == 0) {
        if(block.rows == 1) {
          plot.ma <- tmp.block.row.ma
        } else {   
          plot.ma <- rbind (plot.ma, tmp.block.row.ma)
        }        
        block.rows <- block.rows + 1
      }
      blocks <- blocks + 1
    }
    
    #Log unlogarithmized Data when log=FALSE
    #Careful when already normalized
    if(log == FALSE) {                                      
      plot.ma <- log2(plot.ma)
    }
    
    if(graphics.device=="tiff"){
        #Call plotArrayTiff with the aquired data.
        plotArrayTiff(plot.ma, tag=data.type, idx=idx[i],
          array.label=array.label, block.dim=block.dim, colpal=colpal,
          output.path=output.path)
    }else if(graphics.device=="png"){
        #Call plotArrayPng with the aquired data.
        plotArrayPng(plot.ma, tag=data.type, idx=idx[i],
          array.label=array.label, block.dim=block.dim, colpal=colpal,
          output.path=output.path)
    }else{
    
    }
  }  
}
#+++++++++++++++++++++++++++++++ plotArray +++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++ yellowblue ++++++++++++++++++++++++++++++++++++
bluered <- function(n){
        #colorpanel(n, "blue", "red")
        colorpanel(n, "blue", "white", "red")
}
#++++++++++++++++++++++ yellowblue ++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++ blueyellow ++++++++++++++++++++++++++++++++++++
blueyellow<-function(n){
        colorpanel(n, "blue", "yellow")
}
#++++++++++++++++++++++ blueyellow ++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++ plotFeaturesHeatmap.2 ++++++++++++++++++++++++++++++++++

plotFeaturesHeatmap.2 <- function(features=NULL, elist=NULL, n1=NULL, n2=NULL,
  output.path=NULL, description=FALSE){

    if(is.null(features) || is.null(elist) || is.null(n1) || is.null(n2)) {
        stop("ERROR: Not all mandatory arguments have been defined!")
    }
    
    features <- c(na.exclude(features))
    datamatrix <- elist$E
    brc <- paste(elist$genes$Block,elist$genes$Row,elist$genes$Column)
    rows <- brc %in% features
    datamatrix <- datamatrix[rows,]
    if(description){
        rownames(datamatrix) <- elist$genes$Description[rows]
    }else{
        rownames(datamatrix) <- brc[rows]
    }

    #title <- ""
    if(nrow(datamatrix) > 100) {
        cexRowSize <- 0.1
    } else if(nrow(datamatrix) > 40) {
        cexRowSize <- 0.4
    } else {
        cexRowSize <- 0.8
    }
    #ColSideColors <- c(rep("firebrick",n1), rep("green4",n2))
    #ColSideColors <- c(rep("azure4",n1), rep("gold4",n2))
    ColSideColors <- c(rep("gray25",n1), rep("gray75",n2))
    
    #my.dist <- function(x) as.dist(1 - abs(cor(t(x), method="pearson")))
    my.dist <- function(x) as.dist((1-cor(t(x)))/2)
      #"euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
    #my.hclust <- function(d) hclust(d, method="complete")
      #"ward", "single", "complete", "average", "mcquitty",
      #"median" or "centroid"

    #Defining the dimensions and layout of the heatmap and its color key
    lmat <- rbind(c(5,4),c(0,1),c(3,2))  
    lwid <- c(1.5,4)
    lhei <- c(1.5,0.2,4)
    
    if(!is.null(output.path)){
        png(
            paste0(output.path,"/biomarker_heatmap2.png"),
            width = 2000,
            height = 2000,
            pointsize = 10,
            res = 300
        )   
            heatmap.2(datamatrix,
                distfun=function(x) as.dist((1-cor(t(x)))/2),
                Rowv = TRUE,                  #row and col clustering
                Colv = TRUE,
                ColSideColors=ColSideColors,
                col=blueyellow(300),          #color scheme
                #col=bluered(300),          #color scheme
                scale="row",                  #scale by column
                #main=title,                   #title of the heatmap
                margins = c(15, 15),           #margins around heatmap
                key=TRUE,                     #legend is present
                key.title="",                 #no title for legend
                key.xlab="",                  #no label for legend axis
                key.ylab="",
                keysize=1.1,                  #size of legend strip
                symkey=FALSE,                 #colors not symmetrical around 0
                density.info="density",       #no density information
                trace="none",                 #
                cexRow=cexRowSize,            #column labels' size
                labRow=substr(row.names(datamatrix),1,35), #row labels
                labCol=colnames(datamatrix),   #column labels: limit to 30 chars
                lmat=lmat,
                lhei=lhei,
                lwid=lwid
            )
            #normalized device coordinates (NDC): c(xmin, xmax, ymin, ymax)
            par(fig=c(0,0.975,0,0.94), new=TRUE)
            #xmin <- par("usr")[1]
            #xmax <- par("usr")[2]
            #ymin <- par("usr")[3]
            #ymax <- par("usr")[4]
            legend("topright",      # legend location
                legend = c(unique(elist$targets$Group)), # category labels
                col = c("gray25", "gray75"),  # color key
                pch = c(15,15),
                #cex=1.75,
                bty="n",
                horiz=TRUE
            )   
        dev.off()
    } else {
        heatmap.2(datamatrix,
            distfun=function(x) as.dist((1-cor(t(x)))/2),
            Rowv = TRUE,                  #row and col clustering
            Colv = TRUE,
            ColSideColors=ColSideColors,
            col=blueyellow(300),          #color scheme
            #col=bluered(300),           #color scheme
            scale="row",                  #scale by column
            #main=title,                   #title of the heatmap
            margins = c(15,15),           #margins around heatmap
            key=TRUE,                     #legend is present
            key.title="",                 #no title for legend
            key.xlab="",                  #no label for legend axis
            key.ylab="",
            keysize=1.1,                  #size of legend strip
            symkey=FALSE,                 #colors not symmetrical around 0
            density.info="density",       #no density information
            trace="none",                 #
            cexRow=cexRowSize,            #column labels' size
            labRow=substr(row.names(datamatrix),1,35), #row labels
            labCol=colnames(datamatrix),   #column labels: limit to 30 chars
            lmat=lmat,
            lhei=lhei,
            lwid=lwid
        )
        #normalized device coordinates (NDC): c(xmin, xmax, ymin, ymax)
        par(fig=c(0,0.975,0,0.94), new=TRUE)
        #xmin <- par("usr")[1]
        #xmax <- par("usr")[2]
        #ymin <- par("usr")[3]
        #ymax <- par("usr")[4]
        legend("topright",      # legend location
            legend = c(unique(elist$targets$Group)), # category labels
            col = c("gray25", "gray75"),  # color key
            pch = c(15,15),
            #cex=1.75,
            bty="n",
            horiz=TRUE
        )   
    }

}

#++++++++++++++++++++++ plotFeaturesHeatmap.2 ++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++ batchFilter.anova +++++++++++++++++++++++++++++++++
batchFilter.anova <- function(elist=NULL, log=NULL, p.thresh=0.05,
  fold.thresh=1.5, output.path=NULL){
        
    if(is.null(elist) || is.null(log)) {
        stop("ERROR: Not all mandatory arguments have been defined!")
    }
    
    if(log==FALSE){
        E <- log2(elist$E)
    }else{
        E <- elist$E
    }
    
    row.number <- nrow(E)
    col.number <- ncol(E)
    discard <- c()
    output <- matrix(nrow=row.number, ncol=4)
    colnames(output) <- c("BRC", "ANOVA p-value", "Fold Change", "Description")
    #how many batches and batch combinations are there?
    batches <- unique(elist$targets$Batch)
    Nbatches <- length(batches)
    cat(c("Found batches:", Nbatches, "\n"))
    #if less than two batches, return null
    
    if(Nbatches < 2) {
        return(NULL)
    }
        
    batch.combs<-combn(Nbatches,2)
        
        
    for(zeile in 1:row.number) {
        #make df for anova analysis
        data <- data.frame(Value=E[zeile,], Batch=elist$targets$Batch)
        
        #perform anova on readings vs batch
        p <- oneway.test(Value~Batch,data)$p.value
        if(is.na(p)) {
            p <- 1
        }
          
        #get the means for each batch
        means <- rep(NA,Nbatches)
        for(batchInd in 1:Nbatches){
            means[batchInd]<-mean(data[data$Batch==batches[batchInd],"Value"])
        }
                
        #cycle through all combinations, find fold changes
        #and store them in a separate vector
        FCvector <- rep(NA,ncol(batch.combs))
        for(batchInd in 1:ncol(batch.combs)){
            mean1 <- means[batch.combs[1,batchInd]]
            mean2 <- means[batch.combs[2,batchInd]]
            #determine FC
            FCvector <- 2^(mean1 - mean2)
                        # --> data must be in log2 scale!
                        
        }
        #if max FC value > then 1/minFC,
        #take it as FC, else - otherwise
        FC <- ifelse(max(FCvector) > 1/min(FCvector), max(FCvector),
         min(FCvector))
                
         #prep output table
         output[zeile,] <- c(paste0(elist$genes[zeile,1], ",",
          elist$genes[zeile,3], ",",elist$genes[zeile,2]),p, FC,
          elist$genes[zeile,4])
                
         #determine if the feature (row) is to be rejected
         if(p < p.thresh && (FC > fold.thresh || FC < 1/fold.thresh)) {
            discard <- c(discard,zeile)
         }
    }
    message(paste0("batchFilter.anova - number of features to discard: ",
     length(discard)), "\n")
        
    if(is.null(output.path)) {
        plot(x=log2(as.numeric(output[,3])), 
         y=-log10(as.numeric(output[,2])),
         xlab="log2(fold change)", 
         ylab="-log10(p-value)",
         main="Batch Filter Volcano")
                
        if(length(discard)>0){
            points(x=log2(as.numeric(output[discard,3])),
             y=-log10(as.numeric(output[discard,2])), 
             col="red")
        }
                
        abline(h=-log10(p.thresh), lty=2, col="red")
        abline(v=log2(fold.thresh), lty=2, col="red")                                   
        abline(v=log2(1/fold.thresh), lty=2, col="red")
        if(length(discard) > 0){
            elist <- elist[-discard,]
        }
    }else{
        write.table(x=output, 
         file=paste0(output.path,"/batch_filter.txt"), 
         sep="\t", eol="\n",
         row.names=FALSE,
         quote=FALSE)
                
        write.table(x=output[discard,], 
          file=paste0(output.path,"/batch_filter_discarded.txt"),
          sep="\t", eol="\n",
          row.names=FALSE,
          quote=FALSE)
                
        tiff(paste(output.path,"/batch_filter_before.tiff",sep=""), 
         width=2000,
         height=2000,
         pointsize=15, 
         compression="lzw", 
         res=300)
                
            plot(x=log2(as.numeric(output[,3])),
             y=-log10(as.numeric(output[,2])),
             xlab="log2(fold change)",
             ylab="-log10(p-value)", 
             main="Batch Filter Volcano")
                
            if(length(discard)>0){
                points(x=log2(as.numeric(output[discard,3])),
                y=-log10(as.numeric(output[discard,2])), 
                col="red")
            }
            abline(h=-log10(p.thresh), lty=2, col="red")
            abline(v=log2(fold.thresh), lty=2, col="red")
            abline(v=log2(1/fold.thresh), lty=2, col="red")
        dev.off()
                
        if(length(discard) > 0){
            tiff(paste(output.path,"/batch_filter_after.tiff",sep=""),
             width=2000,
             height=2000, 
             pointsize=15,
             compression="lzw",
             res=300)
                        
                plot(x=log2(as.numeric(output[-discard,3])),
                 y=-log10(as.numeric(output[-discard,2])),
                 xlab="log2(fold change)",
                 ylab="-log10(p-value)",
                 main="Batch Filter Volcano")
            dev.off()
            elist <- elist[-discard,]            
        }
    }
    return(elist)
}
#+++++++++++++++++++++++++++ batchFilter.anova +++++++++++++++++++++++++++++++++