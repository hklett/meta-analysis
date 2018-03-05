### single study CV with inner feature selection ###

#CV
single_CV_wFeature <- function(mat, 
                               feature.select, 
                               CV=c("STRATIFIED", "LOOCV")[1], 
                               sizes=1:5, 
                               metric="AUC", 
                               kfold=10, 
                               SEED=13, 
                               verbose=F, 
                               parameterGrid=-15:15){
  # mat: p x n matrix, with n samples and p features column names should contain "T" for tumor and "N" for normal tissue
  # feature.select: function on how to select features, e.g., logFC.feature.select
  # CV: type of cross-validation. Stratified for balanced cross-validation, i.e., ratio of tumor and normal samples is the same in testing and training sets
  # sizes: number of top genes tested in gene number optimization
  # metric: the metric used for performance evaluation, e.g., AUC (area under the curve of the receiver operating characterisitc), ACC (accuracy), BA (balanced accuracy), SENS (sensitivity), SPEC (specificity)
  # kfold: folds of CV
  # SEED: random seed
  # parameterGrid: grid for C and gamma parameter in SVM classifier. Tested values = 2^(parameterGrid).
  
  if(CV == "STRATIFIED"){
    MAT  <- strat_split(mat, kfold, SEED)
    cores  <- ifelse(kfold < detectCores(), kfold, detectCores())
  } else if(CV == "LOOCV"){
    MAT  <- lapply(1:ncol(mat), function(x) mat[, x, drop=F])
    cores <- detectCores()
  }
  
  # outer kfold-CV: split data list in training sets and testing set
  cv_results <- lapply(1:length(MAT), function(i){
    train  <- do.call(cbind, MAT[-i])
    y_train <- getY(train)
    test  <- MAT[[i]]
    y_test <- getY(test)
    if(verbose) cat("Testing samples:", colnames(test), "\n")
    # obtain gene ranking from training samples according to feature select function
    rank.list <- unlist(feature.select(list(train)))
    
    #Optimize the number of selected genes from rank.list in inner kfold-CV
    #loop over different sizes n of features
    count=1
    out.list <- list()
    for(m in sizes){
      if(verbose) cat("Test size: m\n")
      # select the top m genes from gene ranking
      idx <- na.omit(match(rank.list, rownames(train)))[1:m]
      # evaluate the performance of these genes in a CV of the training samples (inner CV)
      out <- pdac_svm(train[idx, ,drop=F], parameterGrid= parameterGrid, CV = CV, metric = metric, kfold=kfold, SEED=SEED, verbose=verbose)
      # store the results
      mlower <- tolower(metric)
      perf  <-  eval(parse(text=paste0("out$", mlower)))
      out.list[[count]] <- list(out=out, m=m, feats=out$features, perf=perf)   
      count <- count+1
    }
    # select m with the best CV performance
    idx_max <- which.max(sapply(out.list, function(x) x$perf))
    feats <- out.list[[idx_max]]$feats
    # tune parameters
    obj <- tunecont2(train[feats, , drop=F], y_train, parameterGrid=parameterGrid, kfold=kfold, metric=metric, CV=CV)
    if(verbose) print(obj$best.parameters)
    # train SVM classifier on training samples
    SVM <- svm(t(train[feats, , drop=F]), y_train, gamma=obj$best.parameters$gamma, cost=obj$best.parameters$cost, probability=T)
    # predict testing samples using SVM model
    predict_out <- predict(SVM, t(test[feats, , drop=F]), probability = TRUE)
    prob <- attr(predict_out, "probabilities")
  return(list(pred=predict_out, y_test=y_test, prob=prob, feats=feats, out.list=out.list, obj=obj))
  })  
  # store results
  all_prob <- do.call(rbind, lapply(cv_results, function(x) x$prob))
  all_pred <- do.call(c, lapply(cv_results, function(x) as.character(x$pred)))
  all_test <- do.call(c, lapply(cv_results, function(x) as.character(x$y_test)))
  all_feats <- lapply(cv_results, function(x) x$feats)
  ave_acc <- ACC(all_test, all_pred)
  ave_auc <- AUC(all_prob)
  ave_ba <- BA(all_test, all_pred)
  ave_sens <- SENS(all_test, all_pred)
  ave_spec <- SPEC(all_test, all_pred)  
  return(list(acc=ave_acc, auc=ave_auc, ba=ave_ba, sens=ave_sens, spec=ave_spec, feats=all_feats, prob=all_prob, mat=mat, cv_results=cv_results))
}


#################################

