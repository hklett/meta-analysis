#### libraries required for meta-analysis ####
library(limma)
library(org.Hs.eg.db)
library(foreach)
library(doMC)
library(e1071)
require(caret)
require(ROCR)
source("tune.R") # modified tune function from the e1071 package

# get class labels of data set. Here, first letter of sample name indicates class label (N=normal, T=tumor)
getY  <- function(mat){
  Type <- substring(colnames(mat), 1, 1)
  factor(Type, levels=c("N", "T"))
}

# split data in k groups
split_fun <- function(idx, k){ split(idx, sort(rank(idx) %% k)) }

# split data in kfold groups with same ratios between tumor and non-tumor samples (for balanced cross-validation)
strat_split  <- function(mat, kfold, SEED){
  set.seed(SEED)
  y <- getY(mat)
  k <- ifelse(min(table(y)) < kfold, min(table(y)), kfold)
  idx1 <- sample(which(y == levels(y)[1]))
  idx2 <- sample(which(y == levels(y)[2]))
  split1 <- split_fun(idx1, k)
  split2 <- split_fun(idx2, k)
  splits <- mapply(c, split1, split2, SIMPLIFY=FALSE)
  MAT <- lapply(splits, function(x) mat[, x, drop=F])
}

# split for balanced cross-validation in modified tune.R function
strat_tune  <- function(train.y, kfold){
  y <- train.y
  k <- ifelse(min(table(y)) < kfold, min(table(y)), kfold)
  idx1 <- sample(which(y == levels(y)[1]))
  idx2 <- sample(which(y == levels(y)[2]))
  split1 <- split_fun(idx1, k)
  split2 <- split_fun(idx2, k)
  splits <- mapply(function(x, y) setdiff(c(as.vector(unlist(split1)),as.vector(unlist(split2))), c(x,y)), split1, split2, SIMPLIFY=FALSE)
  return(splits)
}

# tunecont <- function(mat, y, parameterGrid, kfold, metric, CV){
#   cost <- 2^parameterGrid
#   ranges <- list(cost=cost)
#   sampling  <- ifelse(CV=="STRATIFIED", "strat", "cross")
#   kfold <- ifelse(sampling=="cross", ncol(mat), kfold)
#   if(sampling == "strat") err.fun <- function(true, pred) 1- eval(parse(text=paste0(metric, "(true, pred)")))
#   if(sampling == "cross") err.fun <- NULL
#   tunecontrol  <-  tune.control(sampling=sampling, nrepeat=1, cross = kfold, error.fun = err.fun)
#   obj <- tune(svm, t(mat), y, ranges=ranges, tunecontrol= tunecontrol, kernel="linear")
#   Cbar <- as.numeric(obj$best.parameters)
#   gamma <- 1/(2*exp(-log(cost)+log(as.numeric(Cbar))))
#   ranges <- list(cost=Cbar, gamma=gamma)
#   obj <- tune(svm, t(mat), y, ranges=ranges, tunecontrol= tunecontrol, kernel="radial")
#   return(obj)
# }

# tune control function for modified tune.R function
tunecont2 <- function(mat, y, parameterGrid, kfold, metric, CV){
  cost <- 2^parameterGrid
  gamma <- 2^parameterGrid
  ranges <- list(cost=cost, gamma=gamma)
  sampling  <- ifelse(CV=="STRATIFIED", "strat", "cross")
  kfold <- ifelse(sampling=="cross", ncol(mat), kfold)
  if(sampling == "strat"){
    if(metric=="AUC"){
      err.fun="AUC"
    }else{
      err.fun <- function(true, pred) 1- eval(parse(text=paste0(metric, "(true, pred)")))
    }
  } 
  if(sampling == "cross"){
    err.fun <-  NULL
  }
  tunecontrol  <-  tune.control(sampling=sampling, nrepeat=1, cross = kfold, error.fun = err.fun)
  obj <- tune(svm, t(mat), y, ranges=ranges, tunecontrol= tunecontrol)
  return(obj)
}

### performance metrics ###

# balanced accuracy
BA <- function(true, pred){
  as.numeric(confusionMatrix(pred, true, positive="T")$byClass[8])
}
# accuracy
ACC <- function(true, pred){
  as.numeric(confusionMatrix(pred, true, positive="T")$overall[1])
}
# sensitivity
SENS <- function(true, pred){
  as.numeric(confusionMatrix(pred, true, positive="T")$byClass[1])
}
# specificity
SPEC <- function(true, pred){
  as.numeric(confusionMatrix(pred, true, positive="T")$byClass[2])
}
# area under curve (AUC) of the receiver operating characteristic (ROC)
AUC <- function(prob){  
  labels <- as.numeric(factor(substring(rownames(prob), 1, 1), labels=c(0, 1)))
  PRED <- prediction(prob[, 2], labels)
  as.numeric(attributes(performance(PRED, measure="auc"))$y.values)
}

# evaluate cross-validation performance of a given data set using SVM classification
pdac_svm <- function(mat, 
         parameterGrid= -15:15, 
         CV = c("LOOCV", "STRATIFIED")[1], 
         metric = c("AUC", "ACC", "BA", "SENS", "SPEC")[1], 
         kfold=10, 
         SEED=13, 
         verbose=F){
  #mat: p x n matrix, with n samples and p features column names should contain "T" for tumor and "N" for normal tissue    
  # y: factor of binary classes
  #parameterGrid: grid for C and gamma parameter in SVM classifier. Tested values = 2^(parameterGrid).
  #CV: type of cross-validation. Stratified for balanced cross-validation, i.e., ratio of tumor and normal samples is the same in testing and training sets
  #metric: performance metric to which CV is evaluated
  #kfold: folds of CV
  #SEED: random seed 
  call <- match.call()
  if(min(table(getY(mat))) < 4) CV = "LOOCV"
  #split data set in training and testing samples
  if(CV == "STRATIFIED"){
    MAT  <- strat_split(mat, kfold, SEED)
    cores  <- ifelse(kfold < detectCores(), kfold, detectCores())
  } else if(CV == "LOOCV"){
    MAT  <- lapply(1:ncol(mat), function(x) mat[, x, drop=F])
    cores <- detectCores()
  }
  cat("CV=", CV, "\n")
  
  cv_results  <- mclapply(1:length(MAT), function(i){
    train  <- do.call(cbind, MAT[-i]); colnames(train) <- unlist(lapply(MAT[-i], colnames))
    y_train <- getY(train)
    test  <- MAT[[i]]
    y_test <- getY(test)
    if(verbose) cat("Testing samples:", colnames(test), "\n")
    obj <- tunecont2(train, y_train, parameterGrid=parameterGrid, kfold=kfold, metric=metric, CV=CV)
    if(verbose) print(obj$best.parameters)
    SVM <- svm(t(train), y_train, gamma=obj$best.parameters$gamma, cost=obj$best.parameters$cost, probability=T)
    svm.weight = t(SVM$coefs) %*% SVM$SV 
    PRED <- predict(SVM, t(test), probability=TRUE)
    prob <- attr(PRED, "probabilities")
    if(verbose) print(prob)
    return(list(pred=PRED, y_test=y_test, probabilities=prob, obj=list(best.parameters=obj$best.parameters, best.performance=obj$best.performance), svm.weight=svm.weight))
  }, mc.cores=cores)
  ### combine results ###
  all_prob <- do.call(rbind, lapply(cv_results, function(x) x$probabilities))
  all_pred <- do.call(c, lapply(cv_results, function(x) as.character(x$pred)))
  all_test <- do.call(c, lapply(cv_results, function(x) as.character(x$y_test)))
  all_para <- do.call(rbind, lapply(cv_results, function(x) x$obj$best.parameters))
  ave_weights  <- colMeans(do.call(rbind, lapply(cv_results, function(x) x$svm.weight)))
  ave_para <- colMeans(all_para)
  ave_acc <- ACC(all_test, all_pred)
  ave_auc <- AUC(all_prob)
  ave_ba <- BA(all_test, all_pred)
  ave_sens <- SENS(all_test, all_pred)
  ave_spec <- SPEC(all_test, all_pred)  
  return(list(acc=ave_acc, auc=ave_auc, ba=ave_ba, sens=ave_sens, spec=ave_spec, prob=all_prob, features=rownames(MAT[[1]]), num.feats = nrow(MAT[[1]]), parameters=ave_para, weights=ave_weights, call=call))
}
