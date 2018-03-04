#######################################
###     META-ANALYSIS ALGORITHM     ###
#######################################

metaCV <- function(data.list, 
                   feature.select, 
                   sizes=c(10, 50, 100), 
                   metric, 
                   rand.samples=F, 
                   num.mismatches=1,
                   ...){
  ## data.list: list of data sets, where within each data set columns represent samples and rows features. Sample names should be labeled such as the first letter of the test group contains a "T" (tumor) and reference group a "N" (normal)
  ## feature.select: function on how to select features, e.g., logFC.feature.select
  ## sizes: number of top genes tested in gene number optimization (e.g., 5, 10, ..., 50)
  ## metric: the metric used for performance evaluation, e.g., AUC (area under the curve of the receiver operating characterisitc), ACC (accuracy), BA (balanced accuracy), SENS (sensitivity), SPEC (specificity)
  ## rand.samples: should be false. Only T for testing for potential bias in method
  ## num.mismatches: number of mismatches allowed, i.e., if num.mismatches = 1, features present in N-1 data sets are still consideredfor further analysis. If features are present in < N-1 data sets, they will be removed from analysis.
  ## ... further arguments passed to pdac_svm() (see SVM_functions.R)
  
  # randomize column names if rand.samples = True
  if(rand.samples) data.list <- lapply(data.list, function(x){ colnames(x)=sample(colnames(x)); return(x)})
  results <- list()
  # outer LOOCV: split data list in training sets and testing set
  for(i in 1:length(data.list)){
    cat("Testing: ", names(data.list)[i], "\n")
    train  <- data.list[-i]
    test  <- data.list[[i]]
    
    #rank features according to feature selection on each single training set
    rank.list <- feature.select(train)
    
    #Optimize the number of selected genes from rank.list in inner CV
    #loop over different sizes n of features
    count=1
    out.list <- list()
    for(m in sizes){
      cat("# features =", m, "\n")
      #inner CV for each n on training set
      innerCV <- lapply(1:length(train), function(j){
        cat("Testing in training set: ", names(train)[j], "\n")
        #get features of size n on inner training set
        feat.list <- feature.rank(rank.list[-j], num.mismatches=num.mismatches)
        inner_test <- train[[j]]
        idx <- na.omit(match(feat.list, rownames(inner_test)))[1:m]
        inner_test <- inner_test[idx, , drop=F]
        #test performance on outer training set
        out <- pdac_svm(inner_test, metric=metric, ...) 
        return(out) 
      })
      #get average performance over testing sets in inner CV
      mlower <- tolower(metric)
      perf  <-  sapply(innerCV, function(x) eval(parse(text=paste0("x$", mlower))))
      cat(perf, "\n")
      ave_perf <- mean(perf)
      cat("# of features = ", m, ":", metric, " = ", ave_perf, "\n")
      names(ave_perf) <- paste0("featNum = ", m)
      out.list[[count]] <- list(perf=ave_perf, m=m, innerCV=innerCV)   
      count <- count+1
    }
    #select top n features with the best average performance in the inner LOOCV
    idx_max <- which.max(sapply(out.list, function(x) x$perf))
    M  <- out.list[[idx_max]]$m
    cat("Loop", i, "# of features = ", M, "\n")
    feats  <- feature.rank(rank.list, num.mismatches=1)
    idx_feats <- na.omit(match(feats, rownames(test)))[1:M]
    test  <-  test[idx_feats, , drop=F]
    #test performance on test set in outer LOOCV
    pdac_out  <- pdac_svm(test, metric=metric, ...)
    results[[i]] <- list(M=M, feats.used = rownames(test), feats.chosen= feats[1:M], test.data = pdac_out, mLoop=out.list)
  }
  ### format results for output 
  names(results) <- names(data.list)
  # list of features
  feats.used <- lapply(results, function(x) x$feats.used)
  feats.chosen <- lapply(results, function(x) x$feats.chosen)
  cat("Total number of unique features: ", length(unique(unlist(feats.used))), "\n")
  cat("Duplicated number of features: ", sum(duplicated(unlist(feats.used))), "\n")
  # number of selected genes for each testing data set in outer LOOCV
  Mvals  <- sapply(results, function(x) x$M) 
  cat("Average number of features = ", mean(Mvals), "\n")
  # Performance
  perf  <- sapply(results, function(x) eval(parse(text=paste0("x$test.data$", mlower))))
  cat("Average performance = ", mean(perf), "\n")
  return(list(feats.used=feats.used, feats.chosen=feats.chosen, Mvals = Mvals, perf=perf, results = results, metric=metric))
}

######################################
###  Feature selection functions   ###
######################################

# rank features according to absolute logFC differences between tumor and non-tumor tissues
logFC.feature.select <- function(data.list){
  limma.lists <- lapply(data.list, function(x){
    design <- create_obj(x, "design")
    deg <- diff_analysis(x, design, ann=F, sort="logFC")
    rownames(deg)
  })  
}

# rank features according to p-value between tumor and non-tumor tissues
pval.feature.select <- function(data.list){
  limma.lists <- lapply(data.list, function(x){
    design <- create_obj(x, "design")
    deg <- diff_analysis(x, design, ann=T, sort="p")
    rownames(deg)
  })  
}

# rank features randomly
rand.feature.select <- function(data.list){
  lapply(data.list, function(x) sample(rownames(x)))
}

# rank features according to SVM weights
svm.weights.feature.select <- function(data.list){
  SVMs <- lapply(data.list, function(x) svm(t(x), y=getY(x), kernel="linear", cost=1))
  w <- lapply(SVMs, function(x) t(x$coefs) %*% x$SV)
  w.sort <- lapply(w, function(x){y <- abs(t(x)); gsub("X", "", rownames(y)[order(y[,1], decreasing=T)])})
  w.sort
}

####################################
###    combine feature ranks     ###
####################################

feature.rank <- function(rank.list, num.mismatches){
  items <- names(table(unlist(rank.list)))[table(unlist(rank.list)) >= (length(rank.list) - num.mismatches)]
  ord <- do.call(cbind, lapply(rank.list, function(x) (match(items, x)-1)/(length(x)-1)))
  ave.rank <- apply(ord, 1, function(x) mean(x, na.rm=T)); names(ave.rank)=items
  sort.rank <- sort(ave.rank)
  names(sort.rank)
}

###############################
###    do LIMMA analysis    ###
###############################

diff_analysis <- function(limma, design, ann=T, method=c("ls", "robust")[1], sort="p", n=nrow(limma)){
  fit <- lmFit(limma, design, method=method)
  fit2 <- eBayes(fit)
  genes <- topTable(fit2, coef=2, adjust="BH",number=n, sort.by=sort)
  if(ann) genes <- add_annot(genes)
  #print(dim(genes))
  return(genes)
}

# create design matrix for limma analysis
create_obj <- function(limma, object=c("Type", "design")[1]){
  Type <- c()
  Type[grep("T", colnames(limma))] <- "T"
  Type[grep("N", colnames(limma))] <- "N"
  Type <- as.factor(Type)
  if(object=="Type") return(Type)
  if(object=="design") return(model.matrix(~Type))
}

########################################
###    add Annotation for features   ###
########################################
add_annot <- function(diffgenes, symbol=T, name=T, entrez=F, Ensembl=F){
  require(org.Hs.eg.db)
  if(is.null(dim(diffgenes))){
    names <- diffgenes
  }else{
    names <- rownames(diffgenes) 
  } 
  if(symbol) diffgenes <- data.frame(symbols= as.character(mget(names, org.Hs.egSYMBOL, ifnotfound=NA)), diffgenes)
  if(name) diffgenes <- data.frame(genename = as.character(mget(names, org.Hs.egGENENAME, ifnotfound=NA)), diffgenes)
  if(Ensembl)diffgenes <- data.frame(ensembl=as.character(sapply(mget(names, org.Hs.egENSEMBL, ifnotfound=NA), function(x) x[1])), diffgenes)
  if(entrez) diffgenes <- data.frame(Entrez_ID=names, diffgenes)
  return(diffgenes)
}

#########################################
###    plot performance vs. number of ###
###     features for meta-analysis    ###
#########################################
# plot.metaCV <- function(out){
#   #plot performance vs. feature size
#   perf <- lapply(out$results, function(x) sapply(x$mLoop, function(y) y$perf))
#   m <- lapply(out$results, function(x) sapply(x$mLoop, function(y) y$m))
#   pdf("Perf_vs_Size.pdf")
#   for ( i in 1:length(perf)){
#     if(i==1)
#       plot(m[[i]], perf[[i]], type="b", ylab = paste0("Performance (", out$metric, ")"), xlab="# of feature genes", pch=19, col=i, ylim=c(min(unlist(perf))-0.01, max(unlist(perf))+0.01))
#     else
#       lines(m[[i]], perf[[i]], type="b", col=i, pch=19)    
#   }
#   legend("bottom", legend=names(out$results), lty=1, col=1:length(perf), cex=0.8)
#   dev.off()
#   #ROC curves
#   probs <- lapply(out$results, function(x) x$test.data$prob)
#   AUCs  <- sapply(probs, AUC)
#   labs <- lapply(probs, function(x) as.numeric(factor(substring(rownames(x), 1, 1), labels=c(0, 1))))
#   pred <- mapply(function(x, y) prediction(x[, 2], y), probs, labs, SIMPLIFY = F)
#   perf2 <- lapply(pred, function(x) performance(x, "tpr", "fpr"))
#   pdf("ROCcurves.pdf")
#   for(i in 1:length(perf2)){
#     plot(perf2[[i]], col=i, add=ifelse(i==1, FALSE, TRUE), lwd=2, lty=ifelse(i==6, 2, 1))
#   }
#   legend("bottomright", legend=paste0(names(probs), " (AUC=", round(AUCs, 2), ")"), lty=1, lwd=2, col=1:length(probs))
#   abline(a=0, b=1, col=1, lty=2)
#   dev.off()
#   
#   # performance boxplots
#   performances <- do.call(c, lapply(out$results, function(x) c(x$test.data$acc, x$test.data$auc, x$test.data$ba, x$test.data$spec, x$test.data$sens)))
#   dat <- as.factor(sapply(names(performances), function(x) substring(x, 1, nchar(x)-1)))
#   groups <- rep(c("acc", "auc", "ba", "spec", "sens"), length(out$results))
#   mycols=c("red", "green", "orange", "blue", "yellow")
#   library(lattice)
#   pdf("performance_plot.pdf")
#   lat <- barchart(performances~dat, groups=groups, col=mycols, key=list(text=list(c("Accuracy", "AUCROC", "Balanced Accuracy", "Specificity", "Sensitivity")), rectangles=list(col=mycols)))
#   print(lat)
#   dev.off()
# }

################################
###    get information about ###
###   secretion of proteins  ###  
################################

# # meta file containing secretion status of proteins is needed
# getSecret <- function(x, meta=meta){
#   require(org.Hs.eg.db)
#   protID <- mget(as.character(x), org.Hs.egUNIPROT)
#   m <- lapply(protID, function(x) meta[match(x, meta$V1), "V3"])
#   tmp <- mapply(function(x, y) paste(x,y), protID, m, SIMPLIFY=F)
#   colprotID <- sapply(tmp, function(x) paste(x, collapse="|"))
# }

#############################################
###    make ROC curves and calculate AUC  ###  
#############################################

# makeROC <- function(perf, folder, name, col=NULL, lty=NULL, lwd=NULL, legend.pos="bottomright", ...){
#   require(ROCR)
#   prob <- lapply(perf, function(x) x$prob)
#   labels <- lapply(prob, function(x) as.numeric(factor(substring(rownames(x), 1, 1), labels=c(0, 1))))
#   pred <- mapply(function(x, y) prediction(x[, 2], y), prob, labels, SIMPLIFY = F)
#   AUC <- sapply(pred, function(x) as.numeric(attributes(performance(x, measure="auc"))$y.values))
#   IntAUC <- mapply(function(x, y) confAUC(x, sum(y==1), sum(y==2)), AUC, labels)
#   perf2 <- lapply(pred, function(x) performance(x, "tpr", "fpr"))
#   pdf(paste0(folder, name, ".pdf"))
#   if(is.null(col)) col=1:length(perf2)
#   if(is.null(lty)) lty = 1
#   if(length(lty)==1) lty = rep(lty, length(perf))
#   if(is.null(lwd)) lwd = 1
#   for(i in 1:length(perf2)){
#     plot(perf2[[i]], col=col[i], add=ifelse(i==1, FALSE, TRUE), lwd=lwd, lty=lty[i], ...)
#   }
#   #legend("bottomright", legend=paste0(names(prob), " (AUC=", round(AUC, 2), ")"), lty=1, lwd=2, col=1:length(prob))
#   legend(legend.pos, legend=paste0(colnames(IntAUC), " AUC=", sprintf("%.2f", round(AUC, 2)), " [", apply(round(IntAUC, 2), 2, function(x) paste(sprintf("%.2f", x), collapse=",")), "]"), lty=lty, lwd=lwd, col=col)
#   abline(a=0, b=1, col="grey50", lty=2)
#   dev.off()
# }

# calculate 95% confidence interval for AUC
# confAUC <- function(AUC, n, m){
#   Q1 = AUC/(2-AUC)
#   Q2 = 2*AUC^2/(1+AUC)
#   std = sqrt((m*n)^-1 * (AUC * (1-AUC) + (m-1)*(Q1 - AUC^2) + (n-1) * (Q2 - AUC^2)))
#   Int = c(AUC - qnorm(0.975) * std, AUC + qnorm(0.975) * std)
#   if(Int[2]>1) Int[2] <- 1
#   Int
# }