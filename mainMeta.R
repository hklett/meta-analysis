# load preprocessed data sets
load("all.RData")
# source functions
source("SVM_functions.R")
source("metaCV.R")

#set seed
set.seed(13)
# calculate meta-analysis and set parameters
out <- metaCV(data.list=all, # data sets
              feature.select=logFC.feature.select, # function for feature rankings
              sizes=seq(5, 50, by=5), # tested sizes of genes in inner LOOCV for number of gene optimization
              metric="AUC", # performance metric according to which parameters and number of genes are optimized
              parameterGrid=-15:15, # parameter space for SVM parameters (C and gamma = 2^(parameterGrid))
              kfold=10, # folds in cross-validation
              CV="STRATIFIED", # type of cross-validation, e.g. LOOCV or STRATIFIED 
              verbose=F)

# annotate entrez IDs with gene symbol and gene name according to org.Hs.eg.db
features=add_annot(unique(unlist(out$feats.used)))

# save output
save(out, features, file="output.RData")