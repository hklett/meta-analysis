# load required functions
source("SVM_functions.R")
source("tune.R")
source("metaCV.R")
source("singleCV.R")

load("all.RData")
# extract single data set M1
M1 <- all[[1]]
# run CV with feature selection on single data set
M1_out <- single_CV_wFeature( M1, 
                              logFC.feature.select,
                              CV="STRATIFIED",
                              sizes= 1:10,
                              metric="AUC", 
                              kfold=10, 
                              SEED=13, 
                              verbose=F, 
                              parameterGrid=-15:15)

save(M1_out, file="M1_out.RData")

