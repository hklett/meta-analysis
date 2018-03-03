# meta-analysis
Meta-analysis for classifier identification between two cases (e.g. tumor vs. non-tumor tissues) of multiple transcriptome data sets

## Run meta-analysis

Install libraries needed for meta-analysis from Bioconductor and cran

```
source("https://bioconductor.org/biocLite.R")
biocLite("limma", "org.Hs.eg.db") # bioconductor packages
install.packages("foreach", "doMC", "e1071", "caret", "ROCR") #cran packages
```

Load preprocessed data into working environment and run meta-analysis, e.g. in main.R 

```
load("all.RData") 
source("main.R")
```