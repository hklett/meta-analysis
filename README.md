# meta-analysis
Meta-analysis for classifier identification between two cases (e.g. tumor vs. non-tumor tissues) of multiple transcriptome data sets


## Initialize
Clone repository from command line
```
git clone https://github.com/hklett/meta-analysis.git
```

In R environment, install R libraries needed for meta-analysis from Bioconductor and cran

```
source("https://bioconductor.org/biocLite.R")
biocLite("limma", "org.Hs.eg.db") # bioconductor packages
install.packages("foreach", "doMC", "e1071", "caret", "ROCR") #cran packages
```

## Run CV with feature selection in single data set
In R environment, load preprocessed data, e.g. M1 into working environment and run singleCV, e.g., in mainSingle.R. See function files for required data format and parameter specifications

```
source("mainSingle.R")
```

## Run meta-analysis of multiple data sets
In R environment, load preprocessed data into working environment and run meta-analysis, e.g. in mainMeta.R. See function files for required data format and parameter specifications.

```
source("mainMeta.R")
```

The output contains selected features and detailed results of cross-validation runs.