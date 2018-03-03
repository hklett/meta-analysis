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

## Run meta-analysis
In R environment, load preprocessed data into working environment and run meta-analysis, e.g. in main.R. See function files for required data format and parameter specifications.

```
load("all.RData") 
source("main.R")
```

The output contains selected features and detailed results of cross-validation runs.