[![DOI](https://zenodo.org/badge/397762008.svg)](https://zenodo.org/badge/latestdoi/397762008)
# Tooth_sciRNAseq
This repository contains the R codes used for the paper (Human iPSC Derived Enamel Organoid Guided by Single Cell Atlas of Human Tooth Development)

## System requirements:
All codes were tested on Windows 10 (R version 4.1.1) and Linux running R (version	3.6.3), and Python (version 3.7).

## Installation:
### Monocle3
To install monocle3 Please follow the instructions provided in the Trapnell-lab website https://cole-trapnell-lab.github.io/monocle3/docs/installation/

### scMLnet 
Please follow the instructions provided by this package authors to install the python and R components https://github.com/SunXQlab/scMLnet

### talklr
```r
library(devtools)
install_github("yuliangwang/talklr")
```

### other dependencies 
```r
install.packages(c('Seurat', 'treemap', 'networkD3','hrbrthemes', 'viridis', 'patchwork', 'circlize','tidyverse','tidyr','rliger','pheatmap','stringr', 'igraph','RColorBrewer','gridExtra','reshape2','ggtext'))
```
### Other Bioconductor dependencies 
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install(version = "3.10")
BiocManager::install(c("DEsingle","ComplexHeatmap","ViSEAGO","simplifyEnrichment","BiocParallel"))
```

## Instructions
The instructions to run each function are written as comments in the each R script file.

### Top_pathway.R
This file contains the main workflow to identify the top pathway activities at each step of ameloblast differentiation.

<img src=https://github.com/Rouhola-Baker-lab/Tooth_sciRNAseq/blob/be678dcd9d7f5da0ae4d16794b86ae8b5057c710/workflow.PNG>

### Customized_heatmap.R
This file contains the instruction to use the customized heatmap function that utilizes the complexHeatmap package, and combines key goterms and age_score calculations per cluster.

### Realtime_heatmap.R
This file contains the instructions to use the realtime heatmap function to show expression of a single gene over timepoints that can be compared across clusters.

### RNAScope_AM.R & RNAScope_OB.R
These files contain the code to generate the RNAScope maps for both ameloblasts and odontoblasts.

