# Tooth_sciRNAseq
This repository contains the R codes used for the paper (Human iPSC derived Ameloblast Differentiation guided by Single Cell Atlas of early Human Oral Development)

## System requirements:
All codes were tested on Windows 10 and Linux running R (version	3.6.3), and Python (version 3.7).

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
install.packages('Seurat', 'treemap', 'networkD3','hrbrthemes', 'viridis', 'patchwork', 'circlize','tidyverse','tidyr','rliger','pheatmap','stringr', 'igraph','RColorBrewer','gridExtra','reshape2','ggtext')
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

### Customized_heatmap.R
This file contains the instruction to use the customized heatmap function that utilizes the complexHeatmap package, and combines key goterms and age_score calculations per cluster.
