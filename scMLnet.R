setwd("~/data_sci_RNA_seq")
library(tidyr)
library(monocle3)
library(ggplot2)
library(dplyr)
library(Matrix.utils)


devtools::install_github("YUZIXD/scMLnet")
library(scMLnet)
library(Matrix)
library(parallel)
library(Seurat)

Amelo_cds_MTi <-readRDS("epi_cds_MTi_70_manhattan_8L_0.3_7k_pesudo_regressed_assigned_cell_type_7.rds")


# import sample data
GCMat <- counts(Amelo_cds_MTi)
row.names(GCMat)<- rowData(Amelo_cds_MTi)$gene_short_name
colnames(GCMat) = as.character(colnames(Amelo_cds_MTi))
rownames(GCMat) = as.character(rowData(Amelo_cds_MTi)$gene_short_name)

class(rownames(GCMat))=="character"


GCMat<- as(GCMat,"dgCMatrix")
GCMat<- as.matrix(GCMat)

# import sample annotation
clustering <- data.frame(id= colnames(Amelo_cds_MTi),Cluster = colData(Amelo_cds_MTi)$assigned_cell_type.4)
write.table(clustering,"./database/barcodetype.txt",sep = "\t")
BarCluFile <- "./database/barcodetype.txt"
BarCluTable <- read.table(BarCluFile,sep = "\t",header = TRUE,stringsAsFactors = FALSE)


types <- unique(BarCluTable$Cluster)

LigClu <- "SI-2"       #types[4]
RecClu <- "PA-1"     #types[8]

pval <- 0.1
logfc <- 1.0
LigRecLib <- "./database/LigRec.txt"
TFTarLib <- "./database/TFTargetGene.txt"
RecTFLib <- "./database/RecTF.txt"

netList <- RunMLnet(GCMat, BarCluFile, RecClu, LigClu, 
                    pval, logfc, 
                    LigRecLib, TFTarLib, RecTFLib)

workdir <- "sample"
DrawMLnet(netList,LigClu,RecClu,workdir,plotMLnet = F)

workdir <- "sample"
PyHome <- "/home/ubuntu/miniconda3/bin/python" #for Window
DrawMLnet(netList,LigClu,RecClu,workdir,PyHome,plotMLnet = T)
