setwd("~/data_sci_RNA_seq")


library(monocle3)
library(DEsingle)

mixed_cds <- readRDS("/home/ubuntu/data_sci_RNA_seq/mixed_cds_simp.rds")

pair.list <- list(c("OE","DE-Prog"),c("DE-Prog","OEE"),c("OEE","PA-1"),c("PA-1","PA-2"),
                  c("PA-2","AM-1"),c("AM-1","AM-2"))

pair.list <- list(c("OEE","PA"),
                  c("PA","Am"))

mixed_cds <- mixed_cds[, colData(mixed_cds)$age_group == "20_22w"]
pair.list <- list(
                  c("PA","Am"))

mixed_cds <- mixed_cds[, colData(mixed_cds)$age_group == "17_19w"]
pair.list <- list(
  c("OEE","PA"))

pair.list <- list(
  c("OE","Am"),
  c("DE-Prog","Am"))

pair.list <- list(
  c("OE","AM-2"),
  c("DE-Prog","AM-2"))

pair.list <- list(
  c("AM-1","AM-2"))
for (pair in pair.list){
  temp_cds <- mixed_cds[,colData(mixed_cds)$assigned_cell_type.simp %in% pair]
  group = factor(colData(temp_cds)$assigned_cell_type.simp)
  results <- DEsingle(counts = counts(temp_cds), group = group, parallel = T)
  results.classified <- DEtype(results = results, threshold = 0.1)
  results.sig <- results.classified[results.classified$pvalue < 0.01, ]
  deg <- results.sig[,c(12,13,11,20,21,23,24)]
  deg$Gene <- rowData(temp_cds[rownames(deg),])$gene_short_name
  deg <- deg[,c(8,1,2,3,4,5,6,7)]
  deg$norm_total_mean_1 <- deg$norm_total_mean_1 *1000
  deg$norm_total_mean_2 <- deg$norm_total_mean_2 *1000
  colnames(deg) <- c("Gene","Cluster1_total_mean", "Cluster2_total_mean", "foldChange" , "pvalue", "pvalue.adj.FDR","Type", "State")
  deg[deg$foldChange > 1,]$State <- "up"
  deg[deg$foldChange < 1,]$State <- "down"
  saveRDS(deg,paste( "deg_0.01_notADJ_",paste(pair,collapse = '_'), ".rds",sep =''))
}

