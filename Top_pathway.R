
# load packages
library(dplyr)
library(talklr)
library(DEsingle)
library(scMLnet)
library(Matrix)
library(parallel)
library(igraph)
library(monocle3)
library(Seurat)

# load the customized functions.
source("top_path_function.R")
source("Run_scMLnet.R")
source("Draw_MLnet.R")
source("plot_lr_wiring.R")


#-------------------------------------------------------------------------------------------#
# This part where you eneter your data specific variables                                   #
#___________________________________________________________________________________________#

#This line is to make a copy of your cds, for simplicity
cds_base <- mixed_cds # replace mixed_cds with your cds name 

# the name of your clustering column or age_group column
clustering_column <- "assigned_cell_type.simp.2"
age_clustering_column <- "age_group"


#define your time point in your data
age_list <-  c("09_11w","12_13w","14_16w","17_19w")


# # I have separate analysis for 20_22w
# age_list <-  c("20_22w")

#this is to pull of the next age of target cell 
age_next_list <-  c("12_13w","14_16w","17_19w","20_22w")


# define your receiver cell of interest at each time point
receiver_cell <- list('09_11w' = "OE",
                      '12_13w' = "DE-Prog",
                      '14_16w' = "OEE",
                      '17_19w' = "PA")

# # I have separate analysis for 20_22w
# receiver_cell <- list('20_22w' = "AM-1")

# define your differentiation target cell type at each time point of analysis
diff_target_cell <- list('09_11w' = "DE-Prog",
                         '12_13w' = "OEE",
                         '14_16w' = "PA",
                         '17_19w' = "Am")

# # I have separate analysis for 20_22w
# diff_target_cell <- list('20_22w' = "AM-2")

# Define the group of cells that are present at each time point
in_groups <- list('09_11w' = c("OE" ,"DE-Prog",  "IK/EK", "D.mes"),
                  '12_13w' = c("DE-Prog", "IK/EK", "SR-1" , "D.mes"),
                  '14_16w' = c("IK/EK", "OEE" , "SI-1", "SI-2",  "D.mes", "Pre.O"),
                  '17_19w' = c( "SR-2" , "SI-2", "PA",  "Pre.O", "OD-1"),
                  '20_22w' = c("Am","SI-3", "OD-2","PA"))


# # # I have separate analysis for 20_22w
# in_groups <- list('20_22w' = c("AM-1", "SI-3", "AM-2","OD-2","PA-2"))



# define the pathways to include in the analsysis. some pathways are not major or not well known or there's no reliabe way to activate them in culture
path_included <- c("TGFb", "BMP", "BMP10", "GDF",   "GDNF", "NODAL", 
                   "ACTIVIN", "WNT", "ncWNT", "EGF", "NRG", "FGF", "PDGF", "VEGF", 
                   "IGF", "INSULIN",  "HH", "EDA", 
                   "NGF",    "NT", "FLT3", "HGF",  
                   "NOTCH", "NRXN", "OCLN", 
                   "ROBO")

# This is the full list of pathways in case you wanted to swap in something of interest. I don't recommed to include so many
# path_included <- c("TGFb", "BMP", "BMP10", "GDF",   "GDNF", "NODAL", 
#                    "ACTIVIN", "WNT", "ncWNT", "EGF", "NRG", "FGF", "PDGF", "VEGF", 
#                    "IGF", "INSULIN",  "HH", "LIFR", "EPO", "GH",  "TNF", "LT", "LIGHT",  "VEGI", "EDA", 
#                    "NGF",   "SPP1", 
#                    "ANGPTL", "ANGPT", "MK", "PTN", "PERIOSTIN", "AGT", "GPR", 
#                    "CCK", "EDN", "MELANOCORTIN", 
#                    "TAC", "UTS2", "AVP", "PRLH", "PARs", "PMCH", 
#                    "PROK", "PACAP", "VIP", "NPR1", "NPR2", "KIT", "GIPR", "FSH", 
#                    "LHB", "TSH", "NT", "FLT3", "HGF", "SEMA3", "CALCR", "ANNEXIN", 
#                    "APJ", "CRH", "ENHO", "GAS", "GCG", "GHRH", "GNRH", "GRN", "GUCA", 
#                    "RLN", "LEP", "GALECTIN", "NPS", "NPVF", "OSTN", "PROS", "PSAP", 
#                    "PTH", "QRFP", "CHEMERIN", "SAA", "UGRP1", "SCT", "SLURP", "BTLA", 
#                    "TRH", "UCN", "UROTENSIN", "BAG", "COLLAGEN", "FN1", "LAMININ", 
#                    "CHAD", "RELN", "THBS", "VTN", "TENASCIN", "NPNT", "DSPP", "VWF", 
#                    "BSP", "DMP1", "AGRN", "HSPG", "ADGRE5", "ALCAM", "ANXA1", "APP", 
#                    "CADM", 
#                    "CDH", "CDH1", 
#                    "CDH5", "CEACAM", "CLDN", "CLEC", "CNTN", "CSPG4", "DESMOSOME", 
#                    "EPGN", "EPHA", "EPHB", "ESAM", "GP1BA", "ICAM", "ICOS", "ITGB2", 
#                    "JAM", "L1CAM", "LCK", "MADCAM", "MAG", "MHC-I", "MHC-II", "MPZ", 
#                    "NCAM", "NECTIN", "NEGR", "NGL", "NKG2D", "NOTCH", "NRXN", "OCLN", 
#                    "PD-L1", "PDL2", "PECAM1", "PTPRM", "PVR", "SELE", "SELL", "SELPLG", 
#                    "SEMA4", "SEMA5", "SEMA6", "SEMA7", "SN", "THY1", "TIGIT", "VCAM", 
#                    "VISTA", "ROBO")

# here assigning colors for each pathway to be consistant. ONLY CHANGE IF YOU MODIFIEED THE DEFAULT PATHWAYs
palette <-c(TGFb = "#DBC6DE", NRG = "#E2E5D5", BMP10 = "#DA3BE0", GDF = "#679C8F", 
            GDNF = "#E2E2A1", NODAL = "#E24C66", ACTIVIN = "#C2E579", WNT = "#DA77E8", 
            ncWNT = "#64B1DF", EGF = "#66E1E1", BMP = "#E88E50", FGF = "#8580DF", 
            PDGF = "#62E4B7", VEGF = "#AAD8E1", IGF = "#D5AB90", INSULIN = "#7945DF", 
            HH = "#6EE07F", EDA = "#E2A1E0", NGF = "#8D984C", NT = "#DE8D9F", 
            FLT3 = "#70E041", HGF = "#D9E83D", NOTCH = "#8D92BF", NRXN = "#CC4E99", 
            OCLN = "#E1C959", ROBO = "#B0E7BE")


#set mode of the analysis:
subset_by_time_point <- T  # TRUE is probably the most accurate, but sometimes you don't have ennough cells, so you skip subsetting and let it use the whole clusters


# Create a directory to save data. 
out.file <- "top_path_output" 
user.file <- "Sesha"
meta.file <- paste(out.file,"/",user.file,"/metadata", sep="")
#dir.create(meta.file, recursive = T) #already created


#-------------------------------------------------------------------------------------------#
# no need to  modify anything below this line. Run as is.  or the heavy chunk first         #
#___________________________________________________________________________________________#

##### this part for differntail expresssion, It takes a lot of time. Might take from 15-1hour####
# only run once with highest computing # 
pair.list <- mapply(c, receiver_cell, diff_target_cell, SIMPLIFY = FALSE)

for (pair in pair.list){
  temp_cds <- cds_base[,colData(cds_base)[[clustering_column]] %in% unlist(pair)]
  group = factor(colData(temp_cds)[[clustering_column]])
  results <- DEsingle(counts = counts(temp_cds), group = group, parallel = T)
  results.classified <- DEtype(results = results, threshold = 0.1)
  results.sig <- results.classified[results.classified$pvalue < 0.05, ]
  deg <- results.sig[,c(12,13,11,20,21,23,24)]
  deg$Gene <- rowData(temp_cds[rownames(deg),])$gene_short_name
  deg <- deg[,c(8,1,2,3,4,5,6,7)]
  deg$norm_total_mean_1 <- deg$norm_total_mean_1 *1000
  deg$norm_total_mean_2 <- deg$norm_total_mean_2 *1000
  colnames(deg) <- c("Gene","Cluster1_total_mean", "Cluster2_total_mean", "foldChange" , "pvalue", "pvalue.adj.FDR","Type", "State")
  # deg[deg$foldChange > 1,]$State <- "up"
  # deg[deg$foldChange < 1,]$State <- "down"
  saveRDS(deg,paste(meta.file,"/","deg_0.05_",paste(unlist(pair),collapse = '_'), ".rds",sep =''))
}


# Get gene enrcichment for all celltype of interest
for (age in age_list){
  temp_cds <- cds_base[, colData(cds_base)[[age_clustering_column]] == age]
  included <- in_groups[age]
  temp_cds <- temp_cds[, colData(temp_cds)[[clustering_column]] %in% unlist(included)]
  if (subset_by_time_point == F){ temp_cds <- cds_base[, colData(cds_base)[[clustering_column]] %in% unlist(included)]}
  colData(temp_cds)[[clustering_column]] <- droplevels(colData(temp_cds)[[clustering_column]])
  
  GCMat <- counts(temp_cds)
  row.names(GCMat)<- rowData(temp_cds)$gene_short_name
  colnames(GCMat) = as.character(colnames(temp_cds))
  rownames(GCMat) = as.character(rowData(temp_cds)$gene_short_name)
  
  clustering <- data.frame(Barcode= colnames(temp_cds),Cluster = colData(temp_cds)[[clustering_column]])
  write.table(clustering,"database/barcodetype.txt",sep = "\t")
  BarCluFile <- "database/barcodetype.txt"
  BarCluTable <- read.table(BarCluFile,sep = "\t",header = TRUE,stringsAsFactors = FALSE)
  
  RecClu <- unlist(receiver_cell[age])    # the reciver
  LigClu <- setdiff(unlist(included), RecClu)     #all clusters except the recevier
  pval <- 0.05
  logfc <-0.15
  #get cores
  cores_num <- detectCores()
  cores = ifelse(cores_num>1,cores_num-1,1)
  RecClus <- getHighExpGene(GCMat,BarCluTable,RecClu,LigClu,pval,logfc,cores)
  saveRDS(RecClus, paste(meta.file,"/",RecClu,"_", age,"_enriched.rds", sep=""))
  
  
  # this is for target cells enrichment 
  next_age <- age_next_list[which(age_list == age)]
  temp_cds <- cds_base[, colData(cds_base)[[age_clustering_column]] == next_age]
  included <- in_groups[next_age]
  temp_cds <- temp_cds[, colData(temp_cds)[[clustering_column]] %in% unlist(included)]
  if (subset_by_time_point == F){ temp_cds <- cds_base[, colData(cds_base)[[clustering_column]] %in% unlist(included)]}
  colData(temp_cds)[[clustering_column]] <- droplevels(colData(temp_cds)[[clustering_column]])
  
  GCMat <- counts(temp_cds)
  row.names(GCMat)<- rowData(temp_cds)$gene_short_name
  colnames(GCMat) = as.character(colnames(temp_cds))
  rownames(GCMat) = as.character(rowData(temp_cds)$gene_short_name)
  
  clustering <- data.frame(Barcode= colnames(temp_cds),Cluster = colData(temp_cds)[[clustering_column]])
  write.table(clustering,"database/barcodetype.txt",sep = "\t")
  BarCluFile <- "database/barcodetype.txt"
  BarCluTable <- read.table(BarCluFile,sep = "\t",header = TRUE,stringsAsFactors = FALSE)
  
  RecClu <- unlist(diff_target_cell[next_age])
  LigClu <- setdiff(unlist(included), RecClu) 
  RecClus <- getHighExpGene(GCMat,BarCluTable,RecClu,LigClu,pval,logfc,cores)
  saveRDS(RecClus, paste(meta.file,"/",RecClu,"_", next_age,"_enriched.rds", sep=""))
  
}


# Top path analysis
for (age in age_list){
  
  temp_cds <- cds_base[, colData(cds_base)[[age_clustering_column]] == age]
  included <- in_groups[age]
  temp_cds <- temp_cds[, colData(temp_cds)[[clustering_column]] %in% unlist(included)]
  if (subset_by_time_point == F){ temp_cds <- cds_base[, colData(cds_base)[[clustering_column]] %in% unlist(included)]}

  colData(temp_cds)[[clustering_column]] <- droplevels(colData(temp_cds)[[clustering_column]])
  cell_group_df <- data.frame(cell=colnames(temp_cds) , group = colData(temp_cds)[[clustering_column]] ) 
  glom_normal <- aggregate_gene_expression(temp_cds,cell_group_df = cell_group_df, norm_method ="size_only" )
  row.names(glom_normal) <- rowData(temp_cds)$gene_short_name
  glom_normal <- glom_normal/as.vector(table(colData(temp_cds)[[clustering_column]]))
  glom_normal <- as.data.frame(glom_normal)

  
  glom_normal$genes <- row.names(glom_normal)
  glom_normal<- glom_normal[,c(ncol(glom_normal),1:(ncol(glom_normal)-1))]
  # min_exprs <- min(glom_normal[!glom_normal == 0])
  # thresh_exprs <- quantile(glom_normal[,-1][na.omit(!glom_normal[,-1] == 0)],0.1)   #1.652309e-07 for 122-13w
  
  
  min_exprs <- 1.713682e-11 # for 9-11w
  thresh_exprs <- .909563e-08
  
  CellChatDB<- readRDS("database/CellChatDB.rds")
  receptor_ligand <- receptor_ligand[receptor_ligand$Ligand.ApprovedSymbol %in% CellChatDB[ CellChatDB$pathway_name %in% path_included,]$ligand,]
  
  lr_glom_normal<-make_expressed_net(glom_normal,expressed_thresh= as.numeric(thresh_exprs),receptor_ligand,KL_method='product',pseudo_count= as.numeric(min_exprs))
  lr_glom_normal<-arrange(lr_glom_normal,desc(KL))
  head(lr_glom_normal)
  
  lr_glom_normal <- lr_glom_normal[lr_glom_normal$Pair.Evidence == "literature supported",]
  
  unique(lr_glom_normal$Ligand.ApprovedSymbol)[1:100]
  unique(lr_glom_normal$Receptor.ApprovedSymbol)[1:100]
  
  minum <- min(lr_glom_normal[,17])    #minum for12-13w 1.000009e-07
  lig_col <- 17:(17+ncol(glom_normal)-2)
  rec_col <- (17+ncol(glom_normal)-1):((17+ncol(glom_normal)-1)+ncol(glom_normal)-2)
  
  talklr::plot_lr_wiring(as.numeric(lr_glom_normal[5,lig_col]),as.numeric(lr_glom_normal[5,rec_col]),cell_labels = colnames(glom_normal)[-1],thresh= thresh_exprs)
  
  
  
  
  lig_mat <- lr_glom_normal[,lig_col]
  colnames(lig_mat) <- colnames(glom_normal)[-1]
  recived_in_sig <- c()
  edgelist_all <- matrix(, nrow = 0, ncol = 2)
  edgelist_rec_only <- matrix(, nrow = 0, ncol = 2)
  for (i in seq(row.names(lr_glom_normal))){
    net <- talklr::plot_lr_wiring(as.numeric(lr_glom_normal[i,lig_col]),as.numeric(lr_glom_normal[i,rec_col]),cell_labels = colnames(glom_normal)[-1],thresh= thresh_exprs)
    edgelist<-  as_edgelist(net)
    edgelist_all <- rbind(edgelist_all, edgelist)
    if (receiver_cell[age] %in% edgelist[,2]){
      recived_in_sig <- c(recived_in_sig,i)
      senders <- edgelist[which(edgelist[,2] %in% receiver_cell[age]),,drop =F][,1]
      lig_mat[i, !colnames(lig_mat) %in% senders] <- 0
      
      edgelist_rec_only <- rbind(edgelist_rec_only, edgelist)
    }
    
  }
  
  
  
  lr_glom_normal <- lr_glom_normal[recived_in_sig,]
  lig_mat <- lig_mat[recived_in_sig,]
  
  
  
  GCMat <- counts(temp_cds)
  row.names(GCMat)<- rowData(temp_cds)$gene_short_name
  colnames(GCMat) = as.character(colnames(temp_cds))
  rownames(GCMat) = as.character(rowData(temp_cds)$gene_short_name)
  
  class(rownames(GCMat))=="character"
  
  
  #GCMat<- as(GCMat,"dgCMatrix")
  #GCMat<- as.matrix(GCMat)
  
  # import sample annotation
  clustering <- data.frame(Barcode= colnames(temp_cds),Cluster = colData(temp_cds)[[clustering_column]])
  write.table(clustering,"database/barcodetype.txt",sep = "\t")
  BarCluFile <- "database/barcodetype.txt"
  BarCluTable <- read.table(BarCluFile,sep = "\t",header = TRUE,stringsAsFactors = FALSE)
  
  
  
  types <- unique(BarCluTable$Cluster)
  
  pair.list <- mapply(c, receiver_cell, diff_target_cell, SIMPLIFY = FALSE)
  next_age <- age_next_list[which(age_list == age)]
  enrich_1<- readRDS(paste(meta.file,"/",receiver_cell[age],"_", age,"_enriched.rds", sep=""))
  enrich_2<- readRDS(paste(meta.file,"/",diff_target_cell[next_age],"_", next_age,"_enriched.rds", sep=""))
  deg <- readRDS(paste(meta.file,"/","deg_0.05_",paste(unlist(pair.list[age]),collapse = '_'), ".rds",sep =''))
  filtered <- c(setdiff(c(subset(deg,deg$State%in%c("down","up"))$Gene),enrich_1),
                         intersect(intersect(enrich_1,enrich_2),subset(deg,deg$State%in%c("down","up"))$Gene))
  RecClu <- unlist(receiver_cell[age])    # the reciver
  LigClu <- setdiff(unlist(included), RecClu)     #all clusters except the recevier
  pval <- 0.05
  logfc <-0.15
  LigRecLib <- "database/LigRec.txt"
  TFTarLib <- "database/TFTargetGene.txt"  #"./database/RecTF.txt"
  RecTFLib <- "database/RecTF.txt"   
  
# saving metadata
  save(list=(c("lr_glom_normal","glom_normal","lig_mat","edgelist_rec_only","deg")), file = paste(meta.file,"/","talklr_",receiver_cell[age],".Rdata",sep=""))
  
  netList <- RunMLnet(GCMat, BarCluFile, RecClu, LigClu, 
                      pval, logfc, 
                      LigRecLib, TFTarLib, RecTFLib,
                      age)
  
  workdir <- paste(user.file,"metadata",age,sep="/")
  PyHome <- "python3" # not sure if it works yet python3 /home/ubuntu/miniconda3/bin/
  DrawMLnet(netList,LigClu,RecClu,workdir,PyHome,plotMLnet = T)
  
  saveRDS(netList, paste("top_path_output",workdir,"netList.rds",sep="/"))
  
  
    top_pathway(
      netList,
      deg,
      TF_fun = "mean",
      main_fun = "sum",
      path_fun = "sum",
      receiver_cell =  receiver_cell[age],
      lr_glom_normal = lr_glom_normal,
      show.sub.pathway = F,
      method = "NULL" , 
      palette = palette,
      output_dir= paste(out.file,user.file,sep= "/")
    )
}
