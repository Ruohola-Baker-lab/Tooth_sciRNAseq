setwd("~/data_sci_RNA_seq")
library(tidyr)
library(monocle3)
library(ggplot2)
library(dplyr)
library(Matrix.utils)


# library(devtools)
# install_github("jokergoo/ComplexHeatmap")

Amelo_cds_MTi <-readRDS("epi_cds_MTi_70_manhattan_8L_0.3_7k_pesudo_regressed_assigned_cell_type_7.rds")






# remove rowsum == 0 (aka genes with zero expression in all cells)
#Amelo_cds_MTi <- Amelo_cds_MTi[rowSums(counts(Amelo_cds_MTi)) !=0 ,]
levels(Amelo_cds_MTi$age_group)[levels(Amelo_cds_MTi$age_group)=="9_11w"] <- "09_11w"

# extract expression matrix
cell_group_df <- data.frame(cell=colnames(Amelo_cds_MTi) , group = Amelo_cds_MTi$assigned_cell_type.4 ) #clusters(Amelo_cds_MTi)
mat_all <- aggregate_gene_expression(Amelo_cds_MTi,cell_group_df = cell_group_df, norm_method ="size_only" )

# extract expression matrix per age group
cell_group_df <- data.frame(cell=colnames(Amelo_cds_MTi) , group = Amelo_cds_MTi$age_group ) #clusters(Amelo_cds_MTi)
mat_age <- aggregate_gene_expression(Amelo_cds_MTi,cell_group_df = cell_group_df, norm_method ="size_only" )
mat_age <- mat_age[,c("09_11w","12_13w","14_16w","17_19w","20_22w")]

# reorder the columns if necessary
# mat_all<- mat_all[,c("OE","DE-Prog","OEE","SI-1","SI-2","SI-3","EK","SR-1","SR-2","PA-1","PA-2","AM-1","AM-2")]

# extract a matrix of percetage of cells expressing each gene per cluster
percentage_mat <-Matrix::Matrix(0, nrow = nrow(Amelo_cds_MTi), ncol  = ncol(mat_all) )
colnames(percentage_mat) <- colnames(mat_all)
row.names(percentage_mat) <- row.names(mat_all)

for (i in colnames(mat_all)){
  non_zero<- Matrix::rowSums(counts(Amelo_cds_MTi[, Amelo_cds_MTi$assigned_cell_type.4 == i])!=0)
  all_sum <- Matrix::rowSums(counts(Amelo_cds_MTi[, Amelo_cds_MTi$assigned_cell_type.4 == i])!=0 | counts(Amelo_cds_MTi[, Amelo_cds_MTi$assigned_cell_type.4 == i]) ==0)
  percentage_mat[,which(i == colnames(mat_all))] <-   (non_zero/all_sum)*100
}



# Scale by row
mat_all <- t(scale(t(mat_all)))  #scale = apply(t(mat_all), 2, sd, na.rm = TRUE)
mat_all[is.nan(mat_all)] <- 0

mat_age <- t(scale(t(mat_age)))  #scale = apply(t(mat_age), 2, sd, na.rm = TRUE)
mat_age[is.nan(mat_age)] <- 0

# determine max
max_rows <-apply(X=mat_all, MARGIN=1, FUN=max)
#mat_all[mat_all ==max_rows] <-  0



# select the max genes and order them diagonal 
df <- c(
  sort(mat_all[mat_all[,1] >= max_rows, 1], decreasing = T),
  sort(mat_all[mat_all[,2] >= max_rows, 2], decreasing = T),
  sort(mat_all[mat_all[,3] >= max_rows, 3], decreasing = T),
  sort(mat_all[mat_all[,4] >= max_rows, 4], decreasing = T),
  sort(mat_all[mat_all[,5] >= max_rows, 5], decreasing = T),
  sort(mat_all[mat_all[,6] >= max_rows, 6], decreasing = T),
  sort(mat_all[mat_all[,7] >= max_rows, 7], decreasing = T),
  sort(mat_all[mat_all[,8] >= max_rows, 8], decreasing = T),
  sort(mat_all[mat_all[,9] >= max_rows, 9], decreasing = T),
  sort(mat_all[mat_all[,10] >= max_rows, 10], decreasing = T),
  sort(mat_all[mat_all[,11] >= max_rows, 11], decreasing = T),
  sort(mat_all[mat_all[,12] >= max_rows, 12], decreasing = T),
  sort(mat_all[mat_all[,13] >= max_rows, 13], decreasing = T)
)

# to keep track of which cluster has the heighest expression per particular gene
highest <-c(
  rep(colnames(mat_all)[1], length(mat_all[mat_all[,1] >= max_rows, 1])),
  rep(colnames(mat_all)[2], length(mat_all[mat_all[,2] >= max_rows, 2])),
  rep(colnames(mat_all)[3], length(mat_all[mat_all[,3] >= max_rows, 3])),
  rep(colnames(mat_all)[4], length(mat_all[mat_all[,4] >= max_rows, 4])),
  rep(colnames(mat_all)[5], length(mat_all[mat_all[,5] >= max_rows, 5])),
  rep(colnames(mat_all)[6], length(mat_all[mat_all[,6] >= max_rows, 6])),
  rep(colnames(mat_all)[7], length(mat_all[mat_all[,7] >= max_rows, 7])),
  rep(colnames(mat_all)[8], length(mat_all[mat_all[,8] >= max_rows, 8])),
  rep(colnames(mat_all)[9], length(mat_all[mat_all[,9] >= max_rows, 9])),
  rep(colnames(mat_all)[10], length(mat_all[mat_all[,10] >= max_rows, 10])),
  rep(colnames(mat_all)[11], length(mat_all[mat_all[,11] >= max_rows, 11])),
  rep(colnames(mat_all)[12], length(mat_all[mat_all[,12] >= max_rows, 12])),
  rep(colnames(mat_all)[13], length(mat_all[mat_all[,13] >= max_rows, 13]))
)



row.names(mat_all) = rowData(Amelo_cds_MTi[row.names(mat_all),])$gene_short_name

percentage_mat <-  percentage_mat[names(df),]

mat_age <- mat_age[names(df),]

ensemble_df <-as.data.frame(names(df))
names(df) = rowData(Amelo_cds_MTi[names(df),])$gene_short_name
row.names(mat_age) = rowData(Amelo_cds_MTi[row.names(mat_age),])$gene_short_name




# create dataframe of the ordered genes
df <- data.frame(id= names(df), var = df)

# label rows by serial numbers for easier slicing
row.names(df) <- 1:nrow(df)
names(highest) <- 1:nrow(df)
row.names(percentage_mat) <- 1:nrow(df)
row.names(mat_age) <- 1:nrow(df)


## using unique  to exclude duplicate maxima
df <- unique(subset(df, var==ave(var, id, FUN=max)) )


# slice out duplicate from highest list 
highest <- highest[row.names(df)]
ensemble_df <- ensemble_df[row.names(df),]

# slice out duplicate from percentage_mat 
percentage_mat <- percentage_mat[row.names(df),]
row.names(percentage_mat) = df$id

mat_age <- mat_age[row.names(df),]
row.names(mat_age) = df$id

# add it to the dataframe
df <- cbind(df,highest )
df <- cbind(df,ensemble_df )

# rename rownames back
row.names(df) <- df$id
#mat_all <- mat_all[df$id,]

cumulative_index <- cumsum(table(df$highest)[colnames(mat_all)])

cell_percentage <- c(percentage_mat[1:cumulative_index[1],1],
                     percentage_mat[(cumulative_index[1]+1):cumulative_index[2],2],
                     percentage_mat[(cumulative_index[2]+1):cumulative_index[3],3],
                     percentage_mat[(cumulative_index[3]+1):cumulative_index[4],4],
                     percentage_mat[(cumulative_index[4]+1):cumulative_index[5],5],
                     percentage_mat[(cumulative_index[5]+1):cumulative_index[6],6],
                     percentage_mat[(cumulative_index[6]+1):cumulative_index[7],7],
                     percentage_mat[(cumulative_index[7]+1):cumulative_index[8],8],
                     percentage_mat[(cumulative_index[8]+1):cumulative_index[9],9],
                     percentage_mat[(cumulative_index[9]+1):cumulative_index[10],10],
                     percentage_mat[(cumulative_index[10]+1):cumulative_index[11],11],
                     percentage_mat[(cumulative_index[11]+1):cumulative_index[12],12],
                     percentage_mat[(cumulative_index[12]+1):cumulative_index[13],13])


# add it to the dataframe
df <- cbind(df, cell_percentage)

df$major_age <- factor(colnames(mat_age)[apply(mat_age,1,which.max)],levels = c("09_11w","12_13w","14_16w","17_19w","20_22w"))

# select genes of interest to display ( order doesn't matter here)
c_1 <- c("KRT13","PARD3", "ANXA1","EGFR","EXOC6B")  # KRT36, PARD3, ANXA1, ,"FAF1","CFTR"
c_2 <- c( "MCHR2","TBCE","DSG2","PITX2","TBX15")# "CLDN34","TBCE", #,"KRT5",krt14, PITX2, "MCHR2","TNFSF10", TBX15
c_3 <-  c("FGF5","SLC14A2","PTN")# TGFBR3L,PTN "CAVIN3",
c_4 <- c( "ROS1","NOTCH1","HOXC4", "MDK", c("ARX","LBX2","EFNA5","ST5"), c("RAET1E","C1QB","DST","ATF3"))# TEAD3
c_5 <- c("FGF4","CUX2","ADGRE1", "MCC")#"ERBB4")"NTRK2" 
c_6 <- c("PEAK1","APP","HAS2","CLDN8","BICC1","CPLX4","SOX5","ZEB2","EMX1")#,"ANTXR1", "SMOC2","RUNX2","TMEM132C","HDAC9","PRRX1","LSAMP","SLIT3","CDH12","HAPLN1","SOX5","DLC1","CLMP","RSPO4","TPST1","TWIST2","GLIS1","ROBO1","VCAN","NRXN1","PTCH1") #,"VAT1L","COL22A1","IBSP","VWDE","C1QB","GNRH1","KIF5C","VWA2","IGSF9","GDF5","CHST1","CKB","STAC3","CCN4","SOX21","GAD1","VIT",'RAB3IL1',"PMCH","CYP2S1") 
c_7 <- c("FGF6","DACH1","SP6","SHH","PAQR3","INSL5","EML6","PATJ") #WNT5A,,"RYR2","RYR2","NXN",
c_8 <- c("DSPP","AMELX","AMBN","ENAM","MMP20","SDK1","FGFR1","TJP1")

# # Here's to group them togther in one list
selected_markers <- c(c_1,c_2,c_3,c_4,c_5,c_6,c_7,c_8)#,c_9,c_10)




# only display the selected labels
df_lab <- df[,c(1,2)]
df_lab$id <- paste("-",df_lab$id)
df_lab[!row.names(df_lab) %in% selected_markers,]<-""
#mat_all[ !duplicated(row.names(mat_all) ),]
setdiff(selected_markers,row.names(df_lab))

df$percentage_of_cells <- cut(df$cell_percentage, 
                              breaks=c(0, 20, 100), 
                              labels=c("<20%",">20%"))

ann_colors = list(Clusters = setNames(
  c(
    "#c8003d",
    "#8e371b",
    "#f6bb77",
    "#5b5f00",
    "#32d058",
    "#01a65b",
    "#92d4b0",
    "#00b0e9",
    "#025bc0",
    "#ff9bff",
    "#ff65d2",
    "#823966",
    "#ec0e84"
  ),
  c(colnames(mat_all)
  )
))
# 
# df$highest <- factor(df$highest, levels = c("OE","DE-Prog","OEE","SI-1","SI-2","SI-3","EK","SR-1","SR-2","PA-1","PA-2","AM-1","AM-2") )
# df_sub <- as.data.frame(c("OE","DE-Prog","OEE","SI-1","SI-2","SI-3","EK","SR-1","SR-2","PA-1","PA-2","AM-1","AM-2"))
# 

# df_sub <- as.data.frame(cut(
#   df$cell_percentage,
#   breaks = c(0, 20, 100),
#   labels = c("<20%", ">20%")
# ))

df_sub <- as.data.frame(colnames(mat_all))
colnames(df_sub) <- "Clusters"
row.names(df_sub) <- df_sub$Clusters
#colnames(df_sub) <- "Clusters"


# add annotation for age group
age_group_df <- table(colData(Amelo_cds_MTi)[,c(9,24)])
age_group_sum <- rowSums(age_group_df)
age_group_df <- (age_group_df/age_group_sum)*100
df_sub$major_age <-row.names(age_group_df)[apply(age_group_df,2,which.max)]
df_sub$major_age <- factor(df_sub$major_age, levels = c("09_11w","12_13w","14_16w","17_19w","20_22w"))


#age_group_df <- table(colData(Amelo_cds_MTi)[,c(9,24)])
age_group_df <- t(t(age_group_df)/colSums(age_group_df))
df_sub$age_score <- sort(scales::rescale(colSums((age_group_df) * c(1:5)), to = c(0,100)))



age_df <- as.data.frame(setNames(df$major_age, df$id))
colnames(age_df) <- "major_age"


mat_all <- mat_all[row.names(df), ]

#plotting

library(ComplexHeatmap)

pdf("heatmap.pdf")
ComplexHeatmap::Heatmap(
  name = "expression",
  mat_all ,
  cluster_columns = F,
  cluster_rows = F,
  show_row_name = F,
  show_column_names = T,
  column_names_side = "top",
  column_names_rot = 45,
  bottom_annotation  = HeatmapAnnotation(age_score =df_sub[, 3] , age_group=t(age_group_df), annotation_name_side="right"), #major_age= df_sub[, 2],
  left_annotation  = rowAnnotation(foo =anno_mark(at = which(row.names(mat_all) %in% selected_markers ), labels =row.names(df)[row.names(df) %in% selected_markers ], side = "left", labels_gp = gpar(fontsize  = 5))),
  row_split  = factor(df$highest, levels = colnames(mat_all)),
  row_title = " ",
  use_raster = T,
  show_heatmap_legend = T
)

dev.off()



library(circlize)
col_fun = colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))
lgd = Legend(col_fun = col_fun, title = "expression", direction = "vertical")
draw(lgd, x = unit(1, "npc"), y = unit(1, "npc"), just = c("right", "top"))



ComplexHeatmap::pheatmap(
  t(mat_all[df$id, ]) ,
  scale = "column",
  cluster_cols = F,
  cluster_rows = F,
  show_rownames = F,
  show_colnames = T,
  border_color = "grey",
  labels_col = df_lab$id,
  fontsize_col = 8,
  annotation_row = df_sub[,c(2,3)] ,
  annotation_colors = ann_colors ,
  annotation_names_col = F,
  gaps_col = cumsum(table(df$highest)[colnames(mat_all)]),

)#, gaps_col = cumsum(c(1,1,1,3,1,2,2,2)))




ComplexHeatmap::pheatmap(
  mat_all ,
  scale = "row",
  cluster_cols = F,
  cluster_rows = F,
  show_rownames = T,
  show_colnames = T,
  border_color = NA,
  labels_row = df_lab$id,
  fontsize_row = 7,
  annotation_col = df_sub[, c(2, 3)] ,
  annotation_colors = ann_colors ,
  annotation_names_col = F,
  annotation_names_row = F,
  annotation_row = age_df,
  gaps_row = cumsum(table(df$highest)[colnames(mat_all)])
)#, gaps_col = cumsum(c(1,1,1,3,1,2,2,2)))




#  go term analysis
#BiocManager::install("topGO")
#BiocManager::install("ViSEAGO")
remotes::install_github("hadley/dplyr@v0.8.2")
remotes::install_github("r-lib/gert")
install.packages("devtools")
devtools::install_github("hadley/devtools")
devtools::install_url("https://cran.r-project.org/src/contrib/Archive/dbplyr/dbplyr_1.3.0.tar.gz")

deb: libgit2-dev

library(topGO)
library(ViSEAGO)





# load genes background
library(org.Hs.eg.db)
background=keys(org.Hs.eg.db, keytype ='ENSEMBL')
background=intersect(background,row.names(Amelo_cds_MTi))
background <- background[background %in% row.names(Amelo_cds_MTi)]

# load gene selection
selection<- df[df$highest== "AM-2",]
selection<- top_n(selection,500,cell_percentage)

selection<- background[background %in% selection$ensemble_df]
setdiff(selection,background) # make sure this is zero



save(list= c("background","selection", "df", "mat_all","myGENE2GO","BP"), file = "basic_heatmap.Rdata")
load("basic_heatmap.Rdata")

# connect to Ensembl
Ensembl<-ViSEAGO::Ensembl2GO()
# Display table of available organisms with Ensembl
ViSEAGO::available_organisms(Ensembl)
# load GO annotations from Ensembl
myGENE2GO<-ViSEAGO::annotate(
  "hsapiens_gene_ensembl",
  Ensembl
)



# create topGOdata for BP
BP<-ViSEAGO::create_topGOdata(
  geneSel=selection,
  allGenes=background,
  gene2GO=myGENE2GO, 
  ont="BP",
  nodeSize=5
)


# perform TopGO test using clasic algorithm
classic<-topGO::runTest(
  BP,
  algorithm ="classic",
  statistic = "fisher"
)


# merge results from topGO
BP_sResults<-ViSEAGO::merge_enrich_terms(
  Input=list(
    condition=c("BP","classic")
  )
)


# display the merged table
ViSEAGO::show_table(BP_sResults)

# print the merged table in a file
ViSEAGO::show_table(
  BP_sResults,
  "BP_sResults.xls"
)


# initialyse 
myGOs<-ViSEAGO::build_GO_SS(
  gene2GO=myGENE2GO,
  enrich_GO_terms=BP_sResults
)

# compute all available Semantic Similarity (SS) measures
myGOs<-ViSEAGO::compute_SS_distances(
  myGOs,
  distance="Wang"
)


# display MDSplot
ViSEAGO::MDSplot(myGOs)

# print MDSplot
ViSEAGO::MDSplot(
  myGOs,
  file="mdsplot1.png"
)

# GOterms heatmap with the default parameters
Wang_clusters_wardD2<-ViSEAGO::GOterms_heatmap(
  myGOs,
  showIC=TRUE,
  showGOlabels=TRUE,
  GO.tree=list(
    tree=list(
      distance="Wang",
      aggreg.method="ward.D2"
    ),
    cut=list(
      dynamic=list(
        pamStage=TRUE,
        pamRespectsDendro=TRUE,
        deepSplit=2,
        minClusterSize =2
      )
    )
  ),
  samples.tree=NULL
)

# Display the clusters-heatmap
ViSEAGO::show_heatmap(
  Wang_clusters_wardD2,
  "GOterms"
)

# print the clusters-heatmap
ViSEAGO::show_heatmap(
  Wang_clusters_wardD2,
  "GOterms",
  "cluster_heatmap_Wang_wardD2.png"
)




################

OE_Am_cds <- readRDS ("data_sci_RNA_seq/OE_Am_cds_may22.rds")

Amelo_cds_MTi <-readRDS("epi_cds_MTi_70_manhattan_8L_0.3_7k_pesudo_regressed_assigned_cell_type_7.rds")


colData(Amelo_cds_MTi)$assigned_cell_type.simp <- dplyr::recode_factor(colData(Amelo_cds_MTi)$assigned_cell_type.4,
                                                              'OE'='OE',
                                                              'DE-Prog'='DE-Prog',
                                                              'EK'='IK/EK',
                                                              'OEE'='OEE',
                                                              'SR-1'="SR-1",
                                                              'SR-2'="SR-2",
                                                              'SI-1'='SI-1',
                                                              'SI-2'='SI-2',
                                                              'SI-3'='SI-3',
                                                              'PA-1' = 'PA-1',
                                                              'PA-2' = 'PA-2',
                                                              'AM-1' = 'AM-1',
                                                              'AM-2' = 'AM-2')

colData(Amelo_cds_MTi)$age_group<- dplyr::recode_factor(colData(Amelo_cds_MTi)$age_group,
                                                        '9_11w'='09_11gw', '12_13w'='12_13gw', '14_16w'='14_16gw', '17_19w'='17_19gw', '20_22w'='20_22gw')


amelo_trajectory <- Amelo_cds_MTi[,colData(Amelo_cds_MTi)$assigned_cell_type.4 %in% c("OE", "DE-Prog", "OEE", "PA-1", "PA-2", "AM-1", "AM-2")]




Amelo_cds_MTi_incisor <- amelo_trajectory[, !colData(amelo_trajectory)$sample_id %in% c("molars_12","molars_17","molars_20")]
Amelo_cds_MTi_molars <- amelo_trajectory[, !colData(amelo_trajectory)$sample_id %in% c("incisors_12","incisors_17","incisors_20")]

colData(Amelo_cds_MTi_incisor)$amelo_trajectory <- dplyr::recode_factor(colData(Amelo_cds_MTi_incisor)$assigned_cell_type.4,
                                                                 "OE" = "OE_incisor", "DE-Prog" = "DE-Prog_incisor", "OEE"="OEE_incisor", "PA-1"="PA-1_incisor", "PA-2"="PA-2_incisor", "AM-1"="AM-1_incisor", "AM-2" = "AM-2_incisor")


colData(Amelo_cds_MTi_molars)$amelo_trajectory <- dplyr::recode_factor(colData(Amelo_cds_MTi_molars)$assigned_cell_type.4,
                                                                "OE" = "OE_molar", "DE-Prog" = "DE-Prog_molar", "OEE"="OEE_molar", "PA-1"="PA-1_molar", "PA-2"="PA-2_molar", "AM-1"="AM-1_molar", "AM-2" = "AM-2_molar")


Amelo_cds_MTi_incisor_molars <- combine_cds(list(Amelo_cds_MTi_incisor,Amelo_cds_MTi_molars))


colData(Amelo_cds_MTi_incisor_molars)$amelo_trajectory <- droplevels(colData(Amelo_cds_MTi_incisor_molars)$amelo_trajectory)

source("/home/ubuntu/data_sci_RNA_seq/custom_Heatmap.R")
custom_Heatmap(
  cds = Amelo_cds_MTi,
  clustering_column = "assigned_cell_type.simp",
  age_clustering_column = "age_group",
  selected_markers = c("auto"),
  number_of_markers = 3, 
  number_of_genes2goterm = 50,
  gene_font_size= 10,
  keywords_font_size= 13, 
  number_of_keywords = 7,
  keywords_blacklist = c(
    "symbiont",
    "head",
    "cardiovascular",
    "transcript",
    "regulated", 
    "cartilage",
    "generation",
    "response",
    "chondrocyte",
    "platelet",
    "sensory",
    "neurons",
    "osteoblast",
    "bone",
    "ureteric",
    "urogenital",
    "renal",
    "relaxation",
    "gland",
    "blood",
    "quality",
    "smooth",
    "ossification",
    "odontogenesis",
    "kidney",
    "head",
    "septum",
    "aortic",
    "taxis",
    "endogenous",
    "valve",
    "locomotion",
    "subcellular",
    "movement",
    "axon",
    "luteolysis",
    "chemical",
    "stimulus",
    "mesonephric",
    "substance",
    "forebrain",
    "homotypic",
    "embryonic",
    "endodermal",
    "drug",
    "nephron",
    "renal",
    "tubule",
    "metanephros",
    "free",
    "ubiquitin",
    "chain",
    "poly",
    "polymerization",
    "synaptic",
    "transmission",
    "spontaneous",
    "neurotransmitter",
    "synapse",
    "vessel",
    "lineage",
    "unfolded",
    "inorganic",
    "estradiol",
    "ketone",
    "lipid",
    "endochondral",
    "topologically",
    "posttranscriptional",
    "cortisol",
    "reverse",
    "signal",
    "cholesterol",
    "interleukin",
    "endoplasmic",
    "incorrect",
    "axonogenesis",
    "sterol",
    "motor",
    "telomere",
    "peptide",
    "folding",
    "anion",
    "reticulum",
    "recognition",
    "programmed",
    "death",
    "guidance",
    "hydroxylation",
    "plasma",
    "transduction",
    "cycle",
    "innervation",
    "organelle",
    "presynaptic",
    "vasculature",
    "heterophilic",
    "chemotaxis",
    "tripeptide",
    "assembly",
    "odontoblast",
    "telencephalon",
    "amyloid",
    "fibril",
    "structure",
    "oligopeptide",
    "artery",
    "pallium",
    "aorta",
    "cerebral",
    "cortex",
    "axonal",
    "fasciculation"
    
    
  ),
  auto_sort_columns = F,
  load_previous = T,
  normalize_by_cluster = F

)



# diff dataset
custom_Heatmap(
  cds = OE_Am_cds,
  clustering_column = "clusters_named",
  age_clustering_column = "diff_time_point",
  selected_markers = c("auto","SOX2", "AMBN","BARX2","ACTA2"),
  number_of_markers = 2, 
  number_of_genes2goterm = 50,
  gene_font_size= 10,
  keywords_font_size= 13, 
  number_of_keywords = 7,
  keywords_blacklist = c(
    "symbiont",
    "transcript",
    "regulated",
    "generation",
    "substrate",
    "response",
    "blood",
    "ion",
    "chromosome",
    "poly",
    "nucleosome",
    "break",
    "splicing",
    "cycle",
    "acid",
    "acid",
    "canonical",
    "peptidase",
    "chondrocyte",
    "chromatin",
    "tail",
    "disassembly",
    "calcium",
    "spliceosome",
    "quality",
    "metanephros",
    "splice",
    "oxygen",
    "axon",
    "organelle",
    "plasma",
    "hydrolase",
    "epidermis",
    "skin",
    "coronary",
    "ossification",
    "metanephros",
    "centrosome",
    "guidance",
    "catalytic",
    "endogenous",
    "neurons",
    "cornification",
    "chemotaxis",
    "structure",
    "vasculature",
    "chordate",
    "stimulus",
    "transesterification",
    "interaction",
    "taxis",
    "programmed",
    "cation",
    "methylation",
    "substance",
    "complex",
    "primary",
    "reactions",
    "ear",
    "death",
    "outflow",
    "activation",
    "proteolysis",
    "bulged",
    "osteoclast",
    "ventricle",
    "circulation",
    "end",
    "tract",
    "ureteric",
    "rate",
    "transmembrane",
    "synaptic",
    "keratinization",
    "axonogenesis",
    "septum",
    "membrane",
    "adenosine",
    "damage",
    "extension",
    "junction",
    "synapse",
    "synapse",
    "nucleophile",
    "raft",
    "bounded",
    "assembly",
    "phase",
    "shortening",
    "repair",
    "remodeling",
    "heterochromatin",
    "extrinsic",
    "dendrite",
    "walking",
    "projection",
    "intermediate",
    "apoptotic",
    "homologous",
    "recombination",
    "NA",
    "filament",
    "cytoskeleton",
    "recombinational",
    "epidermal",
    "adult",
    "exocytosis",
    "extent",
    "biogenesis"

  ),
  auto_sort_columns = T,
  load_previous = T,
  normalize_by_cluster = F,
  disable_legends = T
)



selected_mat <- counts(Amelo_cds_MTi[rowData(Amelo_cds_MTi)$gene_short_name %in% Amelo_cds_MTi_selected_markers,])
row.names(selected_mat)<- rowData(Amelo_cds_MTi[row.names(selected_mat),])$gene_short_name

dput(names(sort(rowSums(selected_mat))))




age= "09_11w"
age= "12_13w"
age= "14_16w"
age= "17_19w"
age= "20_22w"

in_groups <- list('09_11w' = c("OE","IK/EK","DE-Prog") ,
                  '12_13w' = c("DE-Prog","IK/EK","SR-1") ,
                  '14_16w' = c("OEE" , "IK/EK",  "SI-1", "SI-2" ),
                  '17_19w' = c("PA-1","SI-2","SR-2") ,
                  '20_22w' = c("AM-1", "SI-3", "AM-2","PA-2"))

included <- in_groups[age]

temp_cds <- Amelo_cds_MTi[, colData(Amelo_cds_MTi)$assigned_cell_type.simp %in% unlist(included)]
colData(temp_cds)$assigned_cell_type.simp <- droplevels(colData(temp_cds)$assigned_cell_type.simp)

temp_cds <- Amelo_cds_MTi[, colData(Amelo_cds_MTi)$age_group == age]
colData(temp_cds)$age_group <- droplevels(colData(temp_cds)$age_group)


source("/home/ubuntu/functions.R")
plot_heatmap_over_realtime(
  cds= temp_cds,
  clustering_column = "assigned_cell_type.simp",
  age_clustering_column = "age_group",
  selected_genes = c("VCAN","PITX2","NOTCH1","DKK3","IRX1","SFRP5"),
  normalize = T,   # Toggle normalization
  trimm_low = T   # Trimm low cells that can be considered background and missclustered. In another words remove expression in time points that arenot likely to be real
)
