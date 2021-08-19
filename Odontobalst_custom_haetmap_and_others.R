
# load the customized heatmap function.
source("/home/ubuntu/data_sci_RNA_seq/custom_Heatmap.R")

# load your cds
DEendo_122920_15dim_Euclidean_11L_0.05dist_25k <- readRDS("DEendo_122920_15dim_Euclidean_11L_0.05dist_25k.rds")
newclusters_regressed_DE_121820_30dim_Euclidean_11L_0.05dist_25k <- readRDS("newclusters_regressed_DE_121820_30dim_Euclidean_11L_0.05dist_25k.rds")

# custom_Heatmap: only supply your cds,  a clustering_column & age_group_column
# you can specify the markers genes by selected_markers = c("DSPP", "AMBN")
# or leave it to automaticly find markers by selected_markers = "auto", then set the number_of_markers
# or do both selected_markers = c("auto", "DSPP","NES","DMP1","FGF3") to get top genes + yoursupplied genes. Note supllied geneswill be includedin goterm analysis
# number_of_genes2goterm define the number of top genes to  use for goterm analysis. You have to find the sweet spot number. not too many and you dliute the result or too little and you miss important things
# gene_font_size: size of the labeled genes
# keywords_font_size: size of goterm key words
# number_of_keywords: number of top keywords to display.
# keywords_blacklist: this manually remove redundant, meaningless or irrleveant words from the goterm keywords
# load_previous is to load all the saved metadata if present from previous run on same cds, so the second run on the same dataset can be faster. if you made alot of changes then select FALSE and calculate everyting again. Some variables won't do any effect if you cange them unless ou do a complete rerun
# to be more spcific canges to number_of_genes2goterm variable will update only if yoy not load
# all results are saved in custome_Heatmap_output folder
# it will tell you how many minutes it took to run at the end. This data take about 10min in first run

custom_Heatmap(
  cds = newclusters_regressed_DE_121820_30dim_Euclidean_11L_0.05dist_25k,
  clustering_column = "new_clusters",
  age_clustering_column = "Age_group",
  # you can specify the markers genes by selected_markers = c("DSPP", "AMBN")
  #selected_markers = c("DSPP","SP6", "AMBN", "AMELX"),
  # or leave it to automaticly find markers by selected_markers = "auto", then set the number_of_markers
  # or do both selected_markers = c("auto", "DSPP","NES","DMP1","FGF3") to get top genes + yoursupplied genes. Note supllied geneswill be includedin goterm analysis
  selected_markers = c("auto", "DSPP","SP6", "AMBN", "AMELX"),
  number_of_markers = 5, 
  number_of_genes2goterm = 50,
  gene_font_size= 10,
  keywords_font_size= 16, 
  number_of_keywords = 10,
  keywords_blacklist= c("symbiont","head","cardiovascular","transcript"),
  auto_sort_columns = T,
  load_previous = F
)



### REALTIME HEATMAP#####

source("/home/ubuntu/functions.R")

# before combining dental epi and mesenchyme
# make sure the clustering columns and age group column has the same spelling before you combine

amelo_trajectory_cds <- readRDS("/home/ubuntu/data_sci_RNA_seq/amelo_trajectory.rds")
amelo_trajectory_cds <- amelo_trajectory_cds[, colData(amelo_trajectory_cds)$assigned_cell_type.4 %in% c("PA-1","PA-2","AM-1","AM-2")]
colData(amelo_trajectory_cds)$new_clusters <- droplevels(colData(amelo_trajectory_cds)$assigned_cell_type.4)
colData(amelo_trajectory_cds)$Age_group <- droplevels(colData(amelo_trajectory_cds)$age_group)

odonto_trajectory_cds <- readRDS("/home/ubuntu/newclusters_regressed_DE_121820_30dim_Euclidean_11L_0.05dist_25k.rds")
colData(odonto_trajectory_cds)$new_clusters <- dplyr::recode_factor(colData(odonto_trajectory_cds)$new_clusters,
                                                                  "1" = "DE", "2" = "preOB", "3"="subOB_prog", "4"="subOB", "5"="OB")


odonto_trajectory_cds_sub <- odonto_trajectory_cds[, colData(odonto_trajectory_cds)$new_clusters %in% c('DE',"preOB","subOB","OB")]


amelo_odonto_cds <- combine_cds(list(amelo_trajectory_cds,odonto_trajectory_cds_sub))

colData(amelo_odonto_cds)$new_clusters <- droplevels(colData(amelo_odonto_cds)$new_clusters)

#colData(amelo_odonto_cds[,colData(amelo_odonto_cds)$Age_group =="9_11w"])$Age_group <- "09_11w"
colData(amelo_odonto_cds)$Age_group <- factor(colData(amelo_odonto_cds)$Age_group, levels= unique(colData(amelo_odonto_cds)$Age_group))
levels(colData(amelo_odonto_cds)$Age_group)[levels(colData(amelo_odonto_cds)$Age_group)=="9_11w"] <- "09_11w"



source("/home/ubuntu/functions.R")
plot_heatmap_over_realtime(
  cds= amelo_odonto_cds,
  clustering_column = "new_clusters",
  age_clustering_column = "Age_group",
  selected_genes = c("SP6", "DSPP", "AMBN", "AMELX"),
  normalize = T,   # Toggle normalization
  trimm_low = T   # Trimm low cells that can be considered background and missclustered. In another words remove expression in time points that arenot likely to be real
)




#### Incisor vs Molar ####

#odonto_trajectory_cds <- odonto_trajectory_cds[,colData(odonto_trajectory_cds)$new_clusters %in% c(3,4,5)]


odonto_trajectory_cds_incisor <- odonto_trajectory_cds[, colData(odonto_trajectory_cds)$Sample_id %in% c("New_Anterior_9_11w", "Incisors_12_13w", "Incisors_Molars_14_16w", "Incisors_17_19w", "Incisors_20_22w")]
odonto_trajectory_cds_molars <- odonto_trajectory_cds[, colData(odonto_trajectory_cds)$Sample_id %in% c("Posterior_9_11w" , "Molars_12_13w","Incisors_Molars_14_16w","Molars_17_19w", "Molars_20_22w" )]



colData(odonto_trajectory_cds_incisor)$combined_trajectory <- dplyr::recode_factor(colData(odonto_trajectory_cds_incisor)$new_clusters,
                                                                                   "DE" = "DE_incisor", "preOB" = "preOB_incisor", "subOB_prog"="subOB_prog_incisor", "subOB"="subOB_incisor", "OB"="OB_incisor")


colData(odonto_trajectory_cds_molars)$combined_trajectory <- dplyr::recode_factor(colData(odonto_trajectory_cds_molars)$new_clusters,
                                                                                  "DE" = "DE_molar", "preOB" = "preOB_molar", "subOB_prog"="subOB_prog_molar", "subOB"="subOB_molar", "OB"="OB_molar")


odonto_trajectory_cds_incisor_molars <- combine_cds(list(odonto_trajectory_cds_incisor,odonto_trajectory_cds_molars))

colData(odonto_trajectory_cds_incisor_molars)$combined_trajectory <- droplevels(colData(odonto_trajectory_cds_incisor_molars)$combined_trajectory)
colData(odonto_trajectory_cds_incisor_molars)$new_clusters <- droplevels(colData(odonto_trajectory_cds_incisor_molars)$new_clusters)

#colData(odonto_trajectory_cds_incisor_molars[,colData(odonto_trajectory_cds_incisor_molars)$Age_group =="9_11w"])$Age_group <- "09_11w"
colData(odonto_trajectory_cds_incisor_molars)$Age_group <- factor(colData(odonto_trajectory_cds_incisor_molars)$Age_group, levels= unique(colData(odonto_trajectory_cds_incisor_molars)$Age_group))
levels(colData(odonto_trajectory_cds_incisor_molars)$Age_group)[levels(colData(odonto_trajectory_cds_incisor_molars)$Age_group)=="9_11w"] <- "09_11w"





######### bar graph
full_cluster_sum <- table(colData(odonto_trajectory_cds_incisor_molars)[,c("sample","new_clusters")])
sample_id_cluster_sum <- table(colData(odonto_trajectory_cds_incisor_molars)[,c("Sample_id","new_clusters")])
sample_id_sum<- rowSums(sample_id_cluster_sum)
full_cluster_sum <- (sample_id_cluster_sum/sample_id_sum)
subset_full_cluster_sum <- full_cluster_sum[c("Incisors_Molars_14_16w" ),drop=F,]
full_cluster_sum <- rbind(full_cluster_sum,subset_full_cluster_sum)
full_cluster_sum <- full_cluster_sum[c("New_Anterior_9_11w","Incisors_12_13w","Incisors_Molars_14_16w", "Incisors_17_19w", "Incisors_20_22w",  
                                       "Posterior_9_11w","Molars_12_13w","Incisors_Molars_14_16w","Molars_17_19w", "Molars_20_22w"),]  #sort them in proper order where incisors on top
full_cluster_sum <- rbind(full_cluster_sum,colSums((full_cluster_sum[1:5,])))
full_cluster_sum <- rbind(full_cluster_sum,colSums((full_cluster_sum[6:10,])))
full_cluster_sum <- full_cluster_sum[11:12,]
dimnames(full_cluster_sum)<- list(tooth_type= c("Incisor","Molar"),trajectory =colnames(full_cluster_sum))
reshaped <- reshape2::melt(full_cluster_sum, direction = "long")
reshaped$trajectory <- factor(reshaped$trajectory, rev(levels(reshaped$trajectory)))
ggplot(data=(reshaped), aes(y= value, x = trajectory ,fill= tooth_type )) + 
  geom_col(show.legend = T, position = "fill")+theme(text=element_text(size=14, face="bold"), axis.title.y = element_blank()) +theme(legend.position="top",legend.title = element_blank())+coord_flip()+geom_hline(yintercept=0.5, linetype="longdash")+
  scale_y_reverse(labels=c("100%","75%","50%","25%","0%"), name = 'Proportion of cells\n (after normalization)')



### When is the first odonoblats appear ####
full_cluster_sum <- table(colData(odonto_trajectory_cds_incisor_molars)[,c("sample","new_clusters")])
sample_id_cluster_sum <- table(colData(odonto_trajectory_cds_incisor_molars)[,c("Sample_id","new_clusters")])
sample_id_sum<- rowSums(sample_id_cluster_sum)
full_cluster_sum <- (sample_id_cluster_sum/sample_id_sum)
subset_full_cluster_sum <- full_cluster_sum[c("Incisors_Molars_14_16w" ),drop=F,]
full_cluster_sum <- rbind(full_cluster_sum,subset_full_cluster_sum)
full_cluster_sum <- full_cluster_sum[c("New_Anterior_9_11w","Incisors_12_13w","Incisors_Molars_14_16w", "Incisors_17_19w", "Incisors_20_22w",  
                                       "Posterior_9_11w","Molars_12_13w","Incisors_Molars_14_16w","Molars_17_19w", "Molars_20_22w"),]  #sort them in proper order where incisors on top


#full_cluster_sum <- full_cluster_sum[c(1:5,9,6,10,7:8),]

incedence_I <- list()
for (d in seq(colnames(full_cluster_sum))){
  element <- (full_cluster_sum[1:5,d][full_cluster_sum[1:5,d] > sum(full_cluster_sum[1:5,d]) * 0.05]) #quantile(full_cluster_sum[1:5,d], probs = 0.5)])[1]    E# mean(full_cluster_sum[1:5,d]) * 0.5]
  name <- colnames(full_cluster_sum)[d] 
  incedence_I[[name]] <- element
  
}
incedence_I


incedence_M <- list()
for (d in seq(colnames(full_cluster_sum))){
  element <- (full_cluster_sum[6:10,d][full_cluster_sum[6:10,d] >= sum(full_cluster_sum[1:5,d]) * 0.05])
  name <- colnames(full_cluster_sum)[d] 
  incedence_M[[name]] <- element
  
}

incedence_I <- incedence_I[c("OB")]   
incedence_M <- incedence_M[c("OB")]
c(incedence_I,incedence_M)


incedence_tib <- tibble::enframe(unlist(list(incedence_I,incedence_M)))
incedence_tib$Type <- ""
incedence_tib[grep("Incisors", incedence_tib$name ),]$Type <- "Incisor"
incedence_tib[grep("Anterior", incedence_tib$name ),]$Type <- "Incisor"
incedence_tib[grep("Molar", incedence_tib$name ),]$Type <- "Molar"
incedence_tib[grep("Posterior_9", incedence_tib$name ),]$Type <- "Molar"
incedence_tib[grep("Incisors_Molars", incedence_tib$name ),]$Type[1] <- "Incisor"

incedence_tib$Age_group <- ""
incedence_tib[grep("09_11w", incedence_tib$name ),]$Age_group <- "09_11w"
incedence_tib[grep("12_13w", incedence_tib$name ),]$Age_group <- "12_13w"
incedence_tib[grep("14_16w", incedence_tib$name ),]$Age_group <- "14_16w"
incedence_tib[grep("17_19w", incedence_tib$name ),]$Age_group <- "17_19w"
incedence_tib[grep("20_22w", incedence_tib$name ),]$Age_group <- "20_22w"

incedence_tib$name <- gsub("\\..*", "", incedence_tib$name)
incedence_tib$name <- factor(incedence_tib$name, levels = rev(unique(incedence_tib$name)))


ggplot(data=incedence_tib,aes(x=Age_group, y=name,size = value, col=Type))  + geom_point() + facet_wrap("Type")+scale_size_continuous(range = c(5,10))+ theme_bw() +
  theme(legend.position = "none")+ theme(axis.text.y = element_text(face="bold", size=25),axis.text.x = element_text(face="bold", size=20), axis.title=element_blank(),plot.title= element_text(face="bold", size=50, hjust= 0.0))+
  theme(strip.text.x=element_text(size = 20,face = "bold"))






##### DEV scores

# I ran differtial expression between DE and OB, no need to run this again
# library(DEsingle)
# temp_cds <- odonto_trajectory_cds[,colData(odonto_trajectory_cds)$new_clusters %in% c("DE","OD")]
# group = factor(colData(temp_cds)$new_clusters)
# results <- DEsingle(counts = counts(temp_cds), group = group, parallel = T)
# results.classified <- DEtype(results = results, threshold = 0.05)
# results.sig <- results.classified[results.classified$pvalue.adj.FDR < 0.1, ]
# deg <- results.sig[,c(12,13,11,20,21,23,24)]
# deg$Gene <- rowData(temp_cds[rownames(deg),])$gene_short_name
# deg <- deg[,c(8,1,2,3,4,5,6,7)]
# deg$norm_total_mean_1 <- deg$norm_total_mean_1 *1000
# deg$norm_total_mean_2 <- deg$norm_total_mean_2 *1000
# colnames(deg) <- c("Gene","Cluster2_total_mean", "Cluster3_total_mean", "foldChange" , "pvalue", "pvalue.adj.FDR","Type", "State")
# saveRDS(deg,paste( "deg_DE_OD", ".rds",sep =''))
# 


library(Seurat)

odonto_trajectory_cds_incisor_molars_SUB <- odonto_trajectory_cds_incisor_molars[,colData(odonto_trajectory_cds_incisor_molars)$new_clusters %in% c("DE", "preOB", "OB")]
# importing monocle object  to seurat
normalized_data <- normalized_counts(odonto_trajectory_cds_incisor_molars_SUB)
rownames(normalized_data) <- as.character(rowData(odonto_trajectory_cds_incisor_molars_SUB)$gene_short_name)
clusters <- droplevels(colData(odonto_trajectory_cds_incisor_molars_SUB)$combined_trajectory)
cell_id <- colnames(odonto_trajectory_cds_incisor_molars_SUB)
meta_clusters <- data.frame(row.names=cell_id,clusters)
seurat_object <- CreateSeuratObject(counts = normalized_data, meta.data = meta_clusters)

deg_rds <- readRDS("deg_DE_OD.rds")
#deg_rds <- deg_rds[deg_rds$Type %in% c("DEa","DEg")]
deg_rds <-deg_rds %>% filter(Type %in% c("DEs","DEg"))


deg_rds_UP <- deg_rds[deg_rds$State =="down",]
deg_rds_DOWN <- deg_rds[deg_rds$State =="up",]

deg_rds_UP <-  deg_rds_UP[order(!deg_rds_UP$foldChange, deg_rds_UP$pvalue.adj.FDR),]
deg_rds_DOWN <-  deg_rds_DOWN[order(!deg_rds_DOWN$foldChange, deg_rds_DOWN$pvalue.adj.FDR),]
up_in_DE <- list((deg_rds_DOWN$Gene)[1:70])
up_in_OD<- list((deg_rds_UP$Gene)[1:20])

seurat_object <-  AddModuleScore(
  seurat_object,
  features= up_in_DE ,
  pool = NULL,
  nbin = 15,
  ctrl = 1,
  k = FALSE,
  assay = NULL,
  name = "DE_scores",
  seed = 1,
  search = FALSE
)

seurat_object <-  AddModuleScore(
  seurat_object,
  features =up_in_OD,
  pool = NULL,
  nbin = 10,
  ctrl = 1,
  k = FALSE,
  assay = NULL,
  name = "OD_scores",
  seed = 1,
  search = FALSE
)

seurat_object$Dev_Scores<- seurat_object$OD_scores1 - seurat_object$DE_scores1


# get the sort by average scorre
scores_df <- seurat_object[[]][,c("clusters","Dev_Scores")]
scores_df$clusters <- forcats::fct_explicit_na(scores_df$clusters, na_level=NA)
scores_df <- scores_df %>% group_by(clusters) %>% summarise_at(vars(Dev_Scores), funs(mean(., na.rm=F))) %>% arrange(-Dev_Scores)
scores_df$color <- ""
scores_df[grep("molar", scores_df$clusters ),]$color <- "#00BFC4"
scores_df[grep("incisor", scores_df$clusters ),]$color <- "#F8766D"
row.names(scores_df) <- scores_df$clusters

#reorder seurat object
#seurat_object$clusters <- factor(seurat_object$clusters, unique(scores_df$clusters))
sorted_levels  <- c()
incisor_n <-1
molar_n <- (length(levels(seurat_object$clusters))/2)+1

for (v in (1:length(levels(seurat_object$clusters))/2)) {
  # sorted_levels <- c(sorted_levels,
  #                    levels(seurat_object$clusters)[seq.int(v, length(levels(seurat_object$clusters)), (length(levels(seurat_object$clusters))/2)+1 )])
  sorted_levels <- c(sorted_levels, levels(seurat_object$clusters)[molar_n],levels(seurat_object$clusters)[incisor_n])
  incisor_n <- incisor_n +1
  molar_n <- molar_n +1
}
sorted_levels <- unique(sorted_levels[!is.na(sorted_levels)])
seurat_object$clusters <- factor(seurat_object$clusters, levels =  rev(sorted_levels))

library(ggtext)
RidgePlot(seurat_object, features = "Dev_Scores", group.by= "clusters", sort= F) +
  theme(legend.position = "none")+ theme(axis.text.y = element_text(face="bold", size=30),axis.text.x = element_text(face="plain", size=20), axis.title=element_blank(),plot.title= element_text(face="bold", size=50, hjust= 0.0))+
  theme(axis.text.y = element_markdown(color = scores_df[rev(sorted_levels),]$color ))

