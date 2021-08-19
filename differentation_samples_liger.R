

OE_Am_cds <- readRDS("C:/Users/afrit/Desktop/Projects/scRNA-seq/QC_Tooth_data/monocle3_cds/CDS_collecttions/OE_Am_aresidual_100d_15_0,01_40k.rds")
Amelo_cds <- readRDS("C:/Users/afrit/Desktop/Projects/scRNA-seq/QC_Tooth_data/monocle3_cds/combined_amelo_with_trajectory_cds.rds")

Amelo_cds <- readRDS("matched_amelo_trajectory_cds_new_cluster2.rds") ### new matched
Amelo_cds <- readRDS("new_amelo_trajectory_6_30_2020.rds") ### new un matched
Amelo_cds <- readRDS("combined_OE_dental_after_matching_incisors_molars.rds") ### all epi after matching
Amelo_cds <- readRDS("OE_dental_50dim_45L_0.1_10k.rds")### all epi no matching

 

plot <- plot_cells(
  Amelo_cds,
  color_cells_by = "cluster",
  group_cells_by = "cluster",
  cell_size = 1,
  group_label_size = 5,
  alpha = 1,
  label_cell_groups = F,
  genes = c("RYR2","THBD","HEY2","VAT1L","MYLK","SOX5"),
  show_trajectory_graph = F,
) #+ scale_y_reverse()
plot_data <-plot$data 
plot_data <- plot_data[order(plot_data$value, na.last = F),]
plot$data  <- plot_data
print(plot + coord_fixed(ratio = 1))+ theme(
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.title=element_text(size=14))+theme(strip.text.x=element_text(size = 14,face = "bold"),legend.title = element_text(size=14)) #scale_x_reverse()



OE_Am_cds <- preprocess_cds(OE_Am_cds, num_dim = 40, verbose = T)
OE_Am_cds <- align_cds(OE_Am_cds, alignment_group = "diff_time_point") # residual_model_formula_str ="~diff_time_point" ) 
#OE_Am_cds <- align_cds(OE_Am_cds, residual_model_formula_str ="~diff_time_point" ,alignment_group = "diff_time_point") 
#OE_Am_cds <- align_cds(OE_Am_cds, residual_model_formula_str ="~diff_time_point")
OE_Am_cds <- reduce_dimension(OE_Am_cds, umap.metric = "Cosine",   umap.n_neighbors =15L, umap.min_dist = 0.0001, verbose = T, cores =5) 
OE_Am_cds <- cluster_cells(OE_Am_cds, k = 31, verbose = T, partition_qval = 0.05, louvain_iter = 2 )

plot_cells(OE_Am_cds, color_cells_by="cluster", group_cells_by="cluster", cell_size=1, group_label_size= 7, alpha = 1, label_cell_groups = T, show_trajectory_graph = T,label_branch_points=F,label_leaves=F,trajectory_graph_segment_size=1,trajectory_graph_color= "black" ) + coord_fixed(ratio = 1)+ scale_y_reverse()+ scale_x_reverse() +theme(
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank(),
  axis.title=element_text(size=10),
  legend.title=element_text(size=16),
  line = element_line(size = 7)) 


OE_Am_cds <- learn_graph(OE_Am_cds, use_partition = T, close_loop = F, learn_graph_control=list(orthogonal_proj_tip= T, minimal_branch_len= 12,prune_graph=T))
OE_Am_cds <- order_cells(OE_Am_cds)
#plot_cells(epi_cds, color_cells_by="pseudotime", group_cells_by="cluster", cell_size=1.5, group_label_size= 5, label_cell_groups = F, trajectory_graph_color = "black" , label_groups_by_cluster = F,label_leaves = F)+ coord_fixed(ratio = 1)

plot_cells(OE_Am_cds, color_cells_by="pseudotime", group_cells_by="cluster", cell_size=1.5, group_label_size= 5, label_cell_groups = F, trajectory_graph_color = "white" , label_groups_by_cluster = F,label_leaves = F,label_roots = F,label_branch_points = F,trajectory_graph_segment_size=1)+ scale_x_reverse()+ coord_fixed(ratio = 1)+theme(
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank(),
  axis.title=element_text(size=10),
  legend.title=element_text(size=16),
  line = element_line(size = 7))

colData(OE_Am_cds)$new_partition <- dplyr::recode_factor(clusters(OE_Am_cds),
                                                         '7'='1',
                                                         '8'='1',
                                                         '1'='1',
                                                         '5'='1',
                                                         '6'='1',
                                                         '3'='1',
                                                         '2'='1',
                                                         '12'='1',
                                                         '10'='1',
                                                         '4'='1',
                                                         '9'='2',
                                                         '15'='3',
                                                         '14'='4',
                                                         '11'='5',
                                                         '13'='6')
OE_Am_cds@clusters@listData[["UMAP"]][["partitions"]] <- colData(OE_Am_cds)$new_partition
OE_Am_cds@clusters@listData[["UMAP"]][["partitions"]] <- as.factor(rep("1", ncol(OE_Am_cds)))
names(OE_Am_cds@clusters@listData[["UMAP"]][["partitions"]])<- colnames(OE_Am_cds)
partitions(OE_Am_cds)

colData(OE_Am_cds)$new_clusters <- dplyr::recode_factor(clusters(OE_Am_cds),
                                                        '7'='1',
                                                        '8'='2',
                                                        '1'='3',
                                                        '5'='4',
                                                        '6'='5',
                                                        '3'='6',
                                                        '2'='7',
                                                        '12'='8',
                                                        '10'='9',
                                                        '4'='10',
                                                        '9'='11',
                                                        '15'='12',
                                                        '14'='13',
                                                        '11'='14',
                                                        '13'='15')

saveRDS(OE_Am_cds, "OE_Am_residual_50d_15_0.0001_30k_traject_3.rds")
OE_Am_cds <- readRDS("OE_Am_residual_50d_15_0.0001_30k_traject_3.rds")



plot_density(OE_Am_cds, groups_column = "diff_time_point" , nrow = NULL,
             ncol = NULL, x_reverse = T, y_reverse =FALSE) + scale_x_reverse()

plot <- plot_cells(
  OE_Am_cds,
  color_cells_by = "cluster",
  group_cells_by = "cluster",
  cell_size = 1,
  group_label_size = 5,
  alpha = 1,
  label_cell_groups = F,
  genes = c("AMBN"),
  show_trajectory_graph = F,
) #+ scale_y_reverse()
plot_data <-plot$data 
plot_data <- plot_data[order(plot_data$value, na.last = F),]
plot$data  <- plot_data
print(plot + coord_fixed(ratio = 1))+ scale_x_reverse()+theme(
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.title=element_text(size=14))+theme(strip.text.x=element_text(size = 14,face = "bold"),legend.title = element_text(size=14))



#### Highlite the trajectory 
colData(OE_Am_cds)$differentiation_trajectory  <- NA
trajectory <-  colnames(OE_Am_cds[,colData(OE_Am_cds)$new_clusters %in% c(1:8)])
colData(OE_Am_cds[,trajectory])$differentiation_trajectory <-colData((OE_Am_cds[,colData(OE_Am_cds)$new_clusters %in% c(1:8)]))$new_clusters

colData(OE_Am_cds[,colnames(OE_Am_cds[,colData(OE_Am_cds)$new_clusters %in% c(7)])])$differentiation_trajectory <- NA

colData(OE_Am_cds)$differentiation_trajectory <- dplyr::recode(colData(OE_Am_cds)$differentiation_trajectory,
                                                               '4'='3',
                                                               '5'='4',
                                                               '6'='4',
                                                               '8'='5')

plot_cells(OE_Am_cds, color_cells_by="differentiation_trajectory", group_cells_by="cluster", cell_size=1, group_label_size= 7, alpha = 1, label_cell_groups = T, show_trajectory_graph = F,label_branch_points=F,label_leaves=F,trajectory_graph_segment_size=1,trajectory_graph_color= "black" ,labels_per_group=1, label_groups_by_cluster = F) + coord_fixed(ratio = 1)+ scale_x_reverse()+theme(
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank(),
  axis.title=element_text(size=10),
  legend.title=element_text(size=16),
  line = element_line(size = 7))

List_of_MT_genes <-
  OE_Am_cds@rowRanges@elementMetadata@listData$gene_short_name[grep("^MT-",
                                                                    OE_Am_cds@rowRanges@elementMetadata@listData$gene_short_name)]
OE_Am_cds <-
  OE_Am_cds[!OE_Am_cds@rowRanges@elementMetadata@listData[["gene_short_name"]] %in% List_of_MT_genes,]

# removeing side clusters
OE_Am_cds <- OE_Am_cds[,colData(OE_Am_cds)$differentiation_trajectory %in% c(1:5)]

list_of_cell_names <- c()
markers_cds <- OE_Am_cds[rowData(OE_Am_cds)$gene_short_name  %in% c("KRT5","AMBN","ENAM"),]
markers_cds <- detect_genes(markers_cds, min_expr = 0.000001)
markers_cds <- markers_cds[, markers_cds@colData$num_genes_expressed >0]

length_to_match <- length(colnames(Amelo_cds)) - (length(colnames(markers_cds)))

colData(OE_Am_cds)$cluster_time_point <- as.integer(monocle3::clusters(OE_Am_cds))
colData(OE_Am_cds[,colData(OE_Am_cds)$diff_time_point == "day_16"])$cluster_time_point <- colData(OE_Am_cds[,colData(OE_Am_cds)$diff_time_point == "day_16"])$cluster_time_point + 20

# a loop to downsample per cluster
for (n in unique(colData(OE_Am_cds)$cluster_time_point)) {
  r  <- which(unique(colData(OE_Am_cds)$cluster_time_point) == n)
  subset <- OE_Am_cds[, colData(OE_Am_cds)$cluster_time_point == n]
  cell_fraction <- length(colnames(subset)) / length(colnames(OE_Am_cds))
  list_of_cell_names<- c(list_of_cell_names, sample(colnames(subset),cell_fraction*length_to_match))
  
}
list_of_cell_names <- unique(c(list_of_cell_names,colnames(markers_cds)))
OE_Am_cds_down<- OE_Am_cds[,list_of_cell_names]

OE_down <- OE_Am_cds_down[,colData(OE_Am_cds_down)$diff_time_point == "day_10" ]
Am_down <- OE_Am_cds_down[,colData(OE_Am_cds_down)$diff_time_point == "day_16" ]



mat_Amelo <- counts(Amelo_cds)
rownames(mat_Amelo) <- rowData(Amelo_cds)$gene_short_name

mat_OE <- counts(OE_down)
rownames(mat_OE) <- rowData(OE_down)$gene_short_name

mat_Am <- counts(Am_down)
rownames(mat_Am) <- rowData(Am_down)$gene_short_name



library(liger)

Amelo_liger <- createLiger(list(fetal = mat_Amelo, day_10 = mat_OE, day_16 = mat_Am ), take.gene.union= T)
Amelo_liger <- liger::normalize(Amelo_liger)
Amelo_liger <- selectGenes(Amelo_liger,keep.unique =T)
Amelo_liger <- scaleNotCenter(Amelo_liger)
Amelo_liger <- optimizeALS(Amelo_liger, k = 25,  lambda = 10)
Amelo_liger <- quantile_norm(Amelo_liger)
Amelo_liger <- louvainCluster(Amelo_liger, resolution = 0.5,k = 20)
Amelo_liger <- runUMAP(Amelo_liger, distance = 'cosine', n_neighbors = 10, min_dist = 0.05)
all.plots <- plotByDatasetAndCluster(Amelo_liger, axis.labels = c('UMAP 1', 'UMAP 2'), return.plots = T,pt.size = 1, text.size = 7)
all.plots[[1]] + all.plots[[2]]+ coord_fixed(ratio = 1)+scale_x_reverse() 

 all.plots[[1]]$data

 saveRDS(Amelo_liger,"Amelo_liger_all_roorh_OE_unmatched_0.5resolution.rds")
saveRDS(Amelo_liger,"Amelo_liger_all_roorh_OE_unmatched_0.4_10_0.05_cosine.rds")
saveRDS(Amelo_liger,"Amelo_liger_all_roorh_OE_unmatched_25k_10_0.5_10_0.05_cosine.rds")
Amelo_liger <- readRDS("Amelo_liger_all_roorh_OE_unmatched.rds")
Amelo_liger <- readRDS("Amelo_liger_new_amelo_traj_unmatched.rds")
Amelo_liger <- readRDS("Amelo_liger_all_roorh_OE_unmatched_0.5resolution.rds")

plot_liger_density(Amelo_liger, nrow = NULL,
                   ncol = 3, x_reverse = T, y_reverse =F)

 runGSEA(Amelo_liger)
 p_g2 <- c()
p_g2 <- plotGene(Amelo_liger, 'AMELX', return.plots = T)
p_g2 <- append(p_g2,plotGene(Amelo_liger, 'AMELX', return.plots = T)$day_10 +scale_x_reverse())
p_g2 <- append(p_g2,plotGene(Amelo_liger, 'AMELX', return.plots = T)$day_16 +scale_x_reverse())
plot_grid(plotlist = p_g2)


p_1 <- plotGene(Amelo_liger, 'AMELX', return.plots = T)$fetal +scale_x_reverse()
p_2 <- plotGene(Amelo_liger, 'AMELX', return.plots = T)$day_10 +scale_x_reverse()
p_3 <- plotGene(Amelo_liger, 'AMELX', return.plots = T)$day_16 +scale_x_reverse()


p_g2 <- list(p_1,p_2,p_3)
p_g2 <- plotGeneViolin(Amelo_liger, 'AMBN', return.plots = T)
plot_grid(plotlist = p_g2)


p_AMELX <- plotGene(Amelo_liger, 'KRT4', return.plots = T, plot.by = "none", axis.labels=c('UMAP 1', 'UMAP 2'),pt.size = 1, zero.color="#e3e3e3") +scale_x_reverse()
plot_data <-p_AMELX$data 
plot_data <- plot_data[order(plot_data$gene , na.last = F),]
p_AMELX$data  <- plot_data
print(p_AMELX+ coord_fixed(ratio = 1))+theme(
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.title=element_text(size=14))+theme(strip.text.x=element_text(size = 14,face = "bold"),legend.title = element_text(size=20))



plotClusterProportions(Amelo_liger)
plotWordClouds(Amelo_liger)
plotFeature(Amelo_liger, "nUMI")


pdf("word_clouds.pdf")
plotWordClouds(Amelo_liger,dataset1 = "fetal",dataset2 = "day_16")
dev.off()

pdf("gene_loadings.pdf")
plotGeneLoadings(Amelo_liger)
dev.off()


calcAlignment(Amelo_liger)
calcAgreement(Amelo_liger)
# see if certain clusters are more integrated than others
calcAlignmentPerCluster(Amelo_liger)

suggestK(Amelo_liger, num.cores = 5, k.test = seq(5, 70, 5))  #20 or 30 or 50
suggestLambda(Amelo_liger, 25, num.cores=4) #5 - 10



# river
colData(Amelo_cds)$am_trajectory_named<- dplyr::recode_factor(colData(Amelo_cds)$am_trajectory,
                                                              'Am.1'='Oral_epithelium',
                                                              'Am.2'='Dental_epithelium',
                                                              'Am.3'='Inner_enamel_epithelium',
                                                              'Am.4'='Pre-Ameloblasts',
                                                              'Am.5'='Ameloblasts')
names(colData(Amelo_cds)$am_trajectory_named) <- colnames(Amelo_cds)


colData(Amelo_cds)$am_trajectory_named<- dplyr::recode_factor(colData(Amelo_cds)$new_clusters2,
                                                              '1'='OE',
                                                              '2'='DE',
                                                              '3'='EE',
                                                              '4'='PA',
                                                              '5'='Am',
                                                              '6'='EK',
                                                              '7'='SI',
                                                              '8'='SR')
names(colData(Amelo_cds)$am_trajectory_named) <- colnames(Amelo_cds)


colData(OE_Am_cds_down)$clusters_named<- dplyr::recode_factor(clusters(OE_Am_cds_down),
                                                              '4'='diff_4',
                                                              '5'='diff_5',
                                                              '1'='diff_1',
                                                              '2'='diff_2',
                                                              '8'='diff_8',
                                                              '3'='diff_3')

colData(OE_Am_cds)$clusters_named<- dplyr::recode_factor(colData(OE_Am_cds)$new_clusters,
                                                         '1' = 'diff_1',
                                                         '2' = 'diff_2',
                                                         '3' = 'diff_3',
                                                         '4' = 'diff_4',
                                                         '5' = 'diff_5',
                                                         '6' = 'diff_6',
                                                         '7' = 'diff_7',
                                                         '8' = 'diff_8',
                                                         '9' = 'diff_9',
                                                         '10' = 'diff_10',
                                                         '11' = 'diff_11',
                                                         '12' = 'diff_12',
                                                         '13' = 'diff_13',
                                                         '14' = 'diff_14',
                                                         '15' = 'diff_15')

plot_cells(OE_Am_cds_down, color_cells_by="cluster", group_cells_by="cluster", cell_size=1, group_label_size= 5, alpha = 1, label_cell_groups = T) + coord_fixed(ratio = 1)+ scale_x_reverse()


#clusters_prior <- readRDS('~/Downloads/pbmc_alignment/tenx_seqwell_clusters.RDS')

diff_clusters <- droplevels(colData(OE_Am_cds)$modified_clusters[c(rownames(Amelo_liger@H[[2]]),rownames(Amelo_liger@H[[3]]))])
diff_clusters <- droplevels(colData(OE_Am_cds)$clusters_named[c(rownames(Amelo_liger@H[[2]]),rownames(Amelo_liger@H[[3]]))])
diff_clusters <- droplevels(colData(OE_Am_cds)$merge_clusters[c(rownames(Amelo_liger@H[[2]]),rownames(Amelo_liger@H[[3]]))])
diff_clusters <- droplevels(new_liger_cluster[c(rownames(Amelo_liger@H[[2]]),rownames(Amelo_liger@H[[3]]))])

colnames(Amelo_cds) <- gsub('.{2}$', '', colnames(Amelo_cds))  # remove last 2 charcters
names(colData(Amelo_cds)$new_clusters2) <- gsub('.{2}$', '', names(colData(Amelo_cds)$new_clusters2))  # remove last 2 charcters


Amelo_trajectory <- droplevels(colData(Amelo_cds)$assigned_cell_type[rownames(Amelo_liger@H[[1]])])
Amelo_trajectory <- droplevels(colData(Amelo_cds)$am_trajectory_named[rownames(Amelo_liger@H[[1]])])
names(colData(Amelo_cds)$new_clusters)

Amelo_liger@clusters[Amelo_liger@clusters %in% c(0,1,2,3)]


# Set specific node order for cleaner visualization (can also leave to default value for 
# automatic ordering)
set_node_order = list(c(1,2,3,4,5), c(3,2,1,4,5,6,7,8), c(1, 6, 2, 3, 4, 5)) #set order
set_node_order = list(c(5,4,3,2,1), c(8,7,6,5,4,1,2,3), c(5, 3, 4, 2, 6, 1)) #set order
makeRiverplot(Amelo_liger, Amelo_trajectory, diff_clusters, min.frac = 0.1,min.cells=1, river.yscale=4, node.order =set_node_order) #node.order ="auto"
makeRiverplot(Amelo_liger, Amelo_trajectory,cluster_consensus = Amelo_liger@clusters[Amelo_liger@clusters %in% c(0,1,2,3)], diff_clusters, min.frac = 0.1,min.cells=1, river.yscale=3, node.order =set_node_order,label.cex=0.1,river.lty=0, river.usr = c(0.1, 1, 0, 1)) #node.order ="auto"
makeRiverplot(Amelo_liger, Amelo_trajectory, diff_clusters, min.frac = 0.001,min.cells=1, river.yscale=4, node.order ="auto")
river.usr = c(0, 1, -0.6, 1.6)


makeRiverplot(Amelo_liger, Amelo_trajectory, diff_clusters, min.frac = 0.1,min.cells=1, river.yscale=5)
set_node_order = list(c(1,2,3,4,5,6,7,8), c(5,1,2,8,4,3,6,7), c(7, 2,3,4,5,6,8,10,14,15,1,11))#,9,12,13)) #set order
cluster_consensus = droplevels(Amelo_liger@clusters[Amelo_liger@clusters %in% c(0,1,2,3,4,5,6,8)])
makeRiverplot(Amelo_liger, Amelo_trajectory,cluster_consensus=cluster_consensus, diff_clusters, min.frac = 0.25,min.cells=50, river.yscale=5,node.order =set_node_order)
#cluster_consensus = Amelo_liger@clusters[Amelo_liger@clusters %in% c(0,1,2,3)],

clusters_df <- as.data.frame(Amelo_liger@clusters)
clusters_df$dataset <- Amelo_liger@cell.data[["dataset"]]
summary(clusters_df)
colnames(clusters_df)[1] <- "clusters"

clusters_df <- clusters_df[clusters_df$dataset != "fetal", ]


clusters_df$clusters<- factor(clusters_df$clusters,levels = clustmd$Clusters)
clusters_df$clusters<- factor(clusters_df$clusters,levels = c(16,11,12,6,13,7,10,9,3,2,0,14,1,5,8,4,15))

levels(clusters_df$clusters) <- c(16,11,12,6,13,4,10,9,3,2,0,14,1,5,8,7,15)
  
clusters_df_sub <-clusters_df[clusters_df$clusters %in% c(4,0,1,8),]
  
ggplot(data=clusters_df, aes(x=clusters, fill=dataset)) + 
  geom_bar(show.legend = T, position = "fill")+scale_y_continuous(name="Percentage of cells",labels = scales::percent)+theme(text=element_text(size=14, face="bold")) #+coord_flip() 

ggplot(data=clusters_df_sub, aes(x=clusters, fill=dataset)) + 
  geom_bar(show.legend = T, position = "fill")+scale_y_continuous(name="Percentage of cells",labels = scales::percent)+theme(text=element_text(size=14, face="bold")) #+coord_flip() 


names(colData(Amelo_cds)$age_group)<- colnames(Amelo_cds)

clusters_df$fetal <- colData(Amelo_cds)$assigned_cell_type[rownames(clusters_df)]
clusters_df$age_group <- colData(Amelo_cds)$age_group[rownames(clusters_df)]



ggplot(data=clusters_df, aes(x=clusters, fill=fetal)) + 
  geom_bar(show.legend = T, position = "fill")+scale_y_continuous(name="Percentage of cells",labels = scales::percent)+theme(text=element_text(size=14, face="bold")) #+coord_flip() 


clusters_df <- clusters_df[clusters_df$dataset == "fetal", ]

Amelo_cds <- Amelo_cds[,colData(Amelo_cds)$new_clusters2 %in% c(1,2,3,4,5)]

Amelo_cds$am_traj <- dplyr::recode_factor(colData(Amelo_cds)$new_clusters2,
                                                 '1'='OE',
                                                 '2'='DE',
                                                 '3'='EE',
                                                 '4'='PA',
                                                 '5'='Am')

Amelo_cds$am_traj <- dplyr::recode_factor(colData(Amelo_cds)$new_clusters2,
                                          '1'='OE',
                                          '2'='DE',
                                          '3'='EE',
                                          '4'='PA',
                                          '5'='Am',
                                          '6'='EK',
                                          '7'='SI',
                                          '8'='SR')

clusters_df <-clusters_df[colnames(Amelo_cds),]
clusters_df$am_traj <- colData(Amelo_cds)$am_traj[rownames(clusters_df)]

# Density plots with semi-transparent fill
ggplot(clusters_df, aes(x=clusters, fill=am_traj )) + geom_density(alpha=.5)

# Density plots
ggplot(clusters_df, aes(x=clusters, colour=am_traj)) + geom_density()

# Overlaid histograms
ggplot(clusters_df, aes(x=clusters, fill=am_traj)) +
  geom_histogram(binwidth=.5, alpha=.5, position="identity",stat="count")

# Interleaved histograms
ggplot(clusters_df, aes(x=clusters, fill=am_traj)) +
  geom_histogram(binwidth=.5, position="dodge",stat="count")

Timepoints <- factor(clusters_df$am_traj, levels = levels(colData(Amelo_cds)$am_traj))
Timepoints <- droplevels(Timepoints)
Timescore <- as.integer(Timepoints)

clustmd <- data.frame(Timepoints=Timepoints, Clusters =clusters_df$clusters)
clustmd <- reshape2::dcast(clustmd, Clusters~Timepoints, value.var = "Clusters", fun.aggregate = length)
clustmd[, 2:ncol(clustmd)] <- t(apply(clustmd[, 2:ncol(clustmd)], 1, function(x) x/sum(x)))
clustmd$Cluster_time_score <- apply(clustmd[, 2:ncol(clustmd)], 1,
                                    function(x) sum(mapply(function(t, y) t*y, as.numeric(x), sort(unique(Timescore), decreasing = F))))

rank(clustmd$Cluster_time_score)

clustmd$rank <- rank(clustmd$Cluster_time_score)

clustmd <-clustmd[order(clustmd$rank),]

data.frame(clustmd$Clusters[-17],clustmd$Clusters[-1])




# Density plots with semi-transparent fill
ggplot(clustmd, aes(x=rating, fill=Clusters)) + geom_density(alpha=.3)


library(igraph)#  load igraph. Probably not needed unless you get an error

trajectory<- data.frame(clustmd$Clusters[-17],clustmd$Clusters[-1])

#tempora_object@layouts<-  matrix(c(new_postion_for_clusters_on_x_axis, preserve_y_axis),  ncol = 2)
#tempora_object@layouts[,2]<- tempora_object@layouts[,2]*5 # For now this is not necessary ...to scale the graph a bit taller than width, or scale down the width :Note: if you have to do that you'll have to change the sizes of the circles "vertex.size"
#tempora_object@layouts[,1]<- tempora_object@layouts[,1]/2  #  or scale down the width :Note: if you have to do that you'll have to change the sizes of the circles "vertex.size"

edge_graph <- igraph::graph_from_data_frame(d=trajectory, vertices = clustmd, directed = T)
layout <- igraph::layout_as_tree(edge_graph, flip.y =F)
layout <- layout[,c(2,1)]
y_scale <- c(min(layout[,2]-0.1),max(layout[,2]+0.1))
x_scale <- c(0,max(layout[,1]))
colours <- RColorBrewer::brewer.pal(length(levels(Timepoints)), "YlOrRd")

# plot to pdf
cairo_pdf("tempora_liger_am_traj.pdf",50,10)
par(mar=c(0,0,10,0),cex.lab=7)
igraph::plot.igraph(edge_graph, ylim=y_scale,xlim=x_scale, ylab = "", layout = layout, vertex.shape = "pie", vertex.pie = lapply(1:nrow(clustmd), function(x) as.numeric(clustmd[x,2:((length(levels(Timepoints)))+1)])),
                    vertex.pie.color=list(colours), pie.border=list(rep("white", length(levels(Timepoints)))), vertex.frame.color="white",
                    vertex.label.color="black", 
                    vertex.size=50, edge.width =5, edge.arrow.size =2
                    , vertex.label.dist=5,vertex.label.cex=7, vertex.label.degree=pi/2, rescale=F) # label is number of clusters
legend("top",inset=c(0,-0.3) , legend = levels(Timepoints), fill=colours, bty = "n", border = "black", cex = 7, horiz=TRUE)
axis(side=3, at=c(min(layout[,1]-0.1),max(layout[,1]+0.1)), labels=c("Early","Late"), las=1, cex.axis = 7)
dev.off()






## to make them un stacked
cluster_counts <- count(clusters_df,clusters)
dataset_counts <- count(clusters_df,clusters,dataset)
dataset_counts$percent <- dataset_counts$n/cluster_counts$n[dataset_counts$clusters]

ggplot(dataset_counts)+geom_bar(aes(x=clusters,y=percent,fill = dataset), stat="identity",position = position_dodge(width=0.75), width=.75)+scale_y_continuous(name="Percentage of cells",labels = scales::percent)+theme(text=element_text(size=14, face="bold"))
# position = position_dodge(width=0.5), width=.5   can be changedto postion="postion"




meta.data <-
  data.frame(
    row.names = row.names(
      colData(intercalated_acinar_trajectory_3)),
    Clusters = colData(intercalated_acinar_trajectory_3)$new_clusters_2,
    Timepoints = colData(intercalated_acinar_trajectory_3)$age_group)

ggplot(data=meta.data, aes(x=Clusters, fill=Timepoints)) + 
  geom_bar(show.legend = T,position = "fill")+scale_y_continuous(name="Percentage of cells",labels = scales::percent)+theme(text=element_text(size=14, face="bold")) #+coord_flip() 




# clustering day16alone
Am_16d_cds <- readRDS("C:/Users/afrit/Desktop/Projects/scRNA-seq/QC_Tooth_data/monocle3_cds/CDS_collecttions/Am16_diff.RDS")

Am_16d_cds <- preprocess_cds(Am_16d_cds, num_dim = 10, verbose = T)
#Am_16d_cds <- align_cds(Am_16d_cds, alignment_group = "diff_time_point") # residual_model_formula_str ="~diff_time_point" ) 
#Am_16d_cds <- align_cds(Am_16d_cds, residual_model_formula_str ="~diff_time_point" ,alignment_group = "diff_time_point") 
#Am_16d_cds <- align_cds(Am_16d_cds, residual_model_formula_str ="~diff_time_point")
Am_16d_cds <- reduce_dimension(Am_16d_cds, umap.metric = "Cosine",   umap.n_neighbors =70L, umap.min_dist = 0.0000000000001, verbose = T, cores =5) 
Am_16d_cds <- cluster_cells(Am_16d_cds, k = 20, verbose = T, partition_qval = 0.05, louvain_iter = 2 )



plot_cells(Am_16d_cds, color_cells_by="cluster", group_cells_by="cluster", cell_size=1, group_label_size= 7, alpha = 1, label_cell_groups = T, show_trajectory_graph = T,label_branch_points=F,label_leaves=F,trajectory_graph_segment_size=1,trajectory_graph_color= "black" ) + coord_fixed(ratio = 1)+ scale_x_reverse()+theme(
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank(),
  axis.title=element_text(size=10),
  legend.title=element_text(size=16),
  line = element_line(size = 7))

plot <- plot_cells(
  Am_16d_cds, 
  color_cells_by = "cluster",
  group_cells_by = "cluster",
  cell_size = 1,
  group_label_size = 5,
  alpha = 1,
  label_cell_groups = F,
  genes = c("AMBN","SHH","ENAM"), # "SOX2","PITX2","SP6",
  show_trajectory_graph = F,
) #+ scale_y_reverse()
plot_data <-plot$data 
plot_data <- plot_data[order(plot_data$value, na.last = F),]
plot$data  <- plot_data
print(plot + coord_fixed(ratio = 1))+ scale_x_reverse()+theme(
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.title=element_text(size=14))+theme(strip.text.x=element_text(size = 14,face = "bold"),legend.title = element_text(size=14))

saveRDS(Am_16d_cds, "Am_16d_cds_10dim_70L_0.0000000000001dist.rds")


OE_10d_cds <- readRDS("C:/Users/afrit/Desktop/Projects/scRNA-seq/QC_Tooth_data/monocle3_cds/CDS_collecttions/OE_10_diff.RDS")

OE_10d_cds <- preprocess_cds(OE_10d_cds, num_dim = 30, verbose = T)
#OE_10d_cds <- align_cds(OE_10d_cds, alignment_group = "diff_time_point") # residual_model_formula_str ="~diff_time_point" ) 
#OE_10d_cds <- align_cds(OE_10d_cds, residual_model_formula_str ="~diff_time_point" ,alignment_group = "diff_time_point") 
#OE_10d_cds <- align_cds(OE_10d_cds, residual_model_formula_str ="~diff_time_point")
OE_10d_cds <- reduce_dimension(OE_10d_cds, umap.metric = "Cosine",   umap.n_neighbors =20L, umap.min_dist = 0.0001, verbose = T, cores =5) 
OE_10d_cds <- cluster_cells(OE_10d_cds, k = 20, verbose = T, partition_qval = 0.05, louvain_iter = 2 )


plot_cells(OE_10d_cds, color_cells_by="cluster", group_cells_by="cluster", cell_size=1, group_label_size= 7, alpha = 1, label_cell_groups = T, show_trajectory_graph = T,label_branch_points=F,label_leaves=F,trajectory_graph_segment_size=1,trajectory_graph_color= "black" ) + coord_fixed(ratio = 1)+ scale_x_reverse()+theme(
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank(),
  axis.title=element_text(size=10),
  legend.title=element_text(size=16),
  line = element_line(size = 7))

plot <- plot_cells(
  OE_10d_cds, 
  color_cells_by = "cluster",
  group_cells_by = "cluster",
  cell_size = 1,
  group_label_size = 5,
  alpha = 1,
  label_cell_groups = F,
  genes = c( "SOX2","PITX2","SP6"),
  show_trajectory_graph = F,
) #+ scale_y_reverse()
plot_data <-plot$data 
plot_data <- plot_data[order(plot_data$value, na.last = F),]
plot$data  <- plot_data
print(plot + coord_fixed(ratio = 1))+ scale_x_reverse()+theme(
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.title=element_text(size=14))+theme(strip.text.x=element_text(size = 14,face = "bold"),legend.title = element_text(size=14))

OE_10d_cds$diff_clusters <- dplyr::recode_factor(monocle3::clusters(OE_10d_cds),
                                                    '1'='d10_c1',
                                                    '2'='d10_c2',
                                                    '3'='d10_c3',
                                                    '4'='d10_c4',
                                                    '5'='d10_c5',
                                                    '6'='d10_c6',
                                                    '7'='d10_c7',
                                                    '8'='d10_c8',
                                                    '9'='d10_c9',
                                                    '10'='d10_c10')

Am_16d_cds$diff_clusters <- dplyr::recode_factor(monocle3::clusters(Am_16d_cds),
                                                 '1'='d16_c1',
                                                 '2'='d16_c2',
                                                 '3'='d16_c3',
                                                 '4'='d16_c4',
                                                 '5'='d16_c5',
                                                 '6'='d16_c6',
                                                 '7'='d16_c7',
                                                 '8'='d16_c8',
                                                 '9'='d16_c9',
                                                 '10'='d16_c10')

OE_Am_cds  <- combine_cds(list(OE_10d_cds,Am_16d_cds))
names(colData(OE_Am_cds)$diff_clusters) <- colnames(OE_Am_cds)


OE_Am_cds <- preprocess_cds(OE_Am_cds, num_dim = 10, verbose = T)
#OE_Am_cds <- align_cds(OE_Am_cds, alignment_group = "diff_time_point") # residual_model_formula_str ="~diff_time_point" ) 
#OE_Am_cds <- align_cds(OE_Am_cds, residual_model_formula_str ="~diff_time_point" ,alignment_group = "diff_time_point") 
#OE_Am_cds <- align_cds(OE_Am_cds, residual_model_formula_str ="~diff_time_point")
OE_Am_cds <- reduce_dimension(OE_Am_cds, umap.metric = "hamming",   umap.n_neighbors =70L, umap.min_dist = 0.0000000000001, verbose = T, cores =5) 
OE_Am_cds <- cluster_cells(OE_Am_cds, k = 50, verbose = T, partition_qval = 0.05, louvain_iter = 2 )



plot_cells(OE_Am_cds, color_cells_by="cluster", group_cells_by="cluster", cell_size=1, group_label_size= 7, alpha = 1, label_cell_groups = T, show_trajectory_graph = T,label_branch_points=F,label_leaves=F,trajectory_graph_segment_size=1,trajectory_graph_color= "black" ) + coord_fixed(ratio = 1)+ scale_x_reverse()+theme(
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank(),
  axis.title=element_text(size=10),
  legend.title=element_text(size=16),
  line = element_line(size = 7))

plot <- plot_cells(
  OE_Am_cds, 
  color_cells_by = "cluster",
  group_cells_by = "cluster",
  cell_size = 1,
  group_label_size = 5,
  alpha = 1,
  label_cell_groups = F,
  genes = c("AMBN","SHH","ENAM"), # "SOX2","PITX2","SP6",
  show_trajectory_graph = F,
) #+ scale_y_reverse()
plot_data <-plot$data 
plot_data <- plot_data[order(plot_data$value, na.last = F),]
plot$data  <- plot_data
print(plot + coord_fixed(ratio = 1))+ scale_x_reverse()+theme(
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.title=element_text(size=14))+theme(strip.text.x=element_text(size = 14,face = "bold"),legend.title = element_text(size=14))


OE_Am_cds$merge_clusters <- dplyr::recode_factor(monocle3::clusters(OE_Am_cds),
                                                 '1'='diff_c1',
                                                 '2'='diff_c2',
                                                 '3'='diff_c3',
                                                 '4'='diff_c4',
                                                 '5'='diff_c5',
                                                 '6'='diff_c6',
                                                 '7'='diff_c7',
                                                 '8'='diff_c8',
                                                 '9'='diff_c9',
                                                 '10'='diff_c10',
                                                '11'='diff_c11',
                                                '12'='diff_c12',
                                                '13'='diff_c13',
                                                '14'='diff_c14',
                                                '15'='diff_c15',
                                                '16'='diff_c16',
                                                '17'='diff_c17',
                                                '18'='diff_c18',
                                                '19'='diff_c19')

# extracting clusters per sets
clusters_df <- as.data.frame(Amelo_liger@clusters)
clusters_df$dataset <- Amelo_liger@cell.data[["dataset"]]
summary(clusters_df)
colnames(clusters_df)[1] <- "clusters"

clusters_df <- clusters_df[clusters_df$dataset != "fetal", ]

amelo_c <-rownames(clusters_df[clusters_df$clusters =="8",])
oe_c <-rownames(clusters_df[clusters_df$clusters =="4",])
de_c <-rownames(clusters_df[clusters_df$clusters =="1",])
ee_c <-rownames(clusters_df[clusters_df$clusters =="1",])
pa_c <-rownames(clusters_df[clusters_df$clusters =="3",])
ek_c <-rownames(clusters_df[clusters_df$clusters =="4",])
si_c <-rownames(clusters_df[clusters_df$clusters =="0",])



OE_Am_cds$matched<- NA
OE_Am_cds[,oe_c]$matched<- "oe_c"
OE_Am_cds[,amelo_c]$matched<- "amelo_c"
OE_Am_cds[,de_c]$matched<- "de_c"
OE_Am_cds[,pa_c]$matched<- "pa_c"
OE_Am_cds[,ek_c]$matched<- "ek_c"
OE_Am_cds[,ee_c]$matched<- "ee_c"
OE_Am_cds[,si_c]$matched<- "si_c"


sub_ee_c <- sample(ee_c,length(ee_c)*0.2)
sub_oe_c <- sample(oe_c,length(oe_c)*0.2)


# updated "modified"
colData(OE_Am_cds)$modified_clusters<- colData(OE_Am_cds)$clusters_named
OE_Am_cds[,amelo_c]$modified_clusters<- "diff_8"
OE_Am_cds[,ee_c]$modified_clusters<- "diff_6"
OE_Am_cds[,sub_ee_c]$modified_clusters<- "diff_8"

colData(OE_Am_cds)$modified_clusters<-  dplyr::recode(colData(OE_Am_cds)$modified_clusters,
                                                             'diff_1'='diff_7',
                                                             'diff_7'='diff_1')

OE_Am_cds[,sub_oe_c]$modified_clusters<- "diff_2"


names(OE_Am_cds$matched) <- colnames(OE_Am_cds)

plot <- plot_cells(
  OE_Am_cds, 
  color_cells_by = "matched",
  group_cells_by = "cluster",
  cell_size = 1,
  group_label_size = 5,
  alpha = 1,
  label_cell_groups = F,
  show_trajectory_graph = F,
) #+ scale_y_reverse()
plot_data <-plot$data 
plot_data <- plot_data[order(plot_data$cell_color, na.last = F),]
plot$data  <- plot_data
print(plot + coord_fixed(ratio = 1))+ scale_x_reverse()+theme(
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.title=element_text(size=14))+theme(strip.text.x=element_text(size = 14,face = "bold"),legend.title = element_text(size=14))




### liger for only diff samples
OE_down <- OE_Am_cds[,colData(OE_Am_cds)$diff_time_point == "day_10" ]
Am_down <- OE_Am_cds[,colData(OE_Am_cds)$diff_time_point == "day_16" ]


mat_OE <- counts(OE_down)
rownames(mat_OE) <- rowData(OE_down)$gene_short_name

mat_Am <- counts(Am_down)
rownames(mat_Am) <- rowData(Am_down)$gene_short_name



library(liger)

Diff_liger <- createLiger(list(day_10 = mat_OE, day_16 = mat_Am ), take.gene.union= T)
Diff_liger <- liger::normalize(Diff_liger)
Diff_liger <- selectGenes(Diff_liger,keep.unique =T)
Diff_liger <- scaleNotCenter(Diff_liger)
Diff_liger <- optimizeALS(Diff_liger, k = 10,  lambda = 5)
Diff_liger <- quantile_norm(Diff_liger)
Diff_liger <- louvainCluster(Diff_liger, resolution = 0.4,  k = 20)
#Diff_liger <- runTSNE(Diff_liger, use.pca = T, perplexity =50,  method = "fftRtsne",fitsne.path = "C:",)
Diff_liger <- runUMAP(Diff_liger, distance = 'hamming', n_neighbors = 10, min_dist = 0.000001)
all.plots <- plotByDatasetAndCluster(Diff_liger, axis.labels = c('UMAP 1', 'UMAP 2'), return.plots = T,pt.size = 1, clusters = OE_Am_cds$matched)
all.plots[[1]] + all.plots[[2]]+ coord_fixed(ratio = 1)


plot_liger_density(Diff_liger, nrow = NULL,
                   ncol = 1, x_reverse = F, y_reverse =F)


p_g2 <- plotGene(Diff_liger, 'AMBN', return.plots = T)
plot_grid(plotlist = p_g2)

suggestK(Diff_liger, num.cores = 5, k.test = seq(5, 70, 5))  #20 or 30 or 50
suggestLambda(Diff_liger, 20, num.cores=6) #5 - 10


new_liger_cluster <- dplyr::recode_factor(Diff_liger@clusters,
                                          '0'='diff_c0',
                                          '1'='diff_c1',
                                          '2'='diff_c2',
                                          '3'='diff_c3',
                                          '4'='diff_c4',
                                          '5'='diff_c5',
                                          '6'='diff_c6',
                                          '7'='diff_c7',
                                          '8'='diff_c8',
                                          '9'='diff_c9',
                                          '10'='diff_c10',
                                          '11'='diff_c11',
                                          '12'='diff_c12',
                                          '13'='diff_c13',
                                          '14'='diff_c14',
                                          '15'='diff_c15',
                                          '16'='diff_c16',
                                          '17'='diff_c17',
                                          '18'='diff_c18',
                                          '19'='diff_c19') 

OE_Am_cds[, colData(OE_Am_cds)$diff_time_point == "day_16" & colData(OE_Am_cds)$new_clusters ==8 ]
1608/35092




# heat map forfetal and difff

Amelo_traj <- readRDS("new_amelo_trajectory_6_30_2020.rds") ### new un matched
sub_Amelo_traj <- readRDS("ameloblast_only.rds") ### new un matched

OE_Am_traj <- OE_Am_cds[,clusters(OE_Am_cds) %in% c(1,2,3,4,5,6,7,8)]

selected_markers <- c("ENAM","AMELX","AMBN","WDR72","AMTN","ITGB6","ACP4","CLDN16","FAM83H","FAM20C","MMP20","KLK4","SP6","RELT","SLC24A4","DLX3")

cell_group_df <- data.frame(cell=colnames(OE_Am_traj) , group = clusters(OE_Am_traj)) #clusters(OE_Am_traj)
mat_agg <- aggregate_gene_expression(OE_Am_traj,cell_group_df = cell_group_df, norm_method ="size_only" )
row.names(mat_agg) = rowData(OE_Am_traj[row.names(mat_agg),])$gene_short_name

setdiff(selected_markers,row.names(mat_agg))

# This to subset and order the matrix by the selected genes. Note it might throw  an error if you have a gene not in the dataset
mat_agg_selected <- mat_agg[selected_markers,]


#plotting
pheatmap::pheatmap(mat_agg_selected,
                   scale="row",  cluster_cols = F,cluster_rows=F, show_rownames= T, border_color = "grey") #gaps_row = cumsum(gap_indecies), gaps_col = cumsum(c(1,1,1,3,1,2,1,1)))




sub_Amelo_traj@clusters@listData[["UMAP"]][["clusters"]]<- dplyr::recode_factor(clusters(sub_Amelo_traj),
                                          '2'='1',
                                          '1'='2')


 
# number of cells in the right trajectroy
table(colData(OE_Am_cds)$clusters_named) 
38955


###########################33333333
##########################33
###################################3
########################3

OE_Am_cds <- readRDS("OE_Am_cds_may_27.rds")

# creating simplified label
colData(Amelo_cds_MTi)$simplified <- dplyr::recode_factor(colData(Amelo_cds_MTi)$assigned_cell_type.simp,
                                                          "OE" = "OE",
                                                          "DE-Prog" = "DE-Prog",
                                                          "IK/EK" = "IK/EK",
                                                          "OEE" = "OEE", 
                                                          "SR-1" = "SR",
                                                          "SR-2" = "SR", 
                                                          "SI-1" = "SI", 
                                                          "SI-2" = "SI", 
                                                          "SI-3" = "SI", 
                                                          "PA-1" = "PA", 
                                                          "PA-2" = "PA",
                                                          "AM-1" = "AM", 
                                                          "AM-2" = "AM")



diff_clusters <- droplevels(colData(OE_Am_cds)$modified_clusters[c(rownames(Amelo_liger@H[[2]]),rownames(Amelo_liger@H[[3]]))])
diff_clusters <- droplevels(colData(OE_Am_cds)$clusters_named[c(rownames(Amelo_liger@H[[2]]),rownames(Amelo_liger@H[[3]]))])
diff_clusters <- droplevels(colData(OE_Am_cds)$merge_clusters[c(rownames(Amelo_liger@H[[2]]),rownames(Amelo_liger@H[[3]]))])
diff_clusters <- droplevels(new_liger_cluster[c(rownames(Amelo_liger@H[[2]]),rownames(Amelo_liger@H[[3]]))])

colnames(Amelo_cds) <- gsub('.{2}$', '', colnames(Amelo_cds))  # remove last 2 charcters
names(colData(Amelo_cds)$new_clusters2) <- gsub('.{2}$', '', names(colData(Amelo_cds)$new_clusters2))  # remove last 2 charcters


Amelo_trajectory <- droplevels(colData(Amelo_cds_MTi)$assigned_cell_type.simp[rownames(Amelo_liger@H[[1]])])
#Amelo_trajectory <- droplevels(colData(Amelo_cds)$am_trajectory_named[rownames(Amelo_liger@H[[1]])])
names(colData(Amelo_cds)$new_clusters)

Amelo_liger@clusters[Amelo_liger@clusters %in% c(0,1,2,3)]


# Set specific node order for cleaner visualization (can also leave to default value for 
# automatic ordering)
set_node_order = list(c(1,2,3,4,5), c(3,2,1,4,5,6,7,8), c(1, 6, 2, 3, 4, 5)) #set order
set_node_order = list(c(5,4,3,2,1), c(8,7,6,5,4,1,2,3), c(5, 3, 4, 2, 6, 1)) #set order
makeRiverplot(Amelo_liger, Amelo_trajectory, diff_clusters, min.frac = 0.1,min.cells=1, river.yscale=4, node.order =set_node_order) #node.order ="auto"
makeRiverplot(Amelo_liger, Amelo_trajectory,cluster_consensus = Amelo_liger@clusters[Amelo_liger@clusters %in% c(0,1,2,3)], diff_clusters, min.frac = 0.1,min.cells=1, river.yscale=3, node.order =set_node_order,label.cex=0.1,river.lty=0, river.usr = c(0.1, 1, 0, 1)) #node.order ="auto"
makeRiverplot(Amelo_liger, Amelo_trajectory, diff_clusters, min.frac = 0.10,min.cells=1, river.yscale=4, node.order ="auto")
river.usr = c(0, 1, -0.6, 1.6)


makeRiverplot(Amelo_liger, Amelo_trajectory, diff_clusters, min.frac = 0.1,min.cells=1, river.yscale=5)
set_node_order = list(c(1,2,3,4,5,6,7,8), c(5,1,2,8,4,3,6,7), c(7, 2,3,4,5,6,8,10,14,15,1,11))#,9,12,13)) #set order
cluster_consensus = droplevels(Amelo_liger@clusters[Amelo_liger@clusters %in% c(0,1,2,3,4,5,6,8)])
makeRiverplot(Amelo_liger, Amelo_trajectory,cluster_consensus=cluster_consensus, diff_clusters, min.frac = 0.25,min.cells=50, river.yscale=5,node.order =set_node_order)
#cluster_consensus = Amelo_liger@clusters[Amelo_liger@clusters %in% c(0,1,2,3)],






#####NEW CLUSTERING ######

OE_Am_cds<- readRDS("OE_Am_cds_may22.rds")

List_of_MT_genes <-
  OE_Am_cds@rowRanges@elementMetadata@listData$gene_short_name[grep("^MT-",
                                                                    OE_Am_cds@rowRanges@elementMetadata@listData$gene_short_name)]
OE_Am_cds <-
  OE_Am_cds[!OE_Am_cds@rowRanges@elementMetadata@listData[["gene_short_name"]] %in% List_of_MT_genes,]

# removeing side clusters
OE_Am_cds <- OE_Am_cds[,colData(OE_Am_cds)$differentiation_trajectory %in% c(1:5)]

list_of_cell_names <- c()
markers_cds <- OE_Am_cds[rowData(OE_Am_cds)$gene_short_name  %in% c("KRT5","AMBN","ENAM","SP6","PITX2"),]
markers_cds <- detect_genes(markers_cds, min_expr = 0.000001)
markers_cds <- markers_cds[, markers_cds@colData$num_genes_expressed >0]

length_to_match <- length(colnames(Amelo_cds_MTi)) - (length(colnames(markers_cds)))

colData(OE_Am_cds)$cluster_time_point <- as.integer(colData(OE_Am_cds)$new_clusters)  #monocle3::clusters(OE_Am_cds)
colData(OE_Am_cds[,colData(OE_Am_cds)$diff_time_point == "day_16"])$cluster_time_point <- colData(OE_Am_cds[,colData(OE_Am_cds)$diff_time_point == "day_16"])$cluster_time_point + 20

# a loop to downsample per cluster
for (n in unique(colData(OE_Am_cds)$cluster_time_point)) {
  r  <- which(unique(colData(OE_Am_cds)$cluster_time_point) == n)
  subset <- OE_Am_cds[, colData(OE_Am_cds)$cluster_time_point == n]
  cell_fraction <- length(colnames(subset)) / length(colnames(OE_Am_cds))
  list_of_cell_names<- c(list_of_cell_names, sample(colnames(subset),cell_fraction*length_to_match))
  
}
list_of_cell_names <- unique(c(list_of_cell_names,colnames(markers_cds)))
OE_Am_cds_down<- OE_Am_cds[,list_of_cell_names]

OE_down <- OE_Am_cds_down[,colData(OE_Am_cds_down)$diff_time_point == "day_10" ]
Am_down <- OE_Am_cds_down[,colData(OE_Am_cds_down)$diff_time_point == "day_16" ]



mat_Amelo <- counts(Amelo_cds_MTi)
rownames(mat_Amelo) <- rowData(Amelo_cds_MTi)$gene_short_name

mat_OE <- counts(OE_down)
rownames(mat_OE) <- rowData(OE_down)$gene_short_name

mat_Am <- counts(Am_down)
rownames(mat_Am) <- rowData(Am_down)$gene_short_name



library(rliger)

Amelo_liger <- createLiger(list(fetal = mat_Amelo, day_10 = mat_OE, day_16 = mat_Am ), take.gene.union= T)
Amelo_liger <- rliger::normalize(Amelo_liger)
Amelo_liger <- selectGenes(Amelo_liger,var.thresh = 0.0001)  # combine = "union", "intersection"
Amelo_liger <- scaleNotCenter(Amelo_liger)
Amelo_liger <- optimizeALS(Amelo_liger, k = 10,  lambda = 10)
Amelo_liger <- quantile_norm(Amelo_liger)
Amelo_liger <- louvainCluster(Amelo_liger, resolution = 0.3,k = 20)
Amelo_liger <- runUMAP(Amelo_liger, distance = 'cosine', n_neighbors = 10, min_dist = 0.5)
all.plots <- plotByDatasetAndCluster(Amelo_liger, axis.labels = c('UMAP 1', 'UMAP 2'), return.plots = T,pt.size = 1, text.size = 7)
all.plots[[1]] + all.plots[[2]]+ coord_fixed(ratio = 1)+scale_x_reverse() 

all.plots[[1]]$data

saveRDS(Amelo_liger,"Amelo_liger_NEW_FETAL_CLUSTER_0.3resolution_k20_n10_dist0.5_.rds")


plotGene(Amelo_liger, 'AMBN', return.plots = T)

p_1 <- all.plots[[2]]+ coord_fixed(ratio = 1)#+scale_x_reverse()
p_2 <- plotGene(Amelo_liger, 'AMBN', return.plots = T)$fetal + coord_fixed(ratio = 1)#+scale_x_reverse()
p_3 <- plotGene(Amelo_liger, 'ABMN', return.plots = T)$day_10 +coord_fixed(ratio = 1)#+scale_x_reverse()
p_4 <- plotGene(Amelo_liger, 'AMBN', return.plots = T)$day_16 +coord_fixed(ratio = 1)#+scale_x_reverse()



p_g2 <- list(p_1,p_2,p_3,p_4)

plot_grid(plotlist = p_g2)


diff_clusters <- droplevels(colData(OE_Am_cds_down)$clusters_named[c(rownames(Amelo_liger@H[[2]]),rownames(Amelo_liger@H[[3]]))])

diff_clusters <- droplevels(OE_Am_cds_down_updated$simplE_new_cluster[c(rownames(Amelo_liger@H[[2]]),rownames(Amelo_liger@H[[3]]))])
diff_clusters <- droplevels(clusters(OE_Am_cds_down_updated)[c(rownames(Amelo_liger@H[[2]]),rownames(Amelo_liger@H[[3]]))])

diff_clusters <- droplevels(clusters(OE_Am_cds_updated)[c(rownames(Amelo_liger@H[[2]]),rownames(Amelo_liger@H[[3]]))])

diff_clusters <- droplevels(diff_clusters[diff_clusters %in% c("diff_1","diff_2","diff_3","diff_4","diff_5","diff_6","diff_7","diff_8","diff_9","diff_10")])

Amelo_trajectory <- droplevels(colData(Amelo_cds_MTi)$assigned_cell_type.simp[rownames(Amelo_liger@H[[1]])])
Amelo_trajectory <- droplevels(colData(Amelo_cds_MTi)$assigned_cell_type.4[rownames(Amelo_liger@H[[1]])])
#Amelo_trajectory <- droplevels(colData(Amelo_cds)$am_trajectory_named[rownames(Amelo_liger@H[[1]])])

set_node_order = list(c(5,6,7,8,9,3,1,2,4,10,11,12), c(1,2,3,6,5,4,7), c(1:10)) #set o

set_node_order = list(c(1,3,5,6,8,9,4,2,7,10,11,12,13), c(6,3,1,2,5,4,7), c(1,2,4,5,6,7,9,3,10,8)) #set o

set_node_order = list(c(1,3,5,6,8,9,4,2,7,10,11,12,13), c(6,3,1,2,5,4,7), c(1,2,4,5,6)) #set o
set_node_order = list(c(1,3,5,6,8,9,4,2,7,10,11,12,13), c(6,3,1,2,5,4,7), c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)) #set o



makeRiverplot(Amelo_liger, Amelo_trajectory, diff_clusters,cluster_consensus = Amelo_liger@clusters[Amelo_liger@clusters %in% c(0:6)], min.frac = 0.11,min.cells=5, river.yscale=4, node.order =set_node_order)
river.usr = c(0, 1, -0.6, 1.6)

makeRiverplot(Amelo_liger, Amelo_trajectory, diff_clusters,cluster_consensus = Amelo_liger@clusters[Amelo_liger@clusters %in% c(0:6)], min.frac = 0.01,min.cells=5, river.yscale=5, node.order =set_node_order)


Amelo_liger@H



# extracting clusters per sets
clusters_df <- 
clusters_df$dataset <- Amelo_liger@cell.data[["dataset"]]
summary(clusters_df)
colnames(clusters_df)[1] <- "clusters"

clusters_df <- clusters_df[clusters_df$dataset != "fetal", ]

amelo_c <-rownames(clusters_df[clusters_df$clusters =="8",])
oe_c <-rownames(clusters_df[clusters_df$clusters =="4",])
de_c <-rownames(clusters_df[clusters_df$clusters =="1",])
ee_c <-rownames(clusters_df[clusters_df$clusters =="1",])
pa_c <-rownames(clusters_df[clusters_df$clusters =="3",])
ek_c <-rownames(clusters_df[clusters_df$clusters =="4",])
si_c <-rownames(clusters_df[clusters_df$clusters =="0",])

clusters_df$pre_cluster <- c(colData(Amelo_cds_MTi)$assigned_cell_type.4,colData(OE_Am_cds_down)$clusters_named)

clusters_df$projected_cluster <- dplyr::recode_factor(clusters_df$clusters,
                                                      "0" = "SR_OEE",
                                                      "1" = "SI2_3_DE",
                                                      "2" = "EK",
                                                      "3" = "PA_OEE",
                                                      "4" = "DE_SI_",
                                                      "5" = "OE",
                                                      "6" = "AM")

names(clusters_df$projected_cluster) <- row.names(clusters_df)

colData(OE_Am_cds_down)$projected_cluster <- clusters_df$projected_cluster[(row.names(clusters_df) %in% colnames(OE_Am_cds_down))]

names(colData(OE_Am_cds_down)$projected_cluster) <- colnames(OE_Am_cds_down)


plot_cells(OE_Am_cds_down, color_cells_by="projected_cluster", group_cells_by="cluster", cell_size=1, group_label_size= 7, alpha = 1, label_cell_groups = T, show_trajectory_graph = F,label_branch_points=F,label_leaves=F,trajectory_graph_segment_size=1,trajectory_graph_color= "black" ) + coord_fixed(ratio = 1) + scale_x_reverse() +theme(
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank(),
  axis.title=element_text(size=10),
  legend.title=element_text(size=16),
  line = element_line(size = 7)) #+ scale_y_reverse()+ scale_x_reverse()


OE_Am_cds_down_updated <- OE_Am_cds_down
OE_Am_cds_down_updated <- preprocess_cds(OE_Am_cds_down_updated, num_dim = 100, verbose = T)
OE_Am_cds_down_updated <- align_cds(OE_Am_cds_down_updated, alignment_group = "diff_time_point") # residual_model_formula_str ="~diff_time_point" ) 
#OE_Am_cds_down_updated <- align_cds(OE_Am_cds_down_updated, residual_model_formula_str ="~diff_time_point" ,alignment_group = "diff_time_point") 
#OE_Am_cds_down_updated <- align_cds(OE_Am_cds_down_updated, residual_model_formula_str ="~diff_time_point")
OE_Am_cds_down_updated <- reduce_dimension(OE_Am_cds_down_updated, umap.metric = "euclidean",   umap.n_neighbors =10L, umap.min_dist = 0.2, verbose = T, cores =5, preprocess_method = "PCA") 
OE_Am_cds_down_updated <- cluster_cells(OE_Am_cds_down_updated, k = 5, verbose = T, partition_qval = 0.05, louvain_iter = 2 )

plot_cells(OE_Am_cds_down_updated, color_cells_by="cluster", group_cells_by="cluster", cell_size=3, group_label_size= 7, alpha = 1, label_cell_groups = T, show_trajectory_graph = F,label_branch_points=F,label_leaves=F,trajectory_graph_segment_size=1,trajectory_graph_color= "black" ) + coord_fixed(ratio = 1) + scale_x_reverse() +theme(
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank(),
  axis.title=element_text(size=10),
  legend.title=element_text(size=16),
  line = element_line(size = 7)) #+ scale_y_reverse()+ scale_x_reverse()


plot_cells(OE_Am_cds_down_updated, color_cells_by="projected_cluster", group_cells_by="cluster", cell_size=3, group_label_size= 7, alpha = 1, label_cell_groups = F, show_trajectory_graph = F,label_branch_points=F,label_leaves=F,trajectory_graph_segment_size=1,trajectory_graph_color= "black" ) + coord_fixed(ratio = 1) + scale_x_reverse() +theme(
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank(),
  axis.title=element_text(size=10),
  legend.title=element_text(size=16),
  line = element_line(size = 7)) #+ scale_y_reverse()+ scale_x_reverse()


plot <- plot_cells(
  OE_Am_cds_down_updated,
  color_cells_by = "cluster",
  group_cells_by = "cluster",
  cell_size = 1,
  group_label_size = 5,
  alpha = 1,
  label_cell_groups = F,
  genes = c("SOX2","PITX2","SP6","AMBN"),
  show_trajectory_graph = F,
) #+ scale_y_reverse()
plot_data <-plot$data 
plot_data <- plot_data[order(plot_data$value, na.last = F),]
plot$data  <- plot_data
print(plot + coord_fixed(ratio = 1))+ theme(
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.title=element_text(size=14))+theme(strip.text.x=element_text(size = 14,face = "bold"),legend.title = element_text(size=14)) + scale_x_reverse()


OE_Am_cds_down_updated$simplE_new_cluster <- dplyr::recode_factor(clusters(OE_Am_cds_down_updated),
                                                      "1" = "PA",
                                                      "2" = "EPI",
                                                      "3" = "OE",
                                                      "4" = "OEE",
                                                      "5" = "OEE",
                                                      "6" = "DE",
                                                      "7" = "STEM")



#### full diff reclustering

OE_Am_cds_updated <- OE_Am_cds
OE_Am_cds_updated <- preprocess_cds(OE_Am_cds_updated, num_dim = 150, verbose = T)
OE_Am_cds_updated <- align_cds(OE_Am_cds_updated, alignment_group = "diff_time_point") # residual_model_formula_str ="~diff_time_point" ) 
#OE_Am_cds_updated <- align_cds(OE_Am_cds_updated, residual_model_formula_str ="~diff_time_point" ,alignment_group = "diff_time_point") 
#OE_Am_cds_updated <- align_cds(OE_Am_cds_updated, residual_model_formula_str ="~diff_time_point")
OE_Am_cds_updated <- reduce_dimension(OE_Am_cds_updated, umap.metric = "manhattan",   umap.n_neighbors =10L, umap.min_dist = 0.25, verbose = T, cores =5, preprocess_method = "PCA") 
OE_Am_cds_updated <- cluster_cells(OE_Am_cds_updated, k = 30, verbose = T, partition_qval = 0.05, louvain_iter = 2 )

plot_cells(OE_Am_cds_updated, color_cells_by="cluster", group_cells_by="cluster", cell_size=3, group_label_size= 7, alpha = 1, label_cell_groups = T, show_trajectory_graph = F,label_branch_points=F,label_leaves=F,trajectory_graph_segment_size=1,trajectory_graph_color= "black" ) + coord_fixed(ratio = 1) + scale_x_reverse() +theme(
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank(),
  axis.title=element_text(size=10),
  legend.title=element_text(size=16),
  line = element_line(size = 7)) #+ scale_y_reverse()+ scale_x_reverse()


plot_cells(OE_Am_cds_updated, color_cells_by="projected_cluster", group_cells_by="cluster", cell_size=3, group_label_size= 7, alpha = 1, label_cell_groups = F, show_trajectory_graph = F,label_branch_points=F,label_leaves=F,trajectory_graph_segment_size=1,trajectory_graph_color= "black" ) + coord_fixed(ratio = 1) + scale_x_reverse() +theme(
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank(),
  axis.title=element_text(size=10),
  legend.title=element_text(size=16),
  line = element_line(size = 7)) #+ scale_y_reverse()+ scale_x_reverse()


plot <- plot_cells(
  OE_Am_cds_updated,
  color_cells_by = "cluster",
  group_cells_by = "cluster",
  cell_size = 1,
  group_label_size = 5,
  alpha = 1,
  label_cell_groups = F,
  genes = c("SOX2","PITX2","SP6","AMBN"),
  show_trajectory_graph = F,
) #+ scale_y_reverse()
plot_data <-plot$data 
plot_data <- plot_data[order(plot_data$value, na.last = F),]
plot$data  <- plot_data
print(plot + coord_fixed(ratio = 1))+ theme(
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.title=element_text(size=14))+theme(strip.text.x=element_text(size = 14,face = "bold"),legend.title = element_text(size=14)) + scale_x_reverse()


OE_Am_cds_updated$simplE_new_cluster <- dplyr::recode_factor(clusters(OE_Am_cds_updated),
                                                                  "1" = "PA",
                                                                  "2" = "EPI",
                                                                  "3" = "OE",
                                                                  "4" = "OEE",
                                                                  "5" = "OEE",
                                                                  "6" = "DE",
                                                                  "7" = "STEM")



plot_density(OE_Am_cds_updated, groups_column = "diff_time_point" , nrow = NULL,
             ncol = NULL, x_reverse = T, y_reverse =FALSE) + scale_x_reverse()

saveRDS(OE_Am_cds_updated, "OE_Am_cds_updated_manhattan_10L_0.25_25K.rds")
saveRDS(OE_Am_cds_updated, "OE_Am_cds_updated_manhattan_10L_0.25_30K.rds")
