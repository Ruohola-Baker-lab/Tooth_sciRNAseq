


#amelo_trajectory <- readRDS("amelo_trajectory.rds")

#plot_heatmap_over_realtime(amelo_trajectory,"assigned_cell_type.4",age_clustering_column="age_group", selected_genes= c("SP6","DSPP","AMBN"), normalize = T)

plot_heatmap_over_realtime <- function(cds, clustering_column, age_clustering_column, selected_genes,  normalize = T, trimm_low= T){

start_time <- Sys.time()
library(tidyr)
library(monocle3)
library(ggplot2)
library(dplyr)
library(Matrix.utils)
library(ComplexHeatmap)
library(circlize)

assertthat::assert_that(!is.null(cds), msg = "cds object doesn't exist. Check the spelling!")
assertthat::validate_that(!is.null(clustering_column),msg = "please add clustering_column name. Otherwise, the default clustering will be loaded")
if (clustering_column == "cluster"){
  colData(cds)$clustering <- clusters(cds)
  clustering_column<- "clustering"
  
}
assertthat::assert_that(!is.null(age_clustering_column),msg = "please add age_clustering_column name")
assertthat::validate_that(!is.null(colData(cds)[[clustering_column]]), msg = "clustering_column doesn't exist in cds. Check the spelling!")
assertthat::assert_that(!is.null(colData(cds)[[age_clustering_column]]), msg = "age_clustering_column doesn't exist in cds. Check the spelling!")
assertthat::assert_that(!is.null(rowData(cds)$gene_short_name %in% selected_genes), msg = "selected_genes doesn't exist in cds. Check the spelling!")

cds_name <-   deparse(match.call()$cds) # deparse get the charcter nameof an object, while gsub here remove . in the name #gsub("\\.", "_", .deparse(substitute(cds)) doesn't work in function
data.dir <- cds_name  # Create a directory to save data
dir.create("realtime_Heatmap_output", showWarnings = F)
#data_file <- paste("realtime_Heatmap_output","/",cds_name, "_realheatmapData.Rdata", sep="")


print("Running precalculations")
age_group_cluster_tab <- table(colData(cds)[,c(age_clustering_column,clustering_column)])
#sample_id_cluster_sum <- table(colData(cds)[,c("sample_id",clustering_column)])
age_group_sum<- rowSums(age_group_cluster_tab)
cluster_sum <- colSums(age_group_cluster_tab)
age_group_size_factor <- (age_group_cluster_tab/age_group_sum)
n_clusters <- ncol(age_group_cluster_tab)
n_age  <- nrow(age_group_cluster_tab)
n_genes  <- length(selected_genes)
#print(age_group_cluster_tab)
#normalized_cell_count <- t(t(age_group_size_factor)/cluster_sum)
#normalized_cell_count <- t(t(cluster_size_factor)/age_group_sum)

# cluster_size_factor <- (age_group_cluster_tab)/cluster_sum
# normalized_cell_count <-  age_group_cluster_tab *(age_group_size_factor*cluster_size_factor)

 cluster_size_factor <- t(t(age_group_cluster_tab)/cluster_sum)
# normalized_cell_count <-  age_group_cluster_tab *(age_group_size_factor*cluster_size_factor)
 normalized_cell_count <-  (0.5*age_group_size_factor+cluster_size_factor*0.5)
 #print(normalized_cell_count)
# normalized_cell_count <-  (age_group_size_factor*cluster_size_factor)
# scaled <-  scales::rescale(as.vector(normalized_cell_count))
# normalized_cell_count[] <-  matrix(scaled,nrow = nrow(normalized_cell_count),ncol = ncol(normalized_cell_count))
# normalized_cell_count <-  (age_group_size_factor/cluster_size_factor)
# normalized_cell_count <-  (cluster_size_factor/age_group_size_factor)
# normalized_cell_count <- age_group_size_factor
# 
 normalized_cell_count[is.nan(normalized_cell_count)] <- 0
 normalized_cell_count[is.infinite(normalized_cell_count)] <- 0

 
# d= 7
# (normalized_cell_count[1:5,d][normalized_cell_count[1:5,d] > sum(normalized_cell_count[1:5,d]) * 0.01]) #quantile(normalized_cell_count[1:5,d], 
# 
# (cluster_size_factor[1:5,d][cluster_size_factor[1:5,d] > sum(cluster_size_factor[1:5,d]) * 0.1]) #quantile(cluster_size_factor[1:5,d], 
# 
# 
# 
# 
# (age_group_size_factor[1:5,d][age_group_size_factor[1:5,d] > sum(age_group_size_factor[1:5,d]) * 0.05]) #quantile(age_group_size_factor[1:5,d], 
# 
# (age_group_cluster_tab[1:5,d][age_group_cluster_tab[1:5,d] > sum(age_group_cluster_tab[1:5,d]) * 0.15]) #quantile(age_group_cluster_tab[1:5,d], 
# (normalized_cell_count[1:5,d][normalized_cell_count[1:5,d] > sum(normalized_cell_count[1:5,d]) * 0.05]) #quantile(normalized_cell_count[1:5,d], 
# (normalized_cell_count[1:5,d][normalized_cell_count[1:5,d] > sum(normalized_cell_count[1:5,d]) * 0.2]) #quantile(normalized_cell_count[1:5,d], 
# 
# 
# 
# (age_group_size_factor[1:5,d][age_group_size_factor[1:5,d] > sum(age_group_size_factor[1:5,d]) * 0.08]) #quantile(age_group_size_factor[1:5,d], 

normalized_cell_count_drop_low <- normalized_cell_count
age_group_size_factor_drop_low <- age_group_size_factor
cluster_size_factor_drop_low <- cluster_size_factor
binary_mat <- age_group_size_factor
binary_mat[] <- 1


for (d in seq(colnames(normalized_cell_count))){
  age_zero <- names(normalized_cell_count[,d][normalized_cell_count[,d] < sum(normalized_cell_count[,d]) * 0.2]) #quantile(normalized_cell_count[1:5,d], 
  normalized_cell_count_drop_low[age_zero, d] <- 0
  age_group_size_factor_drop_low[age_zero, d] <- 0
  normalized_cell_count_drop_low[age_zero, d] <- 0
  binary_mat[age_zero, d] <- 0
  
}

#print(binary_mat)


# this is not good (cluster_size_factor[1:5,d][cluster_size_factor[1:5,d] > sum(cluster_size_factor[1:5,d]) * 0.05]) #quantile(cluster_size_factor[1:5,d], 





colData(cds)[[clustering_column]] <- droplevels(colData(cds)[[clustering_column]])
multi_cluster <- data.frame('09_11w'=numeric(0),'12_13w'=numeric(0),'14_16w'=numeric(0),'17_19w'=numeric(0),'20_22w'=numeric(0))
multi_cluster <- as.matrix(multi_cluster)
cluster_count <- c()
for (cluster in levels(colData(cds)[[clustering_column]])){
  
  cell_group_df <- data.frame(cell=colnames(cds[,colData(cds)[[clustering_column]] %in% cluster]) , group = cds[,colData(cds)[[clustering_column]] %in% cluster][[age_clustering_column]] ) #clusters(cds)
  mat_dspp <- aggregate_gene_expression(cds[,colData(cds)[[clustering_column]] %in% cluster],cell_group_df = cell_group_df, norm_method ="size_only" )
  row.names(mat_dspp) = rowData(cds[row.names(mat_dspp),])$gene_short_name
  
  # choose  at a time
  # selected_genes <- c("AMELX")
  #selected_genes <- c("AMBN")
  #selected_genes <- c("DSPP")
  #selected_genes <- c("SP6","DSPP","AMBN")
  
  multi_genes <-data.frame('09_11w'=numeric(0),'12_13w'=numeric(0),'14_16w'=numeric(0),'17_19w'=numeric(0),'20_22w'=numeric(0))
  multi_genes <- as.matrix(multi_genes)
  for (gene in selected_genes ){
    
    mat_selected <- mat_dspp[gene,]
    
    
    #create an empty matrix
    temp_m <- Matrix::Matrix(0, 1,length(levels(colData(cds)[[age_clustering_column]])))
    colnames(temp_m) <- levels(colData(cds)[[age_clustering_column]])
    row.names(temp_m) <- gene
    
    
    # this will write to only the clusters that express the genes, while the others remain zero
    temp_m[,match( names(mat_selected), colnames(temp_m)) ] <- mat_selected
    
    
    # temp_s <- as.data.frame(as.matrix(temp_m))/table(colData(cds)[[age_clustering_column]])
    #  temp_s <- as.data.frame(as.matrix(temp_s))/table(colData(cds)[[clustering_column]])[cluster]
    #  
    
    if (normalize){
    #temp_m <- as.data.frame(as.matrix(temp_m))/table(colData(cds)[[clustering_column]])[cluster]
    temp_m <- as.data.frame(as.matrix(temp_m))*normalized_cell_count[,cluster]
 
    }
    #temp_m <- as.data.frame(as.matrix(temp_m))/table(colData(cds)[[age_clustering_column]])
    #temp_s <- as.data.frame(as.matrix(temp_m))/age_group_cluster_tab[,cluster]
    
    if (trimm_low){
    temp_m <- as.data.frame(as.matrix(temp_m))*binary_mat[,cluster]
    }
    #temp_s <- as.data.frame(as.matrix(temp_m))/age_group_sum
    
    
    temp_s <- as.matrix(temp_m)
    #temp_s <- t(scale(t(temp_s))) 
    multi_genes <- rbind(multi_genes,temp_s)
  }
  multi_cluster <- rbind(multi_cluster, multi_genes)
  cluster_count <- c(cluster_count,rep(cluster,n_genes))
}
anno_df <- data.frame(names= row.names(multi_cluster), cluster = cluster_count)
anno_df$cluster <- factor(anno_df$cluster, levels =  unique(anno_df$cluster))
anno_df$names <- factor(anno_df$names, levels =  unique(anno_df$names))

multi_cluster[is.nan(multi_cluster)] <- 0
multi_cluster[is.infinite(multi_cluster)] <- 0
#row.names(multi_cluster) <- seq(anno_df$names)
colnames(multi_cluster) <- sub('.', '', colnames(multi_cluster)) # remove X


# 
# pdf("realtime_test.pdf")
# ComplexHeatmap::Heatmap(
#   name = "expression",
#   multi_cluster ,
#   cluster_columns = F,
#   cluster_rows = F,
#   show_row_name = T,
#   show_column_names = T,
#   column_names_side = "top",
#   column_names_rot = 45,
#   row_order = NULL,
#   # cell_fun = function(j, i, x, y, width, height, fill) {
#   #   grid.text(sprintf("%.1f", multi_cluster[i, j]), x, y, gp = gpar(fontsize = 10))
#   # },
#   
#   left_annotation  = rowAnnotation(foo= anno_block(gp = gpar(fill = randomcoloR::distinctColorPalette(length(unique(anno_df$cluster))) ), #randomcoloR::distinctColorPalette(length(unique(anno_df$cluster)))
#                                                    labels = unique(anno_df$cluster), 
#                                                    labels_gp = gpar(col = "black", fontsize = 10))),
#   row_split = factor(anno_df$cluster, levels =  unique(anno_df$cluster)),
#   row_title = " ",
#   use_raster = T,
#   show_heatmap_legend = T
# )
# h1
# dev.off()
# 
# 
# multi_cluster[is.nan(multi_cluster)] <- 0
# ComplexHeatmap::Heatmap(
#   name = "expression",
#   multi_cluster ,
#   cluster_columns = F,
#   cluster_rows = F,
#   show_row_name = T,
#   show_column_names = T,
#   column_names_side = "top",
#   column_names_rot = 45,
#   row_title = " ",
#   row_split = factor(anno_df$cluster, levels =  unique(anno_df$cluster)),
#   show_heatmap_legend = T
# )

multi_cluster_reordered <- multi_cluster[seq.int(1, length(anno_df$names), length(levels(anno_df$names)) ),]
for (gene_i in 2: n_genes){
  multi_cluster_reordered <- rbind(multi_cluster_reordered,multi_cluster[seq.int(gene_i, length(anno_df$names), length(levels(anno_df$names)) ),])
}


vec <- as.matrix(as.vector( multi_cluster_reordered[1:n_clusters,]))
start<- 1
for (gene_i in 2: n_genes){
vec <- cbind(vec, as.vector( multi_cluster_reordered[(n_clusters+start):(n_clusters*gene_i),]))
start <- n_clusters+start

}

vec <- scale(vec,center = TRUE, scale = TRUE)


ht_list = NULL  ## Heatmap(...) + NULL gives you a HeatmapList object
print("Draw heatmap")
for (gene_i in 1: n_genes){
  matrix_sub <- multi_cluster_reordered[1:n_clusters,]
  matrix_sub[] <- matrix(vec[,gene_i],nrow = n_clusters,ncol = n_age)
  ht_list = ht_list %v% ComplexHeatmap::Heatmap(
    name = selected_genes[gene_i],
    matrix_sub ,
    cluster_columns = F,
    cluster_rows = F,
    show_row_name = T,
    show_column_names = T,
    column_names_side = "top",
    column_names_rot = 90,
    column_names_gp = gpar(fontsize = 15,fontface = "bold"),
    row_order = NULL,
    row_labels = unique(anno_df$cluster),
    row_title = selected_genes[gene_i],
    row_title_gp = gpar(fontsize = 20),
    row_names_gp = gpar(fontsize = 15, fontface = "bold"),
    use_raster = T,
    show_heatmap_legend = F,
    col = circlize::colorRamp2(c( min(matrix_sub),0, max(matrix_sub)), colors = rev(RColorBrewer::brewer.pal(n = 3, name = "RdYlBu")))
                                                                                                                                                        
  )
    

}

#seq.int(from =2,length.out=2 , by=6)
#multi_cluster_reordered[(gene_i*n_clusters):(n_clusters*gene_i),] <- matrix(vec[,gene_i],nrow = n_clusters,ncol = n_age)
# multi_cluster_reordered[1:7,] <- matrix(vec,nrow = 7,ncol = 5)
# 
# vec <- scale(as.vector( multi_cluster_reordered[8:14,]))
# vec <- scales::rescale(vec, to=c(-2,2) )
# multi_cluster_reordered[8:14,] <- matrix(vec,nrow = 7,ncol = 5)
# 
# 
# vec <- scale(as.vector( multi_cluster_reordered[15:21,]))
# vec <- scales::rescale(vec, to=c(-2,2) )
# multi_cluster_reordered[15:21,] <- matrix(vec,nrow = 7,ncol = 5)
# 
# 
# 
# multi_cluster_reordered[1:7,] <- scales::rescale(multi_cluster_reordered[1:7,], to=c(-2,2) )
# multi_cluster_reordered[8:14,] <-scales::rescale(multi_cluster_reordered[8:14,], to=c(-2,2) )
# multi_cluster_reordered[15:21,] <- scales::rescale(multi_cluster_reordered[15:21,] , to=c(-2,2))
# 
# multi_cluster_reordered[1:7,] <- t(scale(t(multi_cluster_reordered[1:7,])))
# multi_cluster_reordered[1:7,] <- (scale((multi_cluster_reordered[1:7,])))
# multi_cluster_reordered[8:14,] <- t(scale(t(multi_cluster_reordered[8:14,])))
# multi_cluster_reordered[8:14,] <- (scale((multi_cluster_reordered[8:14,])))
# multi_cluster_reordered[15:21,] <- t(scale(t(multi_cluster_reordered[15:21,])))
# multi_cluster_reordered[15:21,] <- (scale((multi_cluster_reordered[15:21,])))
# 



# h1 <-  ComplexHeatmap::Heatmap(
#   name = "SP6",
#   multi_cluster_reordered[1:7,] ,
#   cluster_columns = F,
#   cluster_rows = F,
#   show_row_name = T,
#   show_column_names = T,
#   column_names_side = "top",
#   column_names_rot = 90,
#   row_order = NULL,
#   row_labels = unique(anno_df$cluster),
#   row_title = "SP6",
#   use_raster = T,
#   show_heatmap_legend = T,
#   col = circlize::colorRamp2(c( min(multi_cluster_reordered[1:7,]),0, max(multi_cluster_reordered[1:7,])), colors = colorRampPalette(rev(brewer.pal(n = 7, name =
#                                                                                                                                                       "RdYlBu")))(100))
# )
# 
# h2 <-  ComplexHeatmap::Heatmap(
#   name = "DSPP",
#   multi_cluster_reordered[8:14,]  ,
#   cluster_columns = F,
#   cluster_rows = F,
#   show_row_name = T,
#   show_column_names = T,
#   column_names_side = "top",
#   column_names_rot = 90,
#   row_order = NULL,
#   row_labels = unique(anno_df$cluster),
#   row_title = "DSPP",
#   use_raster = T,
#   show_heatmap_legend = T,
#   col = circlize::colorRamp2(c( min(multi_cluster_reordered[8:14,]),0, max(multi_cluster_reordered[8:14,])), c( "lightslateblue","white","red"))
# )
# 
# h3 <-  ComplexHeatmap::Heatmap(
#   name = "AMBN",
#   multi_cluster_reordered[15:21,] ,
#   cluster_columns = F,
#   cluster_rows = F,
#   show_row_name = T,
#   show_column_names = T,
#   column_names_side = "top",
#   column_names_rot = 90,
#   row_order = NULL,
#   row_labels = unique(anno_df$cluster),
#   row_title = "AMBN",
#   use_raster = T,
#   show_heatmap_legend = T,
#   col = circlize::colorRamp2(c( min(multi_cluster_reordered[15:21,]),0, max(multi_cluster_reordered[15:21,])), c( "lightslateblue","white","red"))
# )

pdf(paste("realtime_Heatmap_output","/",cds_name,"_",paste(selected_genes,collapse = "_"),"_realtime_Heatmap.pdf", sep = ""), 4,10)
draw(ht_list, padding = unit(c(2, 2, 2, 25), "mm")) 

col_fun = circlize::colorRamp2(c(-4, 0, 4), color= rev(RColorBrewer::brewer.pal(n = 3, name = "RdYlBu")))
lgd = Legend(col_fun = col_fun, title = "Expression", direction = "vertical",  title_position = "leftcenter-rot", title_gp = gpar(fontsize = 14, fontface = "bold"), at = c(-4, 0, 4), 
             labels = c("low", "med", "high"), labels_gp = gpar(col = "black", font = 4))
draw(lgd, x = unit(0.9, "npc"), y = unit(0.5, "npc"))

dev.off()
print("Done! results are saved in realtime_Heatmap_output folder")
end_time <- Sys.time()

end_time - start_time

}

#plot_heatmap_over_realtime(amelo_trajectory,"assigned_cell_type.4",age_clustering_column="age_group", selected_genes= c("SP6","DSPP","AMBN","PLOD2") )
