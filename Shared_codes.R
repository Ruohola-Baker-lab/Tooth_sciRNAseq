library(monocle3)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)


# density_by time Point #####

if (!requireNamespace("reshape2", quietly = TRUE))
  install.packages('reshape2')
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(monocle3)
library(ggplot2)
library(dplyr)

plot_density <- function(cds, groups_column, nrow = NULL,
                         ncol = NULL, x_reverse = FALSE, y_reverse =FALSE){

assertthat::assert_that(!is.null(cds), msg = "cds object doesn't exist. Check the spelling!")
assertthat::assert_that(!is.null(groups_column),msg = "please add groups_column name")
assertthat::assert_that(!is.null(colData(cds)[[groups_column]]), msg = "groups_column doesn't exist in cds. Check the spelling!")

#updating_multi plot to facets####
plot <-
  plot_cells(
    cds,
    color_cells_by = groups_column ,
    group_cells_by = "cluster",
    cell_size = 0.5,
    group_label_size = 5,
    alpha = 1,
    label_cell_groups = F,
    show_trajectory_graph = F,
    min_expr = 0.1
  ) 

plot_data <- plot$data


if (is.null(levels(colData(cds)[[groups_column]]))){
list_of_grouping <- unique( as.character(colData(cds)[[groups_column]]))

}else list_of_grouping <- levels(colData(cds)[[groups_column]])

for (g in list_of_grouping){
  plot_data_subset<- plot_data[plot_data$cell_color ==g,]
  dplot <-
    densCols(plot_data_subset$data_dim_1,
             plot_data_subset$data_dim_2,
             colramp = colorRampPalette(c("black", "white")))
  plot_data_subset[[g]] <- col2rgb(dplot)[1, ] + 1L
  plot_data_subset[[g]]  = (plot_data_subset[[g]] / 256) * 100
  
  if (is.null(plot_data[[g]])) { 
    plot_data[[g]] <- NA 
  }
  
  plot_data[plot_data$cell_color ==g, ] <- plot_data_subset
}


## Load color palate
jBrewColors <- c("Reds", "Greens","Blues","Oranges", "Purples")


#Mapping density
plot_data_melted <- plot_data[c("data_dim_1", "data_dim_2","cell", list_of_grouping)]
plot_data_melted <- reshape2::melt(plot_data_melted, id.vars = c("data_dim_1", "data_dim_2","cell"), variable.name = "grouping", value.name = "grouping_density")
plot_data_melted$grouping_density <- as.numeric(plot_data_melted$grouping_density)
plot_data_melted$cell_color <- as.numeric(c(1:length(plot_data_melted$grouping_density)))

plot$data <- plot_data_melted[order(plot_data_melted["grouping"],plot_data_melted["grouping_density"], na.last = F), ]
plot[["layers"]] <- NULL
plot[["mapping"]] <- NULL
plot[["guides"]] <- NULL
cols_r <- RColorBrewer::brewer.pal(n = 9, name = jBrewColors[1])[1:9]
plot <-
  plot + scale_colour_gradient(name = "Density per group",
                               low =  cols_r[2], high = cols_r[9],
                               na.value = "#e3e3e3")
plot <- plot + geom_point(
  data = plot_data_melted[order(plot_data_melted["grouping"],plot_data_melted["grouping_density"], na.last = F), ],
  size = 0.5,
  aes(
    x = data_dim_1,
    y = data_dim_2,
    color = grouping_density
  )
)  + coord_fixed(ratio = 1)  + facet_wrap("grouping",ncol=ncol, nrow=nrow)+theme(
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.title=element_text(size=14))+theme(strip.text.x=element_text(size = 6,face = "bold"),legend.title = element_text(size=14, face="bold"))

if (x_reverse== T & y_reverse == F) {
  plot + scale_x_reverse()
}
else if (y_reverse == T & x_reverse == F) {
  plot + scale_y_reverse()
}

else if (x_reverse & y_reverse) {
  plot + scale_x_reverse() + scale_y_reverse()
}
else {
  plot 
}

}
#scale_x_reverse()

#' Duplicates data to create additional facet
#' @param df a dataframe
#' @param col the name of facet column
#'  
CreateAllFacet <- function(df, col){
  df$facet <- df[[col]]
  temp <- df
  temp$facet <- "all"
  merged <-rbind(temp, df)
  
  # ensure the facet value is a factor
  merged[[col]] <- as.factor(merged[[col]])
  
  return(merged)
}

plot_data_melted <- CreateAllFacet(plot_data_melted, "grouping")



#### liger density

plot_liger_density <- function(liger_obj, nrow = NULL,
                         ncol = NULL, x_reverse = FALSE, y_reverse =FALSE){
  
  assertthat::assert_that(!is.null(liger_obj), msg = "liger object doesn't exist. Check the spelling!")

  #updating_multi plot to facets####
  all.plots <- plotByDatasetAndCluster(liger_obj, axis.labels = c('UMAP 1', 'UMAP 2'), return.plots = T)
  all.plots[[1]] + all.plots[[2]]+ coord_fixed(ratio = 1)
  
  
  plot <- all.plots[[1]]
  plot_data <- plot$data
  
  

  list_of_grouping <- unique(plot_data$Dataset)
  
  for (g in list_of_grouping){
    plot_data_subset<- plot_data[plot_data$Dataset ==g,]
    dplot <-
      densCols(plot_data_subset$tsne1,
               plot_data_subset$tsne2,
               colramp = colorRampPalette(c("black", "white")))
    plot_data_subset[[g]] <- col2rgb(dplot)[1, ] + 1L
    plot_data_subset[[g]]  = (plot_data_subset[[g]] / 256) * 100
    
    if (is.null(plot_data[[g]])) { 
      plot_data[[g]] <- NA 
    }
    
    plot_data[plot_data$Dataset ==g, ] <- plot_data_subset
  }
  
  
  ## Load color palate
  jBrewColors <- c("Reds", "Greens","Blues","Oranges", "Purples")
  
  
  #Mapping density
  plot_data_melted <- plot_data[c("tsne1", "tsne2","Cluster", list_of_grouping)]
  plot_data_melted <- reshape2::melt(plot_data_melted, id.vars = c("tsne1", "tsne2","Cluster"), variable.name = "grouping", value.name = "grouping_density")
  plot_data_melted$grouping_density <- as.numeric(plot_data_melted$grouping_density)
  plot_data_melted$cell_color <- as.numeric(c(1:length(plot_data_melted$grouping_density)))
  
  plot$data <- plot_data_melted[order(plot_data_melted["grouping"],plot_data_melted["grouping_density"], na.last = F), ]
  plot[["layers"]] <- NULL
  plot[["mapping"]] <- NULL
  plot[["guides"]] <- NULL
  cols_r <- RColorBrewer::brewer.pal(n = 9, name = jBrewColors[1])[1:9]
  plot <-
    plot + scale_colour_gradient(name = "Density per group",
                                 low =  cols_r[2], high = cols_r[9],
                                 na.value = "#e3e3e3")
  plot <- plot + geom_point(
    data = plot_data_melted[order(plot_data_melted["grouping"],plot_data_melted["grouping_density"], na.last = F), ],
    size = 0.5,
    aes(
      x = tsne1,
      y = tsne2,
      color = grouping_density
    )
  )  + coord_fixed(ratio = 1)  + facet_wrap("grouping", ncol=ncol, nrow=nrow)+theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title=element_text(size=12))+theme(strip.text.x=element_text(size = 30))
  
  if (x_reverse== TRUE & y_reverse == FALSE) {
    plot + scale_x_reverse()
  }
  else if (y_reverse == T & x_reverse == F) {
    plot + scale_y_reverse()
  }
  
  else if (x_reverse & y_reverse) {
    plot + scale_x_reverse() + scale_y_reverse()
  }
  else  {
    plot 
  }
}



## change color function #########
change_colours <- function(p, palette) {
  n <- nlevels(p$data[[deparse(p$mapping$group)]])
  tryCatch(as.character(palette), 
           error=function(e) stop('palette should be a vector of colours', call.=FALSE))
  if(n > length(palette)) stop('Not enough colours in palette.')
  pal <- function(n) palette[seq_len(n)]
  p +  discrete_scale('colour', 'foo', pal)
}




#### Seurat density



plot_seurat_density <- function(Seurat_obj,groups_column, nrow = NULL,
                               ncol = NULL, x_reverse = FALSE, y_reverse =FALSE, drop= NULL, cluster_order = NULL){
  
  if (!requireNamespace("reshape2", quietly = TRUE))
    install.packages('reshape2')
  library(reshape2)
  library(ggplot2)
  library(RColorBrewer)
  library(tidyr)
  library(monocle3)
  library(ggplot2)
  library(dplyr)
  library(Seurat)
  
  assertthat::assert_that(!is.null(Seurat_obj), msg = "liger object doesn't exist. Check the spelling!")
  assertthat::assert_that(!is.null(groups_column),msg = "please add groups_column name")
  assertthat::assert_that(!is.null(Seurat_obj[[groups_column]]), msg = "groups_column doesn't exist in cds. Check the spelling!")
  
  #updating_multi plot to facets####
  

  plot <- DimPlot(refquery, group.by = groups_column, shuffle = F, cols = c("blue","red","green","orange","purple"), pt.size= 1,combine = F)
  
  plot_data <- plot[[1]]$data
  
  
  
  
  if (is.null(levels(Seurat_obj[[groups_column]]))){
    list_of_grouping <- as.factor(unlist((unique(Seurat_obj[[groups_column]]))))
    
  }else list_of_grouping <- levels(Seurat_obj[[groups_column]])
  
  for (g in list_of_grouping){
    plot_data[[g]] <- NA 
  }
  
  if  (!is.null(drop)){
    list_of_grouping <- list_of_grouping[!list_of_grouping %in% drop]
  }
  
  if  (!is.null(cluster_order)){
    list_of_grouping <- droplevels(factor(cluster_order,levels = cluster_order))
  }
  
  for (g in list_of_grouping){
    plot_data_subset<- plot_data[plot_data[[groups_column]] ==g,]
    dplot <-
      densCols(plot_data_subset$UMAP_1,
               plot_data_subset$UMAP_2,
               bandwidth=1,
               colramp = grDevices::colorRampPalette(c("black", "white")))
    plot_data_subset[[g]] <- col2rgb(dplot)[1, ] + 1L
    plot_data_subset[[g]]  = (plot_data_subset[[g]] / 256) * 100
    
    
    plot_data[plot_data[[groups_column]] ==g, ] <- plot_data_subset
  }
  
  
  
  ## Load color palate
  jBrewColors <- c("Reds", "Greens","Blues","Oranges", "Purples")
  
  
  #Mapping density
  plot_data_melted <- plot_data[c("UMAP_1", "UMAP_2",groups_column, as.character(list_of_grouping))]
  plot_data_melted <- reshape2::melt(plot_data_melted, id.vars = c("UMAP_1", "UMAP_2",groups_column), variable.name = "grouping", value.name = "grouping_density")
  plot_data_melted$grouping_density <- as.numeric(plot_data_melted$grouping_density)
  plot_data_melted[[groups_column]] <- as.numeric(c(1:length(plot_data_melted$grouping_density)))
  
  
  plot[[1]]$data <- plot_data_melted[order(plot_data_melted["grouping"],plot_data_melted["grouping_density"], na.last = F), ]
  plot[[1]][["layers"]] <- NULL
  plot[[1]][["mapping"]] <- NULL
  plot[[1]][["guides"]] <- NULL
  cols_r <- RColorBrewer::brewer.pal(n = 9, name = jBrewColors[1])[1:9]
  plot[[1]] <-
    plot[[1]] + scale_colour_gradient(name = "Density per group",
                                 low =  cols_r[2], high = cols_r[9],
                                 na.value = "#e3e3e3")
  plot[[1]] <- plot[[1]] + geom_point(
    data = plot_data_melted[order(plot_data_melted["grouping"],plot_data_melted["grouping_density"], na.last = F), ],
    size = 0.5,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      color = grouping_density
    )
  )  + coord_fixed(ratio = 1)  + facet_wrap("grouping", ncol=ncol, nrow=nrow)+theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title=element_text(size=12))+theme(strip.text.x=element_text(size = 30))
  
  
  
  
  if (x_reverse== TRUE & y_reverse == FALSE) {
    plot + scale_x_reverse()
  }
  else if (y_reverse == T & x_reverse == F) {
    plot + scale_y_reverse()
  }
  
  else if (x_reverse & y_reverse) {
    plot + scale_x_reverse() + scale_y_reverse()
  }
  else  {
    plot 
  }
}
