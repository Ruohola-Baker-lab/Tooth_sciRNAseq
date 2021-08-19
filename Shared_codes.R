library(monocle3)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)



###########
### to plot 3 or less genes with density and expression percent and log10(expression)
###########

genes_list <- c("DSPP","AMELX","PITX2")

# select time point
#time_point<- "20_22w" # ----- 12_13w, 17_18w, 20_22w

plot <-
  plot_cells(
    Incisors_20_cds,
    color_cells_by = "cluster",
    group_cells_by = "cluster",
    cell_size = 0.5,
    group_label_size = 5,
    alpha = 1,
    label_cell_groups = F,
    show_trajectory_graph = F,
    genes = genes_list,
    min_expr = 0.1
  ) + theme(aspect.ratio = 1)

# comment this line out if you want to plot all time points
#plot$data$value[plot$data$age_group != time_point] <- NA

plot2 <- plot
plot3 <- plot
plot_data <- plot$data

plot_data_subset <- plot_data[!is.na(plot_data$value), ]

for (g in genes_list){
  gene_subset<- plot_data_subset[plot_data_subset$feature_label ==g,]
  dplot <-
    densCols(gene_subset$data_dim_1,
             gene_subset$data_dim_2,
             colramp = colorRampPalette(c("black", "white")))
  gene_subset$dens <- col2rgb(dplot)[1, ] + 1L
  #gene_subset$value <- (log10(gene_subset$value)/ log10(max(gene_subset$value, na.rm = T))) * 100
  #gene_subset$value  <- quantile(gene_subset$value , c(0, 0.995))
  gene_subset$value  <- percent_rank(gene_subset$value)*100
  
  if (is.null(plot_data_subset$dens)) { 
    plot_data_subset$dens <- NA 
  }
  
  plot_data_subset[plot_data_subset$feature_label ==g, ] <- gene_subset
}

plot_data$dens <- NA
plot_data[!is.na(plot_data$value), ] <- plot_data_subset
plot_data$dens = (plot_data$dens  / 256) * 100



## Load color palate
jBrewColors <-
  c(
    "Reds",
    "Greens",
    "Blues",
    "Oranges",
    "Purples",
    "Reds",
    "Greens",
    "Blues",
    "Oranges",
    "Purples"
  )

#Mapping density
cols_r <- brewer.pal(n = 9, name = jBrewColors[1])[1:9]
plot <-
  plot + scale_colour_gradient(name = "Expression Density",
                                
                                low = "#e3e3e3", high = cols_r[9],
                                na.value = "#e3e3e3")
plot <-
  plot + geom_point(data = plot_data[order(plot_data$dens, na.last = F), ],
                    size = 0.7,
                    aes(x = data_dim_1, y = data_dim_2, color = dens))  + theme(aspect.ratio = 1)

#Mapping expression

cols_g <- brewer.pal(n = 9, name = jBrewColors[2])[1:8]
plot2 <-
  plot2 + scale_colour_gradient(
    name = "Expression percentile",
    low = "#000300", high = "#33d23b",
    na.value = "#e3e3e3", limits = c(0, 95.5), labels =c("0th",">99.5th"), breaks= c(0, 95.5)
  )
plot2 <-
  plot2 + geom_point(data = plot_data[order(plot_data$value, na.last = F), ],
                     size = 0.7,
                     aes(x = data_dim_1, y = data_dim_2, color = value))  + theme(aspect.ratio = 1)


# 10Log(expression)

plot_data3 <- plot3$data
plot_data3$value = log10(plot_data3$value)

cols_b <- brewer.pal(n = 9, name = jBrewColors[3])[1:8]
plot3 <-
  plot3 + scale_colour_gradientn(
    name = "log10(Expression)",
    colors = cols_b,
    na.value = "#e3e3e3"
  )
plot3 <-
  plot3 + geom_point(data = plot_data3[order(plot_data3$value, na.last = F), ],
                     size = 0.7,
                     aes(x = data_dim_1, y = data_dim_2, color = value))  + theme(aspect.ratio = 1)


grid.arrange(plot, plot2, ncol = 1)


############
######  Temp
############

plot <-
  plot_cells(
    Inciosrs_epi,
    color_cells_by = "cluster",
    group_cells_by = "cluster",
    cell_size = 0.5,
    group_label_size = 5,
    alpha = 1,
    label_cell_groups = F,
    show_trajectory_graph = F,
    genes = "AMELX",
    min_expr = 0.1
  ) + theme(aspect.ratio = 1)
plot_data <- plot$data
plot_data_subset <- plot_data[!is.na(plot_data$value), ]

## Use densCols() output to get density at each point
dplot <-
  densCols(plot_data_subset$data_dim_1,
           plot_data_subset$data_dim_2,
           colramp = colorRampPalette(c("black", "white")))
plot_data_subset$dens <- col2rgb(dplot)[1, ] + 1L
plot_data$dens <- NA
plot_data[!is.na(plot_data$value), ] <- plot_data_subset
plot_data$dens = (plot_data$dens  / 256) * 100
plot_data$value = (plot_data$value  / max(plot_data$value, na.rm = T)) * 100

## Load color palate
jBrewColors <-
  c(
    "Reds",
    "Greens",
    "Blues",
    "Oranges",
    "Purples",
    "Reds",
    "Greens",
    "Blues",
    "Oranges",
    "Purples"
  )

#Mapping density
cols_r <- brewer.pal(n = 9, name = jBrewColors[1])[1:8]
plot <-
  plot + scale_colour_gradientn(name = "Expression Density",
                                colors = cols_r,
                                na.value = "#e3e3e3")
plot <-
  plot + geom_point(data = plot_data[order(plot_data$dens, na.last = F), ],
                    size = 0.7,
                    aes(x = data_dim_1, y = data_dim_2, color = dens))  + theme(aspect.ratio = 1)

#Mapping expression
cols_g <- brewer.pal(n = 9, name = jBrewColors[2])[1:8]
plot2 <-
  plot + scale_colour_gradientn(
    name = "Expression score",
    colors = cols_g,
    na.value = "#e3e3e3",
  )
plot2 <-
  plot2 + geom_point(data = plot_data[order(plot_data$value, na.last = F), ],
                     size = 0.7,
                     aes(x = data_dim_1, y = data_dim_2, color = value))  + theme(aspect.ratio = 1)

grid.arrange(plot, plot2, ncol = 1)




###########
### to plot multiple genes
###########

# use this if your genes are in txt file
wnt_ligands <- read.delim("Wnt/Wnt_ligands.txt", header = F)
wnt_ligands <- as.character(wnt_ligands[,1])

Activated_canonical_wnt <-
  read.delim("Wnt/Wnt_postive_regulation_of_canonical.txt", header = F)
Activated_canonical_wnt <- as.character(Activated_canonical_wnt[, 2])
Activated_canonical_wnt <-
  Activated_canonical_wnt[!Activated_canonical_wnt %in% wnt_ligands]

# use this to add your genes manually 
wnt_ligands <- c("WNT1",	"WNT2",	"WNT2B",	"WNT3",	"WNT3A",	"WNT4",	"WNT5A",	"WNT5B",	"WNT6",	"WNT7A",
                 "WNT7B",	"WNT8A",	"WNT8B",	"WNT9A",	"WNT9B",	"WNT10A",	"WNT10B",	"WNT11",	"WNT16")

Activated_canonical_wnt <- c("DAPK3",	"ZBED3",	"PSMA4","PSMA3"	,"PSMA2",	"PSMA1",	"PSMC2", "SCEL",	"CDH3",
                             "LGR5",	"AMER1",	"EGFR",	"AXIN2",	"CTDNEP1",	"SFRP2",	"USP8",	"RNF220",	"URS0000527F89_9606",
                             "EGF",	"PSMA5",	"PSMB4",	"CSNK1D",	"PSMB9"	,"PSMB8",	"PSMD14",	"FRAT1",	"UBE2B","CSNK1G1",	"TRPM4",
                             "JUP",	"PSMD6",	"LGR4",	"NRARP",	"EDA",	"GID8",	"ROR2",	"PSMB11",	"LGR6",	"USP34",	"FGFR2",	"SFRP4",
                             "CSNK1G3",	"RSPO2",	"ASPM",	"LRRK2",	"CCAR2",	"FGF10"	,"PSME2",	"ILK"	,"CAPRIN2",	"GPRC5B",	"PSME3",
                             "CAV1",	"TBL1XR1",	"TBL1X",	"BAMBI",	"TMEM198",	"ZBED2",	"ADGRA2",	"BIRC8",	"RNF146",	"XIAP",	
                             "COL1A1",	"AXIN1",	"VPS35",	"SRC",	"PTK7",	"DACT1",	"NFKB1",	"FZD9",	"SEMA5A",	"CSNK1G2",
                             "RSPO3",	"PIN1",	"DAAM2",	"WLS",	"PSMA8","SFRP1",	"PSMC4",	"VCP",	"ARNTL",	"LYPD6",
                             "PSMD1",	"PSMD7",	"GPC3",	"PSMB7",	"YAP1",	"DKK2",	"PSMA7",	"PRDM15",	"PPM1A",	"PSMC5",	"PSMC1",
                             "UBR5",	"PSMA6",	"FAM53B",	"PSMF1",	"PSMB10",	"RSPO1",	"PSMD13",	"RBPJ",	"PSME1",	"RUVBL1",	"PSMB1",
                             "PSMB2",	"PSMB3",	"WNK1",	"USP47",	"CSNK1E",	"TNKS2",	"SOX4",	"TNKS",	"PSMC3",	"PLEKHA4",	"DDX3X",
                             "LRRK1",	"GSKIP",	"PSMD11",	"GPC5",	"PSMD8",	"PSMB5",	"PSMB6",	"PSMD4",	"PSME4",	"PSMD10",	"WNK2",	
                             "RSPO4",	"RECK",	"NLE1",	"SMURF2",	"PSMD5",	"JRK",	"DLX5",	"KANK1",	"FGF9",	"PSMD9",	"PSMD12",	"SULF2",
                             "PSMD2",	"URS00002C6949_9606",	"URS000024463E_9606",	"PSMD3",	"SMAD3",	"PSMC6"
)


control <- sample(fData(Inciosrs_epi)$id , size= length(wnt_ligands), replace =F)
control2 <- sample(fData(Inciosrs_epi)$id, size= length(Activated_canonical_wnt), replace =F)

genes1 <- as.data.frame(rowData(Inciosrs_epi)[rowData(Inciosrs_epi)$gene_short_name %in% wnt_ligands,])
genes2 <- as.data.frame(rowData(Inciosrs_epi)[rowData(Inciosrs_epi)$gene_short_name %in% Activated_canonical_wnt,])
control <- as.data.frame(rowData(Inciosrs_epi)[rowData(Inciosrs_epi)$id %in% control,])
control2 <- as.data.frame(rowData(Inciosrs_epi)[rowData(Inciosrs_epi)$id %in% control2,])


genes1$pathway <- "Wnt_ligands"
genes2$pathway <- "Activated_canonical_wnt"
control$pathway <- "Random_control-1"
control2$pathway <- "Random_control-2"

genes <- rbind(genes1, genes2, control, make.row.names=F)

plot <- plot_cells(Inciosrs_epi, color_cells_by="cluster", group_cells_by="cluster", cell_size=0.5, 
                   group_label_size= 5, alpha = 1, label_cell_groups = F, show_trajectory_graph = F, genes= genes[,c(1,3)], 
                   min_expr=0.1) + theme(aspect.ratio = 1)


plot2 <- plot
plot3 <- plot
plot_data <- plot$data

plot_data_subset <- plot_data[!is.na(plot_data$value), ]

for (g in unique(genes$pathway)){
  gene_subset<- plot_data_subset[plot_data_subset$feature_label ==g,]
  dplot <-
    densCols(gene_subset$data_dim_1,
             gene_subset$data_dim_2,
             colramp = colorRampPalette(c("black", "white")))
  gene_subset$dens <- col2rgb(dplot)[1, ] + 1L
  gene_subset$value <- (gene_subset$value)/ (max(gene_subset$value, na.rm = T)) * 100
  if (is.null(plot_data_subset$dens)) { 
    plot_data_subset$dens <- NA 
  }
  
  plot_data_subset[plot_data_subset$feature_label ==g, ] <- gene_subset
}

plot_data$dens <- NA
plot_data[!is.na(plot_data$value), ] <- plot_data_subset
plot_data$dens = (plot_data$dens  / 256) * 100
plot_data$average = ((plot_data$dens *0.25) + (plot_data$value * 0.75))


## Load color palate
jBrewColors <-
  c(
    "Reds",
    "Greens",
    "Blues",
    "Oranges",
    "Purples",
    "Reds",
    "Greens",
    "Blues",
    "Oranges",
    "Purples"
  )

#Mapping density
cols_r <- brewer.pal(n = 9, name = jBrewColors[1])[1:8]
plot <-
  plot + scale_colour_gradientn(name = "Expression Density",
                                colors = cols_r,
                                na.value = "#e3e3e3")
plot <-
  plot + geom_point(data = plot_data[order(plot_data$dens, na.last = F), ],
                    size = 0.7,
                    aes(x = data_dim_1, y = data_dim_2, color = dens))  + theme(aspect.ratio = 1)

#Mapping expression

cols_g <- brewer.pal(n = 9, name = jBrewColors[2])[1:8]
plot2 <-
  plot2 + scale_colour_gradientn(
    name = "Expression percent",
    colors = cols_g,
    na.value = "#e3e3e3",
  )
plot2 <-
  plot2 + geom_point(data = plot_data[order(plot_data$value, na.last = F), ],
                     size = 0.7,
                     aes(x = data_dim_1, y = data_dim_2, color = value))  + theme(aspect.ratio = 1)


# 10Log(expression)

#plot_data3 <- plot3$data
#plot_data3$value = log10(plot_data3$value)

cols_b <- brewer.pal(n = 9, name = jBrewColors[3])[1:8]
plot3 <-
  plot3 + scale_colour_gradientn(
    name = "Overlay averaged",
    colors = cols_b,
    na.value = "#e3e3e3",
  )
plot3 <-
  plot3 + geom_point(data = plot_data[order(plot_data$average, na.last = F), ],
                     size = 0.7,
                     aes(x = data_dim_1, y = data_dim_2, color = average))  + theme(aspect.ratio = 1)


grid.arrange(plot, plot2, plot3, ncol = 1)




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
