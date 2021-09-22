### REALTIME HEATMAP#####


# load the customized heatmap function.
source("R/functions.R")

# load your cell dataset object (cds)
amelo_odonto_cds <- readRDS("sample_datasets/ameloblast_odontoblast_lineage_cds.rds")

#  only supply your cds,  a clustering_column & age_group_column
# you can specify the markers genes by selected_genes = c("DSPP", "AMBN")

plot_heatmap_over_realtime(
  cds= amelo_odonto_cds,
  clustering_column = "assigned_cell_type",
  age_clustering_column = "age_group",
  selected_genes = c("SP6", "DSPP", "AMBN", "AMELX"),
  normalize = T,   # Toggle normalization
  trimm_low = T   # Trimm low cells that can be considered background and missclustered. In another words remove expression in time points that arenot likely to be real
)
