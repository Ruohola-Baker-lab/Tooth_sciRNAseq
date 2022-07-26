suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  library(viridis)
  library(ggnewscale)
})


w_directory = "/Users/maryregier/Desktop/HRB/RNAscope_ammar_sesha/final_image_set"
setwd(dir=w_directory)

key = c("VWDE","SALL1","FGF4","IGFBP5","DSPP","FGF10","PRRX1","FBN2","ENAM","PCDH7","SOX5","KRT5")

densities = F

### CHOOSE SAMPLE ----

samplename = "80d_B2"

## 117d_H3 
if (samplename == "117d_H3") {
  full_det = read.delim(paste0("/Users/maryregier/Desktop/HRB/RNAscope_ammar_sesha/final_image_set/117d_H3_central_incisor_section1/QuPath_117d_H3_ci_s1_full/117d_H3_full_detections_bc.txt"), header = TRUE)
  missing_probes = c(1,4)
  ptsize = 1.5
  time = 117
  xdim = 8986/8986
  ydim = 8986/8986
}

## 117d_K3 
if (samplename == "117d_K3") {
  full_det = read.delim(paste0("/Users/maryregier/Desktop/HRB/RNAscope_ammar_sesha/final_image_set/117d_K3_left_incisor_section1/QuPath_117d_K3/117d_K3_full_detections.txt"), header = TRUE)
  full_det$Centroid.Y.px = max(full_det$Centroid.Y.px) - full_det$Centroid.Y.px
  missing_probes = c(1)
  ptsize = 1.5
  time = 117
  xdim = 6816/7476
  ydim = 8136/7476
}

## 80d_E1 
if (samplename == "80d_E1") {
  full_det = read.delim(paste0("/Users/maryregier/Desktop/HRB/RNAscope_ammar_sesha/final_image_set/80d_E1_lower_Ant/QuPath_80d_E1/80d_E1_full_detections.txt"), header = TRUE)
  missing_probes = c(1)
  ptsize = 2.5
  time = 80
  xdim = 6221/6221
  ydim = 6221/6221
  
}

## 80d_B2 
if (samplename == "80d_B2") {
  full_det = read.delim(paste0("/Users/maryregier/Desktop/HRB/RNAscope_ammar_sesha/final_image_set/80d_B2_U.Incisor/80d_B2_UI/80d_B2_UI_full_detections.txt"), header = TRUE)
  full_det$Centroid.X.px = max(full_det$Centroid.X.px) - full_det$Centroid.X.px
  missing_probes = c(1)
  ptsize = 2.5
  time = 80
  xdim = 8180/8180
  ydim = 8180/8180
}



### PROCESS FULL DATASET ----

# pull coordinate column names (may be in px or um)
x_coord <- colnames(full_det)[6]
y_coord <- colnames(full_det)[7]

# number each cell
full_det$cellID = NA
cell_num = 1
for (i in 1:length(full_det[,1])) {
  if (full_det$Class[i] == ""){
    full_det$cellID[i] <- paste0("cell_",cell_num)
    cell_num = cell_num+1
  }
}

# remove rows for subcellular detections
det_cells_full = full_det %>% filter(!is.na(cellID))

# pull cell IDs
cells_full = as.data.frame(det_cells_full$cellID)
names(cells_full)[1] <- "cellID"

# make columns for probes
for (i in 1:length(key)){
  cells_full[1+i] <- NA
  names(cells_full)[1+i] <- key[i]
}


# place number of spots estimated for each probe in probe's column for each cell
for (i in 1:length(det_cells_full[,1])){
  for (k in 2:(length(key)+1)){
    if (k %in% missing_probes){
    } else {
      kcol = paste0("Subcellular..Channel.",k,"..Num.spots.estimated")
      cells_full[i,k] = det_cells_full[i,kcol]
    }
  }
}

# merge cells IDs and their number of spots estimated with their coordinates
positions_full <- dplyr::select(det_cells_full,c(length(det_cells_full[1,]),6,7))
cells_full <- merge(positions_full, cells_full, by = "cellID")
cells_full[cells_full == 0] <- NA

# mean center and unit variance scale by probe
for (channel_col in 4:length(cells_full[1,])){
  cells_full[,channel_col] <- scale(cells_full[,channel_col], scale = T, center = T)
}


### Plot by transcript ----

if (densities == T){
  
  for (k in c(1,2,3,4,5,6,7,8,9,10,11,12)){
    pos = cells_full
    column = key[k]
    pos_data = pos %>% filter(!is.na(pos[[column]]))
    print(ggplot() +
            geom_point(data = pos %>% filter(is.na(pos$KRT5)),  aes(x = .data[[x_coord]], y = .data[[y_coord]]), size = 0.5,  color = "gray90")+
            ylim(min(pos[3]),max(pos[3]))+
            geom_point(data = pos %>% filter(!is.na(pos$KRT5)),  aes(x = .data[[x_coord]], y = .data[[y_coord]]), size = 0.5,  color = "gray85")+
            geom_point(data = pos_data, aes(x = .data[[x_coord]], y = .data[[y_coord]], colour = .data[[column]]), size = 0.5) +
            theme_void() +
            ggtitle(key[k]) +
            theme(plot.title = element_text(hjust = 0.5)))
    ggsave(paste0(samplename,"_all_",key[k],".png"),
           width = xdim*6,
           height = ydim*6+2)
  }
  
}



### Epithelial-derived cell type analysis ----

alpha = 0.75


#### 117d - Epithelial ----

if (time == 117) {
  
  cells <- cells_full
  
# select probes to analyze
  subset = c(1,4,5,8,9,10,12)
  probeset <- cells[,c(1,2,3,subset+3)]
  
# set cell type clusters, alternative IDs, and color values  
  clusters = c( "SRI", "SRO", "OEE","PA", "IEE", "SII", "SIO",  "CL", "eAM", "sAM")
  noID = c("Ambiguous", "Unannotated")
  color_values = c( "SRI" = "#8b0000" , "SRO" = "#fbc02d","OEE" = "#a1887f", "IEE" = "#0000ff", "PA" = "#ff8f00", "SII" = "#4caf50","SIO" = "#4b0082",  "CL" = "#81d4f4", "eAM" = "#6495ed", "sAM" = "#ff69b4","Ambiguous" = "gray20", "Unannotated" = "gray90")
  
# use only KRT5 negative cells
  clusterset <- probeset %>% filter(!is.na(probeset$KRT5))
  
# make category columns for each probe  
  col_offset = length(probeset[1,])
  j=1
  for (i in subset){
    cat_label = paste0(key[i],"_category")
    clusterset[col_offset+j] <- 0
    names(clusterset)[col_offset+j] <- cat_label
    j = j+1
  }
  
# categorize expression of each probe using alpha to set percentile cutoff 
  for (j in 1:length(subset)) {
    high = as.numeric(quantile(clusterset[[j+3]],prob=alpha, na.rm = T))
    for (k in 1:length(clusterset[,1])) {
      if (is.na(clusterset[k,j+3])) {}
      else {
        if (clusterset[k,j+3] >= high ) {
          clusterset[k,col_offset+j] = "High"
        }
        if (!is.na(clusterset[k,j+3] ) & (clusterset[k,j+3] < high)) {
          clusterset[k,col_offset+j] = "Low"
        }
      }
    }
  }

 
# set logic table for cell type-specific combinations of expression levels
  
  PA = list("KRT5" = c("Low"), "IGFBP5" = c("Low"), "PCDH7" = c(0, "Low"), "FBN2" = c(0, "Low"), "VWDE" = c("High"))
  IEE = list("KRT5" = c("Low"), "IGFBP5" = c(0), "PCDH7" = c(0, "Low"), "FBN2" = c("Low"), "VWDE" = c("Low"))
  SII = list("KRT5" = c("High"), "IGFBP5" = c("Low"), "PCDH7" = c(0, "Low"), "FBN2" = c("High"), "VWDE" = c("High"))
  SIO = list("KRT5" = c("High"), "IGFBP5" = c("Low","High"), "PCDH7" = c(0, "Low"), "FBN2" = c("Low", "High"), "VWDE" = c(0, "Low"))
  OEE = list("KRT5" = c("Low","High"), "IGFBP5" = c("High"), "PCDH7" = c(0, "Low", "High"), "FBN2" = c(0, "Low"), "VWDE" = c(0))
  CL = list("KRT5" = c("Low","High"), "IGFBP5" = c(0), "PCDH7" = c(0, "Low"), "FBN2" = c("High"), "VWDE" = c("Low"), "DSPP" = c(0), "ENAM" = c(0))
  SRI = list("KRT5" = c("Low"), "IGFBP5" = c(0, "Low"), "PCDH7" = c(0,"Low","High"), "FBN2" = c(0), "VWDE" = c("Low"))
  SRO = list("KRT5" = c("Low"), "IGFBP5" = c(0, "Low"), "PCDH7" = c(0,"Low","High"), "FBN2" = c(0), "VWDE" = c(0))
  sAM = list("KRT5" = c("Low","High"), "DSPP" = c(0, "Low"), "ENAM" = c("High"))
  eAM = list("KRT5" = c("Low","High"), "DSPP" = c("High"), "ENAM" = c(0, "Low"))
  
  clusterset$Cluster = 0
  clusterset$Combo = 0
  
  
# Determine which criteria each cell meets
  
  for (i in 1:length(clusterset[,1])){
    
    if ((clusterset$PCDH7_category[i] %in% PA$PCDH7) & (clusterset$FBN2_category[i] %in% PA$FBN2) & (clusterset$KRT5_category[i] %in% PA$KRT5) & (clusterset$VWDE_category[i] %in% PA$VWDE) & (clusterset$IGFBP5_category[i] %in% PA$IGFBP5)){ 
      clusterset$Cluster[i] = "PA"
      clusterset$Combo[i] = ("PA")
    }
    
    if ((clusterset$PCDH7_category[i] %in% IEE$PCDH7) & (clusterset$FBN2_category[i] %in% IEE$FBN2) & (clusterset$KRT5_category[i] %in% IEE$KRT5) & (clusterset$VWDE_category[i] %in% IEE$VWDE) & (clusterset$IGFBP5_category[i] %in% IEE$IGFBP5)){ 
      if (clusterset$Cluster[i] != 0){
        clusterset$Cluster[i] = "Ambiguous" 
        clusterset$Combo[i] = paste0(clusterset$Combo[i],"+IEE")
      }
      else {clusterset$Cluster[i] = "IEE"
      clusterset$Combo[i] = ("IEE")
      }
    }
    
    if ((clusterset$PCDH7_category[i] %in% SII$PCDH7) & (clusterset$FBN2_category[i] %in% SII$FBN2) & (clusterset$KRT5_category[i] %in% SII$KRT5) & (clusterset$VWDE_category[i] %in% SII$VWDE) & (clusterset$IGFBP5_category[i] %in% SII$IGFBP5)){ 
      if (clusterset$Cluster[i] != 0){
        clusterset$Cluster[i] = "Ambiguous"
        clusterset$Combo[i] = paste0(clusterset$Combo[i],"+SII")
      }
      else {clusterset$Cluster[i] = "SII"
      clusterset$Combo[i] = "SII"
      }
    }
    
    if ((clusterset$PCDH7_category[i] %in% SIO$PCDH7) & (clusterset$FBN2_category[i] %in% SIO$FBN2) & (clusterset$KRT5_category[i] %in% SIO$KRT5) & (clusterset$VWDE_category[i] %in% SIO$VWDE) & (clusterset$IGFBP5_category[i] %in% SIO$IGFBP5)){ 
      if (clusterset$Cluster[i] != 0){
        clusterset$Cluster[i] = "Ambiguous" 
        clusterset$Combo[i] = paste0(clusterset$Combo[i],"+SIO")
      }
      else {clusterset$Cluster[i] = "SIO"
      clusterset$Combo[i] = "SIO"
      }
    }
    
    if ((clusterset$PCDH7_category[i] %in% OEE$PCDH7) & (clusterset$FBN2_category[i] %in% OEE$FBN2) & (clusterset$KRT5_category[i] %in% OEE$KRT5) & (clusterset$VWDE_category[i] %in% OEE$VWDE) & (clusterset$IGFBP5_category[i] %in% OEE$IGFBP5)){ 
      if (clusterset$Cluster[i] != 0){
        clusterset$Cluster[i] = "Ambiguous"
        clusterset$Combo[i] = paste0(clusterset$Combo[i],"+OEE")
      }
      else {clusterset$Cluster[i] = "OEE"
      clusterset$Combo[i] = "OEE"
      }
    }
    
    if ((clusterset$PCDH7_category[i] %in% CL$PCDH7) & (clusterset$FBN2_category[i] %in% CL$FBN2) & (clusterset$KRT5_category[i] %in% CL$KRT5) & (clusterset$VWDE_category[i] %in% CL$VWDE) & (clusterset$IGFBP5_category[i] %in% CL$IGFBP5) & (clusterset$DSPP_category[i] %in% CL$DSPP) & (clusterset$ENAM_category[i] %in% CL$ENAM) ){ 
      if (clusterset$Cluster[i] != 0){
        clusterset$Cluster[i] = "Ambiguous"
        clusterset$Combo[i] = paste0(clusterset$Combo[i],"+CL")
      }
      else {clusterset$Cluster[i] = "CL"
      clusterset$Combo[i] = "CL"
      }
    }
    
    if ((clusterset$PCDH7_category[i] %in% SRI$PCDH7) & (clusterset$FBN2_category[i] %in% SRI$FBN2) & (clusterset$KRT5_category[i] %in% SRI$KRT5) & (clusterset$VWDE_category[i] %in% SRI$VWDE) & (clusterset$IGFBP5_category[i] %in% SRI$IGFBP5)){ 
      if (clusterset$Cluster[i] != 0){
        clusterset$Cluster[i] = "Ambiguous"
        clusterset$Combo[i] = paste0(clusterset$Combo[i],"+SRI")
      }
      else {clusterset$Cluster[i] = "SRI"
      clusterset$Combo[i] = "SRI"
      }
    }
    
    if ((clusterset$PCDH7_category[i] %in% SRO$PCDH7) & (clusterset$FBN2_category[i] %in% SRO$FBN2) & (clusterset$KRT5_category[i] %in% SRO$KRT5) & (clusterset$VWDE_category[i] %in% SRO$VWDE) & (clusterset$IGFBP5_category[i] %in% SRO$IGFBP5)){ 
      if (clusterset$Cluster[i] != 0){
        clusterset$Cluster[i] = "Ambiguous"
        clusterset$Combo[i] = paste0(clusterset$Combo[i],"+SRO")
      }
      else {clusterset$Cluster[i] = "SRO"
      clusterset$Combo[i] = "SRO"
      }
    }
    
    if ((clusterset$DSPP_category[i] %in% eAM$DSPP) & (clusterset$ENAM_category[i] %in% eAM$ENAM) & (clusterset$KRT5_category[i] %in% eAM$KRT5)){
      if (clusterset$Cluster[i] != 0){
        clusterset$Cluster[i] = "Ambiguous"
        clusterset$Combo[i] = paste0(clusterset$Combo[i],"+eAM")
      }
      else {clusterset$Cluster[i] = "eAM"
      clusterset$Combo[i] = "eAM"
      }
    }
    
    if ((clusterset$DSPP_category[i] %in% sAM$DSPP) & (clusterset$ENAM_category[i] %in% sAM$ENAM) & (clusterset$KRT5_category[i] %in% sAM$KRT5)){
      if (clusterset$Cluster[i] != 0){
        clusterset$Cluster[i] = "Ambiguous"
        clusterset$Combo[i] = paste0(clusterset$Combo[i],"+sAM")
      }
      else {clusterset$Cluster[i] = "sAM"
      clusterset$Combo[i] = "sAM"
      }
    }
  }

# Designate cells that don't meet criteria sets as unannotated
    
  for (i in 1:length(clusterset[,1])){
    if (clusterset$Cluster[i] == 0 & clusterset$KRT5_category[i] !=0){
      clusterset$Cluster[i] = "Unannotated"
    }
  }
  
 # Plot KRT5- and KRT5+ with each cell type separately
  for (k in 1:length(clusters)){
    pos = clusterset %>% filter(Cluster == clusters[k])
    
    ggplot() +
      geom_point(data = cells %>% filter(is.na(cells$KRT5)),  aes(x = .data[[x_coord]], y = .data[[y_coord]]), size = ptsize/4,  color = "gray95")+
      geom_point(data = cells %>% filter(!is.na(cells$KRT5)),  aes(x = .data[[x_coord]], y = .data[[y_coord]]), size = ptsize/4,  color = "gray90")+
      ylim(min(cells[3]),max(cells[3]))+
      geom_point(data = pos, aes(x = .data[[x_coord]], y = .data[[y_coord]], color = Cluster), size = ptsize/4) +
      scale_color_manual(values = color_values) + 
      theme_void() + 
      ggtitle(clusters[k]) +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
    ggsave(paste0(samplename,"_analysis3_",clusters[k],".png"),
           width = xdim*6,
           height = ydim*6+.5)
  }
  
# Plot, KRT5-, KRT5+ Unannotated, KRT5+ Ambiguous, and all cell types
  pos = clusterset %>% filter(Cluster %in% c(clusters,"Ambiguous"))
  
  
  print(ggplot() +
          geom_point(data = cells %>% filter(is.na(cells$KRT5)),  aes(x = .data[[x_coord]], y = .data[[y_coord]]), size = ptsize/4,  color = "gray95")+
          geom_point(data = cells %>% filter(!is.na(cells$KRT5)),  aes(x = .data[[x_coord]], y = .data[[y_coord]]), size = ptsize/4,  color = "gray90")+
          new_scale_color() +
          geom_point(data = pos, aes(x = .data[[x_coord]], y = .data[[y_coord]], color = Cluster), size = ptsize/4) +
          scale_color_manual(values = color_values) + 
          theme_void())
  ggsave(paste0(samplename,"_analysis3_all.png"),
         width = xdim*6+1,
         height = ydim*6)
  
  
# Plot KRT5- and KRT5+ with specific combinations of cell types  
  print(ggplot() +
          geom_point(data = cells %>% filter(is.na(cells$KRT5)),  aes(x = .data[[x_coord]], y = .data[[y_coord]]), size = ptsize/4,  color = "gray95")+
          geom_point(data = cells %>% filter(!is.na(cells$KRT5)),  aes(x = .data[[x_coord]], y = .data[[y_coord]]), size = ptsize/4,  color = "gray90")+
          new_scale_color() +
          geom_point(data = pos %>% filter(pos$Cluster %in% c("PA", "SII", "SIO")), aes(x = .data[[x_coord]], y = .data[[y_coord]], color = Cluster), size = ptsize/4) +
          scale_color_manual(values = color_values) + 
          theme_void())
  ggsave(paste0(samplename,"_analysis3_PA_SI.png"),
         width = xdim*6+1,
         height = ydim*6)
  
  print(ggplot() +
          geom_point(data = cells %>% filter(is.na(cells$KRT5)),  aes(x = .data[[x_coord]], y = .data[[y_coord]]), size = ptsize/4,  color = "gray95")+
          geom_point(data = cells %>% filter(!is.na(cells$KRT5)),  aes(x = .data[[x_coord]], y = .data[[y_coord]]), size = ptsize/4,  color = "gray90")+
          new_scale_color() +
          geom_point(data = pos %>% filter(pos$Cluster %in% c("SRO", "SRI", "OEE", "CL")), aes(x = .data[[x_coord]], y = .data[[y_coord]], color = Cluster), size = ptsize/4) +
          scale_color_manual(values = color_values) + 
          theme_void())
  ggsave(paste0(samplename,"_analysis3_SR_OEE_CL.png"),
         width = xdim*6+1,
         height = ydim*6)
  
  
  print(ggplot() +
          geom_point(data = cells %>% filter(is.na(cells$KRT5)),  aes(x = .data[[x_coord]], y = .data[[y_coord]]), size = ptsize/4,  color = "gray95")+
          geom_point(data = cells %>% filter(!is.na(cells$KRT5)),  aes(x = .data[[x_coord]], y = .data[[y_coord]]), size = ptsize/4,  color = "gray90")+
          new_scale_color() +
          geom_point(data = pos %>% filter(pos$Cluster %in% c("PA", "IEE", "eAM", "sAM")), aes(x = .data[[x_coord]], y = .data[[y_coord]], color = Cluster), size = ptsize/4) +
          scale_color_manual(values = color_values) + 
          theme_void())
  ggsave(paste0(samplename,"_analysis3_PA_AM.png"),
         width = xdim*6+1,
         height = ydim*6)
}


#### 80d - Epithelial ----

if (time == 80){
  
  cells <- cells_full
  
# select probes to analyze
  subset = c(1,3,4,5,8,9,10,12)
  probeset <- cells[,c(1,2,3,subset+3)]
  
# set cell type clusters, alternative IDs, and color values    
  clusters = c( "EK", "IEE", "SII", "OEE", "CL", "SRI")
  noID = c("Ambiguous", "Unannotated")
  color_values = c("EK" = "#9575cd", "IEE" = "#0000ff", "SII" = "#4caf50", "OEE" = "#a1887f", "CL" = "#81d4f4", "SRI" = "#8b0000", "Ambiguous" = "gray20", "Unannotated" = "gray90")
  noID_values = c("Ambiguous" = "gray60", "Unannotated" = "gray90")
  
  
# use only KRT5 negative cells 
  clusterset <- probeset %>% filter(probeset$KRT5 != 0)
  clusterset[clusterset == 0] <- NA
  
  
# make category columns for each probe    
  col_offset = length(probeset[1,])
  j=1
  for (i in subset){
    cat_label = paste0(key[i],"_category")
    clusterset[col_offset+j] <- 0
    names(clusterset)[col_offset+j] <- cat_label
    j = j+1
    
  }
  
# categorize expression of each probe using alpha to set percentile cutoff 
  for (j in 1:length(subset)) {
    high = as.numeric(quantile(clusterset[[j+3]],prob=alpha, na.rm = T))
    for (k in 1:length(clusterset[,1])) {
      if (is.na(clusterset[k,j+3])) {}
      else {
        if (clusterset[k,j+3] >= high ) {
          clusterset[k,col_offset+j] = "High"
        }
        if ((clusterset[k,j+3] > 0 ) & (clusterset[k,j+3] < high)) {
          clusterset[k,col_offset+j] = "Low"
        }
      }
    }
  }
  
# set logic table for cell type-specific combinations of expression levels  
  
  EK = list(PCDH7 = c(0, "Low"), FBN2 = c(0,"Low", "High"), KRT5 = c("Low", "High"), VWDE = c(0, "Low"), IGFBP5 = c(0), FGF4 = c("High", "Low"))
  IEE = list(PCDH7 = c(0, "Low"), FBN2 = c(0,"Low"), KRT5 = c("Low", "High"), VWDE = c("Low", "High"), IGFBP5 = c(0, "Low"), FGF4 = c(0, "Low"))
  SRI = list(PCDH7 = c(0, "Low", "High"), FBN2 = c(0), KRT5 = c("Low"), VWDE = c(0,"Low"), IGFBP5 = c(0, "Low"), FGF4 = c(0, "Low"))
  OEE = list(PCDH7 = c(0, "Low", "High"), FBN2 = c(0,"Low", "High"), KRT5 = c("Low","High"), VWDE = c(0), IGFBP5 = c("Low", "High"), FGF4 = c(0, "Low"))
  CL = list(PCDH7 = c(0, "Low", "High"), FBN2 = c("Low","High"), KRT5 = c("Low", "High"), VWDE = c(0,"Low"), IGFBP5 = c(0, "Low"), FGF4 = c(0, "Low"))
  SII = list(PCDH7 = c(0, "Low"), FBN2 = c("Low", "High"), KRT5 = c("Low","High"), VWDE = c("Low", "High"), IGFBP5 = c(0, "Low"), FGF4 = c(0, "Low"))
  
  
  clusterset$Cluster = 0
  clusterset$Combo = 0
 
  
# Determine which criteria each cell meets
  
  for (i in 1:length(clusterset[,1])){
    if ((clusterset$PCDH7_category[i] %in% EK$PCDH7) & (clusterset$FBN2_category[i] %in% EK$FBN2) & (clusterset$KRT5_category[i] %in% EK$KRT5) & (clusterset$VWDE_category[i] %in% EK$VWDE) & (clusterset$IGFBP5_category[i] %in% EK$IGFBP5) & (clusterset$FGF4_category[i] %in% EK$FGF4)){ 
      if (clusterset$Cluster[i] != 0){
        clusterset$Cluster[i] = "Ambiguous"
        clusterset$Combo[i] = paste0(clusterset$Combo[i],"+EK")
      }
      else {clusterset$Cluster[i] = "EK"
      clusterset$Combo[i] = "EK"}
    }
    if ((clusterset$PCDH7_category[i] %in% IEE$PCDH7) & (clusterset$FBN2_category[i] %in% IEE$FBN2) & (clusterset$KRT5_category[i] %in% IEE$KRT5) & (clusterset$VWDE_category[i] %in% IEE$VWDE) & (clusterset$IGFBP5_category[i] %in% IEE$IGFBP5) & (clusterset$FGF4_category[i] %in% IEE$FGF4)){ 
      if (clusterset$Cluster[i] != 0){
        clusterset$Cluster[i] = "Ambiguous" 
        clusterset$Combo[i] = paste0(clusterset$Combo[i],"+IEE")
      }
      else {clusterset$Cluster[i] = "IEE"
      clusterset$Combo[i] = "IEE"}
    }
    if ((clusterset$PCDH7_category[i] %in% SRI$PCDH7) & (clusterset$FBN2_category[i] %in% SRI$FBN2) & (clusterset$KRT5_category[i] %in% SRI$KRT5) & (clusterset$VWDE_category[i] %in% SRI$VWDE) & (clusterset$IGFBP5_category[i] %in% SRI$IGFBP5) & (clusterset$FGF4_category[i] %in% SRI$FGF4)){ 
      if (clusterset$Cluster[i] != 0){
        clusterset$Cluster[i] = "Ambiguous" 
        clusterset$Combo[i] = paste0(clusterset$Combo[i],"+SRI")
      }
      else {clusterset$Cluster[i] = "SRI"
      clusterset$Combo[i] = "SRI"}
    }
    if ((clusterset$PCDH7_category[i] %in% OEE$PCDH7) & (clusterset$FBN2_category[i] %in% OEE$FBN2) & (clusterset$KRT5_category[i] %in% OEE$KRT5) & (clusterset$VWDE_category[i] %in% OEE$VWDE) & (clusterset$IGFBP5_category[i] %in% OEE$IGFBP5)& (clusterset$FGF4_category[i] %in% OEE$FGF4)){ 
      if (clusterset$Cluster[i] != 0){
        clusterset$Cluster[i] = "Ambiguous" 
        clusterset$Combo[i] = paste0(clusterset$Combo[i],"+OEE")
      }
      else {clusterset$Cluster[i] = "OEE"
      clusterset$Combo[i] = "OEE"}
    }
    if ((clusterset$PCDH7_category[i] %in% CL$PCDH7) & (clusterset$FBN2_category[i] %in% CL$FBN2) & (clusterset$KRT5_category[i] %in% CL$KRT5) & (clusterset$VWDE_category[i] %in% CL$VWDE) & (clusterset$IGFBP5_category[i] %in% CL$IGFBP5) & (clusterset$FGF4_category[i] %in% CL$FGF4)){ 
      if (clusterset$Cluster[i] != 0){
        clusterset$Cluster[i] = "Ambiguous"
        clusterset$Combo[i] = paste0(clusterset$Combo[i],"+CL")
      }
      else {clusterset$Cluster[i] = "CL"
      clusterset$Combo[i] = "CL"}
    }
    if ((clusterset$PCDH7_category[i] %in% SII$PCDH7) & (clusterset$FBN2_category[i] %in% SII$FBN2) & (clusterset$KRT5_category[i] %in% SII$KRT5) & (clusterset$VWDE_category[i] %in% SII$VWDE) & (clusterset$IGFBP5_category[i] %in% SII$IGFBP5) & (clusterset$FGF4_category[i] %in% SII$FGF4)){ 
      if (clusterset$Cluster[i] != 0){
        clusterset$Cluster[i] = "Ambiguous" 
        clusterset$Combo[i] = paste0(clusterset$Combo[i],"+SII")
      }
      else {clusterset$Cluster[i] = "SII"
      clusterset$Combo[i] = "SII"}
    }
  }
  
  # Designate cells that don't meet criteria sets as unannotated  
  for (i in 1:length(clusterset[,1])){
    if (clusterset$Cluster[i] == 0 & clusterset$KRT5_category[i] !=0){
      clusterset$Cluster[i] = "Unannotated"
    }
    
  }
  
  # Plot KRT5- and KRT5+ with each cell type separately
  
  for (k in 1:length(clusters)){
    pos = clusterset %>% filter(Cluster == clusters[k])
    
    ggplot() +
      geom_point(data = probeset %>% filter(is.na(probeset$KRT5)),  aes(x = .data[[x_coord]], y = .data[[y_coord]]), size = ptsize/2,  color = "gray95")+
      geom_point(data = probeset %>% filter(!is.na(probeset$KRT5)),  aes(x = .data[[x_coord]], y = .data[[y_coord]]), size = ptsize/2,  color = "gray90")+
      new_scale_color() +
      ylim(min(cells[3]),max(cells[3]))+
      geom_point(data = pos, aes(x = .data[[x_coord]], y = .data[[y_coord]], color = Cluster), size = ptsize/2) +
      scale_color_manual(values = color_values) + 
      theme_void() +
      ggtitle(clusters[k]) +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
    

      ggsave(paste0(samplename,"_analysis3_80d",clusters[k],".png"),
             width = xdim*6,
             height = ydim*6+.5)
  }

   
  # Plot, KRT5-, KRT5+ Unannotated, KRT5+ Ambiguous, and all cell types
  pos = clusterset %>% filter(Cluster %in% c(clusters,"Ambiguous"))
  
  
  
  print(ggplot() +
          geom_point(data = probeset %>% filter(is.na(probeset$KRT5)),  aes(x = .data[[x_coord]], y = .data[[y_coord]]), size = ptsize/2,  color = "gray95")+
          geom_point(data = probeset %>% filter(!is.na(probeset$KRT5)),  aes(x = .data[[x_coord]], y = .data[[y_coord]]), size = ptsize/2,  color = "gray90")+
          new_scale_color() +
          geom_point(data = pos, aes(x = .data[[x_coord]], y = .data[[y_coord]], color = Cluster), size = ptsize/2) +
          scale_color_manual(values = color_values) +
          theme_void()
  )
  
    
    ggsave(paste0(samplename,"_80d_analysis3_all.png"),
           width = xdim*6+1,
           height = ydim*6)

}
