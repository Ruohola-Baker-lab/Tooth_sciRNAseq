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

samplename = "80d_E1"

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


### Mesenchyme-derived cell type analysis ----

# Set alpha 
a = .75


#### 117d - Mesenchymal -----

if (time == 117) {
  
  cells <- cells_full
  
# select probes to analyze
  subset = c(2,4,5,6,7,8,11,12)
  
# set colors for cluster cell types
  color_vector = c("#F8766D","#B79F00","#00BA38","#00BFC4","#619CFF","#F564E3")
  clusters = c("DP","POB","DEM","DF","SOB","OB","Epithelial") 
  
  # Select probes in set  
  probeset <- cells[,c(1,subset+3)]
  probeset <- filter(probeset, is.na(probeset$KRT5))
  probeset <- merge(probeset, positions_full, by = "cellID")

# Make columns for cell 
  col_offset = length(subset)+3
  for (i in 1:length(clusters)){
    probeset[col_offset+i] <- NA
    names(probeset)[col_offset+i] <- clusters[i]
  }
  
# Score each cell for its cell type marker combo expression
  for (i in 1:length(probeset[,1])){
    if ( !is.na(probeset$FGF10[i]) & !is.na(probeset$SOX5[i]) & !is.na(probeset$SALL1[i]) & is.na(probeset$KRT5[i])) {
      probeset$DP[i] = (probeset$SOX5[i]+probeset$SALL1[i]+probeset$FGF10[i])/3
    }
    if ((!is.na(probeset$DSPP[i])) & (is.na(probeset$KRT5[i]))) {
      probeset$OB[i] = probeset$DSPP[i]
    }
    if ((!is.na(probeset$IGFBP5[i])) & (!is.na(probeset$SALL1[i]))  & (is.na(probeset$KRT5[i]))) {
      probeset$SOB[i] = (probeset$IGFBP5[i] + probeset$SALL1[i])/2
    }
    if ((!is.na(probeset$IGFBP5[i])) & (!is.na(probeset$PRRX1[i])) & (is.na(probeset$SALL1[i])) & (is.na(probeset$KRT5[i]))) {
      probeset$DF[i] = (probeset$IGFBP5[i]+probeset$PRRX1[i])/2
    }
    if ((!is.na(probeset$FBN2[i])) & (!is.na(probeset$SALL1[i]))  & (is.na(probeset$KRT5[i]))) {
      probeset$POB[i] = (probeset$FBN2[i]+probeset$SALL1[i])/2
    }
    if ((!is.na(probeset$PRRX1[i])) & (is.na(probeset$IGFBP5[i])) & (is.na(probeset$SALL1[i])) & (is.na(probeset$KRT5[i]))) {
      probeset$DEM[i] = probeset$PRRX1[i]
    }
    if (!is.na(probeset$KRT5[i])) {
      probeset$Epithelial[i] = 1
    }
  }
 
# Remove the scores for below alpha percentile for each cell type
  probeset$DP[probeset$DP < quantile(probeset$DP,prob=a, na.rm=TRUE)] <- NA
  probeset$OB[probeset$OB < quantile(probeset$OB,prob=a, na.rm=TRUE)] <- NA
  probeset$SOB[probeset$SOB < quantile(probeset$SOB,prob=a, na.rm=TRUE)] <- NA
  probeset$DF[probeset$DF < quantile(probeset$DF,prob=a, na.rm=TRUE)] <- NA
  probeset$POB[probeset$POB < quantile(probeset$POB,prob=a, na.rm=TRUE)] <- NA
  probeset$DEM[probeset$DEM < quantile(probeset$DEM,prob=a, na.rm=TRUE)] <- NA
  
  
# Select only cells scoring in top alpha percentile for at least one cell type 
  pos = probeset %>% filter(!is.na(probeset$DP)  | !is.na(probeset$OB)  | !is.na(probeset$SOB)  | !is.na(probeset$POB)  | !is.na(probeset$DF)  | !is.na(probeset$DEM))
  
# set color ranges
  probecolors_high = color_vector
  probecolors_low = color_vector
  for (i in 1:length(probecolors_high)){
    gradfunc = colorRampPalette(c("white", probecolors_high[i]))
    probecolors_low[i] = gradfunc(6)[3]
  }
  

#Plot all cell types and unannotated KRT- cells, and KRT5+ cells 
  print(ggplot() +
          geom_point(data = cells_full %>% filter(is.na(cells_full$KRT5)),  aes(x = .data[[x_coord]], y = .data[[y_coord]]), size = ptsize,  color = "gray95")+
          geom_point(data = cells_full %>% filter(!is.na(cells_full$KRT5)),  aes(x = .data[[x_coord]], y = .data[[y_coord]]), size = ptsize,  color = "gray90")+
          geom_point(data = pos %>% filter(!is.na(DP)), aes(x = .data[[x_coord]], y = .data[[y_coord]], color = DP), size = ptsize) +
          scale_color_gradient(low = probecolors_low[1], high = probecolors_high[1], limits = c(min(pos$DP),max(pos$DP)),  oob = scales::squish) +
          new_scale_color() +
          geom_point(data = pos %>% filter(!is.na(POB)), aes(x = .data[[x_coord]], y = .data[[y_coord]],color = POB), size = ptsize) +
          scale_color_gradient(low = probecolors_low[2], high = probecolors_high[2], limits = c(min(pos$POB),max(pos$POB)),  oob = scales::squish) +
          new_scale_color() +
          geom_point(data = pos %>% filter(!is.na(DEM)), aes(x = .data[[x_coord]], y = .data[[y_coord]],color = DEM), size = ptsize) +
          scale_color_gradient(low = probecolors_low[3], high = probecolors_high[3], limits = c(min(pos$DEM),max(pos$DEM)),  oob = scales::squish) +
          new_scale_color() +
          geom_point(data = pos %>% filter(!is.na(DF)), aes(x = .data[[x_coord]], y = .data[[y_coord]],color = DF), size = ptsize) +
          scale_color_gradient(low = probecolors_low[4], high = probecolors_high[4], limits = c(min(pos$DF),max(pos$DF)),  oob = scales::squish) +
          new_scale_color() +
          geom_point(data = pos %>% filter(!is.na(SOB)), aes(x = .data[[x_coord]], y = .data[[y_coord]],color = SOB), size = ptsize) +
          scale_color_gradient(low = probecolors_low[5], high = probecolors_high[5], limits = c(min(pos$SOB),max(pos$SOB)),  oob = scales::squish) +
          new_scale_color() +
          geom_point(data = pos %>% filter(!is.na(OB)), aes(x = .data[[x_coord]], y = .data[[y_coord]] ,color = OB), size = ptsize) +
          scale_color_gradient(low = probecolors_low[6], high = probecolors_high[6], limits = c(min(pos$OB),max(pos$OB)),  oob = scales::squish) +
          new_scale_color() + 
          ylim(min(cells_full[3]),max(cells_full[3]))+
          theme_void()+
          theme(legend.position="none"))
  ggsave(paste0(samplename,"_117d_OB.png"),
         width = xdim*8,
         height = ydim*8)
  
  
#Plot each cell type individually
  
  for (k in 1:length(clusters)){
    ggplot() +
      geom_point(data = cells_full %>% filter(is.na(cells_full$KRT5)),  aes(x = .data[[x_coord]], y = .data[[y_coord]]), size = ptsize/2,  color = "gray95")+
      geom_point(data = cells_full %>% filter(!is.na(cells_full$KRT5)),  aes(x = .data[[x_coord]], y = .data[[y_coord]]), size = ptsize/2,  color = "gray90")+
      geom_point(data = pos %>% filter(!is.na(.data[[clusters[k]]])), aes(x = .data[[x_coord]], y = .data[[y_coord]] ,color = .data[[clusters[k]]]), size = ptsize/2) +
      scale_color_gradient(low = probecolors_low[k], high = probecolors_high[k], limits = c(min(pos[[clusters[k]]]),max(pos[[clusters[k]]])),  oob = scales::squish) +
      theme_void()+
      theme(legend.position="bottom")
    ggsave(paste0(samplename,"_117d_",clusters[k],"_OB.png"),
           width = xdim*6,
           height = ydim*6+.5)
  }
}

#### 80d - Mesenchyme ----

if (time == 80) {
  
  cells <- cells_full
  
  # select probes to analyze
  subset = c(2,4,5,6,7,8,11,12)
  
  # set colors for cluster cell types 
  color_vector = c("#F8766D","#00BA38","#00BFC4")
  clusters = c("DP","DEM","DF","Epithelial")
  

  # Select probes in set  
  probeset <- cells[,c(1,subset+3)]
  probeset <- filter(probeset, is.na(probeset$KRT5))
  probeset <- merge(probeset, positions_full, by = "cellID")
  
  # Make columns for cell 
  col_offset = length(subset)+3
  for (i in 1:length(clusters)){
    probeset[col_offset+i] <- NA
    names(probeset)[col_offset+i] <- clusters[i]
  }
  
  
  # Score each cell for its cell type marker combo expression
  for (i in 1:length(probeset[,1])){
    if ( !is.na(probeset$SOX5[i]) & !is.na(probeset$SALL1[i] ) & !is.na(probeset$FGF10[i]) & (is.na(probeset$KRT5[i]))) {
      probeset$DP[i] = (probeset$SOX5[i]+probeset$SALL1[i]+probeset$FGF10[i])/3
    }
    
    if ((!is.na(probeset$IGFBP5[i])) & (!is.na(probeset$PRRX1[i])) & (is.na(probeset$DSPP[i])) & (is.na(probeset$KRT5[i]))) {
      probeset$DF[i] = (probeset$IGFBP5[i]+probeset$PRRX1[i])/2
    }
    
    if ((!is.na(probeset$PRRX1[i])) & (is.na(probeset$IGFBP5[i])) & (is.na(probeset$SALL1[i])) & (is.na(probeset$KRT5[i]))) {
      probeset$DEM[i] = probeset$PRRX1[i]
    }
    if (!is.na(probeset$KRT5[i])) {
      probeset$Epithelial[i] = 1
    }
  }
  
  # Remove the scores for below alpha percentile for each cell type 
  probeset$DP[probeset$DP < quantile(probeset$DP,prob=a, na.rm=TRUE)] <- NA
  probeset$DF[probeset$DF < quantile(probeset$DF,prob=a, na.rm=TRUE)] <- NA
  probeset$DEM[probeset$DEM < quantile(probeset$DEM,prob=a, na.rm=TRUE)] <- NA

  # Select only cells scoring in top alpha percentile for at least one cell type   
  pos = probeset %>% filter(!is.na(probeset$DP)  |  !is.na(probeset$DF)  | !is.na(probeset$DEM))
  
  # set color ranges
  probecolors_high = color_vector
  probecolors_low = color_vector
  for (i in 1:length(probecolors_high)){
    gradfunc = colorRampPalette(c("white", probecolors_high[i]))
    probecolors_low[i] = gradfunc(6)[3]
  }
  
  #Plot all cell types and unannotated KRT- cells, and KRT5+ cells 
  print(ggplot() +
          geom_point(data = cells_full %>% filter(is.na(cells_full$KRT5)),  aes(x = .data[[x_coord]], y = .data[[y_coord]]), size = ptsize,  color = "gray95")+
          geom_point(data = cells_full %>% filter(!is.na(cells_full$KRT5)),  aes(x = .data[[x_coord]], y = .data[[y_coord]]), size = ptsize,  color = "gray90")+
          geom_point(data = pos %>% filter(!is.na(DP)), aes(x = .data[[x_coord]], y = .data[[y_coord]], color = DP), size = ptsize) +
          scale_color_gradient(low = probecolors_low[1], high = probecolors_high[1], limits = c(min(pos$DP),max(pos$DP)),  oob = scales::squish) +
          new_scale_color() +
          geom_point(data = pos %>% filter(!is.na(DEM)), aes(x = .data[[x_coord]], y = .data[[y_coord]],color = DEM), size = ptsize) +
          scale_color_gradient(low = probecolors_low[2], high = probecolors_high[2], limits = c(min(pos$DEM),max(pos$DEM)),  oob = scales::squish) +
          new_scale_color() +
          geom_point(data = pos %>% filter(!is.na(DF)), aes(x = .data[[x_coord]], y = .data[[y_coord]],color = DF), size = ptsize) +
          scale_color_gradient(low = probecolors_low[3], high = probecolors_high[3], limits = c(min(pos$DF),max(pos$DF)),  oob = scales::squish) +
          new_scale_color() +
          ylim(min(cells_full[3]),max(cells_full[3]))+
          theme_void()+
          theme(legend.position="bottom"))
  
    ggsave(paste0(samplename,"_80d_OB.png"),
           width = xdim*8,
           height = ydim*8+0.5)
  
  
  
  
  
  
  for (k in 1:length(clusters)){
    ggplot() +
      geom_point(data = cells_full %>% filter(is.na(cells_full$KRT5)),  aes(x = .data[[x_coord]], y = .data[[y_coord]]), size = ptsize*3/4,  color = "gray95")+
      geom_point(data = cells_full %>% filter(!is.na(cells_full$KRT5)),  aes(x = .data[[x_coord]], y = .data[[y_coord]]), size = ptsize*3/4,  color = "gray90")+
      geom_point(data = pos %>% filter(!is.na(.data[[clusters[k]]])), aes(x = .data[[x_coord]], y = .data[[y_coord]] ,color = .data[[clusters[k]]]), size = ptsize*3/4) +
      scale_color_gradient(low = probecolors_low[k], high = probecolors_high[k], limits = c(min(pos[[clusters[k]]]),max(pos[[clusters[k]]])),  oob =scales::squish)+
      theme_void()+
      theme(legend.position="bottom")

      ggsave(paste0(samplename,"_80d_",clusters[k],"_OB.png"),
             width = xdim*6,
             height = ydim*6+.5)
    
    
  }
  
  
}