setwd("~/data_sci_RNA_seq")

#devtools::install_github("sqjin/CellChat")
library(tidyr)
library(monocle3)
library(ggplot2)
library(dplyr)
library(Matrix.utils)

pacman::p_unlock()
#devtools::install_github("sqjin/CellChat")
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(igraph)
options(stringsAsFactors = FALSE)

mixed_cds <- readRDS("/home/ubuntu/data_sci_RNA_seq/mixed_cds_simp.rds")

colData(mixed_cds)$assigned_cell_type.simp.2 <- dplyr::recode_factor(colData(mixed_cds)$assigned_cell_type.simp,
                                                              'OE'='OE',
                                                              'DE-Prog'='DE-Prog',
                                                              'EK'='IK/EK',
                                                              'OEE'='OEE',
                                                              'SR-1'="SR-1",
                                                              'SR-2'="SR-2",
                                                              'SI-1'='SI-1',
                                                              'SI-2'='SI-2',
                                                              'SI-3'='SI-3',
                                                            'PA-1' = 'PA',
                                                            'PA-2' = 'PA',
                                                            'AM-1' = 'Am',
                                                            'AM-2' = 'Am',
                                                            'D.mes' = 'D.mes',
                                                            'Pre.O' = 'Pre.O',
                                                            'OD-1' = 'OD-1',
                                                            'OD-2' = 'OD-2')

colData(mixed_cds)$assigned_cell_type.simp.2 <- dplyr::recode_factor(colData(mixed_cds)$assigned_cell_type.simp,
                                                              'OE'='OE',
                                                              'DE-Prog'='DE-Prog',
                                                              'EK'='IK/EK',
                                                              'OEE'='OEE',
                                                              'SR-1'="SR-1",
                                                              'SR-2'="SR-2",
                                                              'SI-1'='SI-1',
                                                              'SI-2'='SI-2',
                                                              'SI-3'='SI-3',
                                                            'PA-1' = 'PA',
                                                            'PA-2' = 'PA',
                                                            'AM-1' = 'Am',
                                                            'AM-2' = 'Am',
                                                            'D.mes' = 'D.mes',
                                                            'Pre.O' = 'Pre.O',
                                                            'OD-1' = 'OD-1',
                                                            'OD-2' = 'OD-2')

colData(mixed_cds)$assigned_cell_type.simp.3 <- dplyr::recode_factor(colData(mixed_cds)$assigned_cell_type.simp,
                                                              'OE'='OE',
                                                              'DE-Prog'='DE-Prog',
                                                              'EK'='IK/EK',
                                                              'OEE'='OEE',
                                                              'SR-1'="SR-1",
                                                              'SR-2'="SR-2",
                                                              'SI-1'='SI-1',
                                                              'SI-2'='SI-2',
                                                              'SI-3'='SI-3',
                                                              'PA-1' = 'PA',
                                                              'PA-2' = 'PA',
                                                              'AM-1' = 'Am',
                                                              'AM-2' = 'Am',
                                                              'D.mes' = 'D.mes',
                                                              'Pre.O' = 'Pre.O',
                                                              'OD-1' = 'Od',
                                                              'OD-2' = 'Od')

saveRDS(mixed_cds, "/home/ubuntu/data_sci_RNA_seq/mixed_cds_simp3.rds")
#levels(mixed_cds$age_group)[levels(mixed_cds$age_group)=="9_11w"] <- "09_11w"

in_groups <- list('09_11w' = c("OE" ,"DE-Prog",  "EK", "D.mes"),
                  '12_13w' = c("DE-Prog", "EK", "SR-1" , "D.mes"),
                  '14_16w' = c("EK", "OEE" , "SI-1", "SI-2",  "D.mes", "Pre.O"),
                  '17_19w' = c( "SR-2" , "SI-2", "PA-1",  "Pre.O", "OD-1"),
                  '20_22w' = c("SI-3", "PA-2", "AM-1", "AM-2",  "OD-2"))


in_groups <- list('09_11w' = c("OE" ,"DE-Prog",  "IK/EK", "D.mes"),
                  '12_13w' = c("DE-Prog", "IK/EK", "SR-1" , "D.mes"),
                  '14_16w' = c("IK/EK", "OEE" , "SI-1", "SI-2",  "D.mes", "Pre.O"),
                  '17_19w' = c( "SR-2" , "SI-2", "PA",  "Pre.O", "OD-1"),
                  '20_22w' = c("SI-3", "PA", "Am",  "OD-2"))

in_groups <- list('09_11w' = c("OE" ,"DE-Prog",  "IK/EK", "D.mes"),
                  '12_13w' = c("DE-Prog", "IK/EK", "SR-1" , "D.mes"),
                  '14_16w' = c("IK/EK", "OEE" , "SI-1", "SI-2",  "D.mes", "Pre.O"),
                  '17_19w' = c( "SR-2" , "SI-2", "PA",  "Pre.O", "Od"),
                  '20_22w' = c("SI-3", "PA", "Am",  "Od"))

#levels(colData(mixed_cds)$assigned_cell_type.simp.2)
#age = "09_11w"
#age = "17_19w"
#age = "20_22w"
#cellchat <- readRDS("/home/ubuntu/data_sci_RNA_seq/09_11w/cellchat_09_11w.rds")
for (age in levels(colData(mixed_cds)$age_group)) {
  setwd("~/data_sci_RNA_seq")
  
  if (age == "09_11w"){
    temp_cds <-  mixed_cds[, !(colData(mixed_cds)$age_group  != age & !colData(mixed_cds)$assigned_cell_type.simp.2 %in% c("OE","IK/EK"))]
  }
  if (age == "20_22w"){
    temp_cds <-  mixed_cds[, !(colData(mixed_cds)$age_group  != age & !colData(mixed_cds)$assigned_cell_type.simp.2 %in% c("Am","SI-3"))]
  }

  if (age != "09_11w" |  age != "20_22w"){
    temp_cds <-  mixed_cds[,colData(mixed_cds)$age_group  == age]
  }
  
  #temp_cds <-  mixed_cds
  #table( colData(temp_cds)$assigned_cell_type.simp.2)
  #included <- names(table( colData(temp_cds)$assigned_cell_type.simp.2)[table( colData(temp_cds)$assigned_cell_type.simp.2) > 20])
  included <- in_groups[age]
  temp_cds <- temp_cds[, colData(temp_cds)$assigned_cell_type.simp.2 %in% unlist(included)]
  

  data.input = normalized_counts(temp_cds, norm_method = "log") # normalized data matrix
  row.names(data.input) = rowData(temp_cds)$gene_short_name
  identity = data.frame(group = droplevels(colData(temp_cds)$assigned_cell_type.simp.2), row.names = names(colData(temp_cds)$assigned_cell_type.simp.2)) # create a dataframe consisting of the cell labels
  unique(identity$group) # check the cell labels
  data.input = normalized_counts(temp_cds, norm_method = "log") # normalized data matrix
  row.names(data.input) = rowData(temp_cds)$gene_short_name
  identity = data.frame(group = droplevels(colData(temp_cds)$assigned_cell_type.simp.2), row.names = names(colData(temp_cds)$assigned_cell_type.simp.2)) # create a dataframe consisting of the cell labels
  unique(identity$group) # check the cell labels
  cellchat <- createCellChat(object = data.input, meta = identity,  group.by = "group")
  groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
  CellChatDB <- CellChatDB.human # use CellChatDB.human if running on human data
  CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling")) # use Secreted Signaling for cell-cell communication analysis
  CellChatDB.use <- subsetDB(CellChatDB, search = c( "FGF inhibition receptor", "INSULIN inhibition receptor", 
                                                        "HH inhibition receptor", "ANGPT inhibition receptor 1", "ANGPT inhibition receptor 2"), key = "co_I_receptor")
  CellChatDB.use <- subsetDB(CellChatDB, search = c("TGFb antagonist", "BMP antagonist", "NODAL antagonist", 
                                                    "ACTIVIN antagonist", "WNT antagonist", "PDGF antagonist", "IGF antagonist", 
                                                    "IL antagonist", "GHRL antagonist", "POMC antagonist"), key = "antagonist")
  
  # CellChatDB.use <- CellChatDB # simply use the default CellChatDB
  cellchat@DB <- CellChatDB.use # set the used database in the object
  # use all CellChatDB for cell-cell communication analysis
  
  cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
  parallel::detectCores()
  cores_num <- parallel::detectCores()
  cores = ifelse(cores_num>1,cores_num-1,1)
  future::plan("multiprocess", workers = cores) # do parallel
  cellchat <- identifyOverExpressedGenes(cellchat) # if error, check cellchat@data.signaling
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)
  
  cellchat <- computeCommunProb(cellchat, type = "triMean", trim = 0.05, population.size = T ) #truncatedMean, triMean
  #if (age == "09_11w" | age == "20_22w"){  cellchat <- computeCommunProb(cellchat, type = "triMean", trim = 0.05, population.size = F )}
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  
  cellchat <- aggregateNet(cellchat)
  
  data.dir <- age  # Create a directory to save figures
  #dir.create(data.dir)
  setwd(data.dir)
  
  groupSize <- as.numeric(table(cellchat@idents))
  pdf(paste("antagonist_netCircle_weight_",age,".pdf", sep = ""))
  print(netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength",vertex.label.cex =3, arrow.size = 0.6))
  dev.off()
  pdf(paste("antagonist_netCircle_count_",age,".pdf", sep = ""))
  print(netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions",vertex.label.cex =3, arrow.size = 0.6 ))
  dev.off()
  # pdf(paste("antagonist_netCircle_both_",age,".pdf", sep = ""), 10,8)
  # par(mfrow = c(1,2), xpd=TRUE)
  # print(netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions",vertex.label.cex =3, arrow.size = 0.6 ))
  # print(netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength",vertex.label.cex =3, arrow.size = 0.6))
  # dev.off()
  # To examine the signaling sent from each cell group
  # mat <- cellchat@net$weight
  # par(mfrow = c(3,4), xpd=TRUE)
  # pdf("CellChat.pdf")
  # for (i in 1:nrow(mat)) {
  #   mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  #   mat2[i, ] <- mat[i, ]
  #   netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  # }
  # dev.off()
  # 
  
  
  # pathways.show <- c("TGFb") 
  # levels(cellchat@idents)
  # vertex.receiver = c(1,2) # a numeric vector
  # # Hierarchy plot
  # netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, vertex.size = groupSize)
  # # Circle plot
  # par(mfrow=c(1,1))
  # netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle", vertex.size = groupSize)
  # # Chord diagram
  # par(mfrow=c(1,1))
  # netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
  #netAnalysis_contribution(cellchat, signaling = pathways.show)
  # Access all the signaling pathways showing significant communications
  pathways.show.all <- cellchat@netP$pathways
  # Compute the network centrality scores
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
  # Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
  # netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
  # 
  
  #Visualize the dominant senders (sources) and receivers (targets) in a 2D space
  # Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
  pdf(paste("2D_space_",age,".pdf", sep = ""))
  gg1 <- netAnalysis_signalingRole_scatter(cellchat)
  # Signaling role analysis on the cell-cell communication networks of interest
  # gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("TGFb", "BMP"))
  print(gg1)
  dev.off()
  
  pdf(paste("pathway_heatmap_",age,".pdf", sep = ""), 10,8)
  # Identify signals contributing most to outgoing or incoming signaling of certain cell groups
  # Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
  ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
  ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
  print(ht1 + ht2)
  dev.off()
  saveRDS(cellchat, file = paste("cellchat_",age, ".rds",sep = ""))
  
  if (age == "09_11w"){
    df.net <- subsetCommunication(cellchat, sources.use = c(1,3,4), targets.use = c(1))
    write.csv(df.net,paste("df.net.OE_",age,".csv", sep = ""))
  }
  if (age == "12_13w"){
    df.net <- subsetCommunication(cellchat, sources.use = c(1,2,3,4), targets.use = c(1))
    write.csv(df.net,paste("df.net.DE_",age,".csv", sep = ""))
    
  }
  if (age == "14_16w"){
    df.net <- subsetCommunication(cellchat, sources.use = c(1,2,3,4,6), targets.use = c(2))
    write.csv(df.net,paste("df.net.OEE_",age,".csv", sep = ""))
    # df.net <- subsetCommunication(cellchat, sources.use = c(1,3,4,5,6,7), targets.use = c(5))
    # write.csv(df.net,paste("df.net.PA-1_",age,".csv", sep = ""))
    
  }
  if (age == "17_19w"){
    df.net <- subsetCommunication(cellchat, sources.use = c("SR-2","SI-2", "PA","Pre.O", "OD-1" ), targets.use = c("PA"))
    write.csv(df.net,paste("antagonist_df.net.PA_",age,".csv", sep = ""))
    # df.net <- subsetCommunication(cellchat, sources.use = c(1,2,3,5,6,7,8), targets.use = c(5))
    # write.csv(df.net,paste("df.net.PA-2_",age,".csv", sep = ""))
    
  }
  if (age == "20_22w"){
    df.net <- subsetCommunication(cellchat, sources.use = c(1,2,3,4), targets.use = c(3))
    write.csv(df.net,paste("df.net.Am_",age,".csv", sep = ""))
    df.net <- subsetCommunication(cellchat, sources.use = c(1,2,4), targets.use = c(2))
    write.csv(df.net,paste("df.net.PA_",age,".csv", sep = ""))
    
  }
  
  
}

age= "09_11w"
age= "12_13w"
age= "14_16w"
age= "17_19w"
age= "20_22w"

for (age in levels(colData(mixed_cds)$age_group)) {
  
  cellchat <- readRDS(paste("/home/ubuntu/data_sci_RNA_seq/",age,"/cellchat_",age,".rds", sep=""))
  setwd("~/data_sci_RNA_seq")
  data.dir <- age  # Create a directory to save figures
  #dir.create(data.dir)
  setwd(data.dir)
  
  if (age == "09_11w"){
  pdf(paste("bubble_oe_",age,".pdf", sep = ""),5,14)
  print(netVisual_bubble(cellchat, sources.use = c(1,3,4), targets.use = c(1), remove.isolate = FALSE, font.size =15))
  
  dev.off()
  
  groupSize <- as.numeric(table(cellchat@idents))
  pdf(paste("netCircle_count_",age,".pdf", sep = ""))
  print(netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions", margin = 0.1, vertex.label.cex = 3, , arrow.size = 0.8))
  dev.off()
  }
  
  if (age == "12_13w"){
    pdf(paste("bubble_de_",age,".pdf", sep = ""),5,14)
    print(netVisual_bubble(cellchat, sources.use = c(1,2,3,4), targets.use = c(1), remove.isolate = FALSE, font.size =15))
    
    dev.off()
    
    
    groupSize <- as.numeric(table(cellchat@idents))
    pdf(paste("netCircle_count_",age,".pdf", sep = ""))
    print(netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions", vertex.label.cex = 3, , arrow.size = 0.8))
    dev.off()
  }
  
  if (age == "14_16w"){
    pdf(paste("bubble_oee_",age,".pdf", sep = ""),5,14)
    print(netVisual_bubble(cellchat, sources.use = c(1,2,3,4,6), targets.use = c(2), remove.isolate = FALSE, font.size =15))
    
    dev.off()
    
    groupSize <- as.numeric(table(cellchat@idents))
    pdf(paste("netCircle_count_",age,".pdf", sep = ""))
    print(netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions", vertex.label.cex = 3, , arrow.size = 0.8))
    dev.off()
  }
  if (age == "17_19w"){
    pdf(paste("antagonist_bubble_pa1_",age,".pdf", sep = ""),5,14)
    print(netVisual_bubble(cellchat, sources.use = c("SR-2","SI-2", "PA","Pre.O", "OD-1" ), targets.use = c("PA"),  remove.isolate = FALSE, font.size =15))
    
    dev.off()
    
    groupSize <- as.numeric(table(cellchat@idents))
    pdf(paste("netCircle_count_",age,".pdf", sep = ""))
    print(netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions", vertex.label.cex = 3, , arrow.size = 0.8))
    dev.off()
  }
  if (age == "20_22w"){
    pdf(paste("bubble_PA_",age,".pdf", sep = ""),5,14)
    print(netVisual_bubble(cellchat, sources.use = c(1,2,4), targets.use = c(2),  remove.isolate = FALSE, font.size =15))
    
    dev.off()
    pdf(paste("bubble_Am_",age,".pdf", sep = ""),5,14)
    print(netVisual_bubble(cellchat, sources.use = c(1,2,3,4), targets.use = c(3),  remove.isolate = FALSE, font.size =15))
    
    dev.off()
    
    groupSize <- as.numeric(table(cellchat@idents))
    pdf(paste("netCircle_count_",age,".pdf", sep = ""))
    print(netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions", vertex.label.cex = 3, , arrow.size = 0.8))
    dev.off()
  }
  
  

}



###############
mixed_cds_9 <- mixed_cds[, colData(mixed_cds)$age_group =="09_11w"]

included <- names(table( colData(mixed_cds_9)$assigned_cell_type.simp)[table( colData(mixed_cds_9)$assigned_cell_type.simp) > 20])
mixed_cds_9 <- mixed_cds_9[, colData(mixed_cds_9)$assigned_cell_type.simp %in% included]


data.input = normalized_counts(mixed_cds_9, norm_method = "log") # normalized data matrix
row.names(data.input) = rowData(mixed_cds_9)$gene_short_name
identity = data.frame(group = droplevels(colData(mixed_cds_9)$assigned_cell_type.simp), row.names = names(colData(mixed_cds_9)$assigned_cell_type.simp)) # create a dataframe consisting of the cell labels
unique(identity$group) # check the cell labels




cellchat <- createCellChat(object = data.input, meta = identity,  group.by = "group")

# cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
# cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group


CellChatDB <- CellChatDB.human # use CellChatDB.human if running on human data
showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling","Cell-Cell Contact")) # use Secreted Signaling for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat@DB <- CellChatDB.use # set the used database in the object
# use all CellChatDB for cell-cell communication analysis


cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
parallel::detectCores()
cores_num <- parallel::detectCores()
cores = ifelse(cores_num>1,cores_num-1,1)
future::plan("multiprocess", workers = cores) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat) # if error, check cellchat@data.signaling
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)


cellchat <- computeCommunProb(cellchat, type = "triMean", trim = 0.1, population.size = TRUE ) #truncatedMean
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 20)

# {We provide a function subsetCommunication to easily access the inferred cell-cell communications of interest. For example,
#   
#   df.net <- subsetCommunication(cellchat) returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways
#   
#   df.net <- subsetCommunication(cellchat, sources.use = c(4), targets.use = c(2)) gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.
#   
#   df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb")) gives the inferred cell-cell communications mediated by signaling WNT and TGFb.
# }

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# To examine the signaling sent from each cell group
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
pdf("CellChat.pdf")
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

data.dir <- 'chat_results'  # Create a directory to save figures
dir.create(data.dir)
setwd(data.dir)


pathways.show <- c("TGFb") 
levels(cellchat@idents)
vertex.receiver = c(1,2) # a numeric vector
# Hierarchy plot
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, vertex.size = groupSize)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle", vertex.size = groupSize)
# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")


netAnalysis_contribution(cellchat, signaling = pathways.show)



# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = seq(1,2)
data.dir <- 'chat_results'  # Create a directory to save figures
dir.create(data.dir)
setwd(data.dir)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "circle")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}



# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use = seq(1,4), targets.use = c(1,2), remove.isolate = FALSE)
#> Comparing communications on a single object

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# show all the interactions sending from Inflam.FIB
pdf("CellChat_chord.pdf")
netVisual_chord_gene(cellchat, sources.use = seq(1,4), targets.use = c(1), lab.cex = 0.5,legend.pos.y = 30)
dev.off()


# Plot the signaling gene expression distribution using violin/dot plot
plotGeneExpression(cellchat, signaling = "TGFb")



# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)


#Visualize the dominant senders (sources) and receivers (targets) in a 2D space
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("TGFb", "BMP"))
gg1 + gg2



# Identify signals contributing most to outgoing or incoming signaling of certain cell groups
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

# Signaling role analysis on the cell-cell communication networks of interest
ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("TGFb", "BMP"))

# cellchat <- netAnalysis_signalingRole(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# 
# netVisual_signalingRole(cellchat, signaling = pathways.show)

nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat,  slot.name = "netP", pattern = "outgoing", k = nPatterns)
# river plot
netAnalysis_river(cellchat, pattern = "outgoing")
# dot plot
netAnalysis_dot(cellchat, pattern = "outgoing")


nPatterns = 1
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
# river plot
netAnalysis_river(cellchat, pattern = "incoming")
# dot plot
netAnalysis_dot(cellchat, pattern = "incoming")


saveRDS(cellchat, file = "cellchat_example.rds")




df.net <- subsetCommunication(cellchat, sources.use = c(4), targets.use = c(2))
write.csv(df.net,"df.net.mesen_de.csv")
