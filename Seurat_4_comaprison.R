#install.packages('Seurat')
library(Seurat)
library(monocle3)

OE_Am_cds <- readRDS("OE_Am_cds_may_27.rds")
Amelo_cds_MTi <-readRDS("epi_cds_MTi_70_manhattan_8L_0.3_7k_pesudo_regressed_assigned_cell_type_7.rds")


#simplify the labels
colData(Amelo_cds_MTi)$simplified <- dplyr::recode_factor(colData(Amelo_cds_MTi)$assigned_cell_type.4,
                                                          "OE" = "OE",
                                                          "DE-Prog" = "DE-Prog",
                                                          "EK" = "IK/EK",
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

# remove MT genes
List_of_MT_genes <-
  OE_Am_cds@rowRanges@elementMetadata@listData$gene_short_name[grep("^MT-",
                                                                    OE_Am_cds@rowRanges@elementMetadata@listData$gene_short_name)]
OE_Am_cds <-
  OE_Am_cds[!OE_Am_cds@rowRanges@elementMetadata@listData[["gene_short_name"]] %in% List_of_MT_genes,]


OE_Am_cds_sub <- OE_Am_cds_down[,colData(OE_Am_cds_down)$new_clusters %in% c(1:10)]
OE_Am_cds_sub <- OE_Am_cds[,colData(OE_Am_cds)$diff_time_point  %in% c("day_16")]

mat_Amelo <- counts(Amelo_cds_MTi)
rownames(mat_Amelo) <- rowData(Amelo_cds_MTi)$gene_short_name

mat_OE_Am <- counts(OE_Am_cds)
rownames(mat_OE_Am) <- rowData(OE_Am_cds)$gene_short_name

mat_OE_Am_sub <- counts(OE_Am_cds_sub)
rownames(mat_OE_Am_sub) <- rowData(OE_Am_cds_sub)$gene_short_name

mat_OE_Am_down <- counts(OE_Am_cds_down)
rownames(mat_OE_Am_down) <- rowData(OE_Am_cds_down)$gene_short_name


preAnnotation <-colData(Amelo_cds_MTi)$simplified
cell_id <- colnames(Amelo_cds_MTi)
meta_clusters <- data.frame(row.names=cell_id,preAnnotation)
meta_clusters$experiment<- "Fetal"

preAnnotation <-colData(OE_Am_cds)$clusters_named
cell_id <- colnames(OE_Am_cds)
meta_clusters2 <- data.frame(row.names=cell_id,preAnnotation)
meta_clusters2$experiment<- "InVitro"

preAnnotation <-colData(OE_Am_cds_down)$clusters_named
cell_id <- colnames(OE_Am_cds_down)
meta_clusters2 <- data.frame(row.names=cell_id,preAnnotation)
meta_clusters2$experiment<- "InVitro"

preAnnotation <-colData(OE_Am_cds_sub)$clusters_named
cell_id <- colnames(OE_Am_cds_sub)
meta_clusters3 <- data.frame(row.names=cell_id,preAnnotation)
meta_clusters3$experiment <- "InVitro"
meta_clusters3$time_point <- colData(OE_Am_cds_sub)$diff_time_point 


seurat_Amelo <-CreateSeuratObject(counts = mat_Amelo, project = "integration", min.cells = 1, min.features = 1, meta.data = meta_clusters)
seurat_OE_Am <-CreateSeuratObject(counts = mat_OE_Am_down, project = "integration", min.cells = 1, min.features = 1, meta.data = meta_clusters2)
seurat_OE_Am_sub <-CreateSeuratObject(counts = mat_OE_Am_sub, project = "integration", min.cells = 1, min.features = 1, meta.data = meta_clusters3)


integration.list <- list(seurat_Amelo,seurat_OE_Am)


### building reference
### I'll use inicors & molar and mixed as 3 sets

jaw_AMelo_cds <- Amelo_cds_MTi[,colData(Amelo_cds_MTi)$sample_id %in% c("Jaw_9w") ]
incisor_AMelo_cds <- Amelo_cds_MTi[,colData(Amelo_cds_MTi)$sample_id %in% c("incisors_12", "incisors_17","incisors_20") ]
molar_AMelo_cds <- Amelo_cds_MTi[,colData(Amelo_cds_MTi)$sample_id %in% c("molars_12", "molars_17","molars_20") ]
mixed_AMelo_cds <- Amelo_cds_MTi[,colData(Amelo_cds_MTi)$sample_id %in% c("teeth_14") ]


mat_jaw <- counts(jaw_AMelo_cds)
rownames(mat_jaw) <- rowData(jaw_AMelo_cds)$gene_short_name

preAnnotation <-colData(jaw_AMelo_cds)$simplified
cell_id <- colnames(jaw_AMelo_cds)
meta_jaw <- data.frame(row.names=cell_id,preAnnotation)
meta_jaw$sample <- "jaw"

mat_incisor <- counts(incisor_AMelo_cds)
rownames(mat_incisor) <- rowData(incisor_AMelo_cds)$gene_short_name

preAnnotation <-colData(incisor_AMelo_cds)$simplified
cell_id <- colnames(incisor_AMelo_cds)
meta_incisor <- data.frame(row.names=cell_id,preAnnotation)
meta_incisor$sample <- "incisor"

mat_molar <- counts(molar_AMelo_cds)
rownames(mat_molar) <- rowData(molar_AMelo_cds)$gene_short_name

preAnnotation <-colData(molar_AMelo_cds)$simplified
cell_id <- colnames(molar_AMelo_cds)
meta_molar <- data.frame(row.names=cell_id,preAnnotation)
meta_molar$sample <- "molar"

mat_mixed <- counts(mixed_AMelo_cds)
rownames(mat_mixed) <- rowData(mixed_AMelo_cds)$gene_short_name

preAnnotation <-colData(mixed_AMelo_cds)$simplified
cell_id <- colnames(mixed_AMelo_cds)
meta_mixed <- data.frame(row.names=cell_id,preAnnotation)
meta_mixed$sample <- "mixed"


seurat_jaw <-CreateSeuratObject(counts = mat_jaw, project = "reference", min.cells = 1, min.features = 1, meta.data = meta_jaw)
seurat_incisor <-CreateSeuratObject(counts = mat_incisor, project = "reference", min.cells = 1, min.features = 1, meta.data = meta_incisor)
seurat_molar <-CreateSeuratObject(counts = mat_molar, project = "reference", min.cells = 1, min.features = 1, meta.data = meta_molar)
seurat_mixed <-CreateSeuratObject(counts = mat_mixed, project = "reference", min.cells = 1, min.features = 1, meta.data = meta_mixed)

reference.list <- list(seurat_jaw,seurat_incisor,seurat_molar,seurat_mixed)


for (i in 1:length(reference.list)) {
  reference.list[[i]] <- NormalizeData(reference.list[[i]], verbose = FALSE)
  reference.list[[i]] <- Seurat::FindVariableFeatures(reference.list[[i]], selection.method = "vst", 
                                             verbose = FALSE, nfeatures = 20000 ) #nrow(seurat_incisor)
}

#reference.list <- reference.integrated <- IntegrateData(anchorset = reference.anchors, dims = 1:30).list[c("seurat_jaw", "seurat_incisor", "seurat_molar","seurat_mixed" )] no need to subset
reference.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:50, anchor.features = unique(rowData(Amelo_cds_MTi)$gene_short_name))#, anchor.features = unique(rownames(mat_incisor)), assay = rep("RNA",4))
reference.integrated <- IntegrateData(anchorset = reference.anchors,normalization.method = "LogNormalize", dims = 1:50, features.to.integrate = unique(rowData(Amelo_cds_MTi)$gene_short_name), features = unique(rowData(Amelo_cds_MTi)$gene_short_name))


library(ggplot2)
library(cowplot)
library(patchwork)
# switch to integrated assay. The variable features of this assay are automatically set during
# IntegrateData
DefaultAssay(reference.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
reference.integrated <- ScaleData(reference.integrated, verbose = FALSE)
reference.integrated <- RunPCA(reference.integrated, npcs = 50, verbose = FALSE)
reference.integrated <- RunUMAP(reference.integrated, reduction = "pca", dims = 1:50, verbose = FALSE)
p1 <- DimPlot(reference.integrated, reduction = "umap", group.by = "sample")+ scale_x_reverse() +scale_y_reverse() 
p2 <- DimPlot(reference.integrated, reduction = "umap", group.by = "preAnnotation", label = TRUE, repel = TRUE) + 
   scale_x_reverse() +scale_y_reverse() #NoLegend()+
p1 + p2 

saveRDS(reference.integrated,"reference.integrated_amelo_set_Seurat_close_to_monocle.rds")



 # projection
# diff data  seurat_OE_Am is the query 
OE_Am.anchors <- FindTransferAnchors(reference = reference.integrated, query = seurat_OE_Am, 
                                        dims = 1:10, reference.reduction = "pca")
predictions <- TransferData(anchorset = OE_Am.anchors, refdata = reference.integrated$preAnnotation, 
                            dims = 1:30)
seurat_OE_Am <- AddMetaData(seurat_OE_Am, metadata = predictions)


reference.integrated <- RunUMAP(reference.integrated, dims = 1:30, reduction = "pca", return.model = TRUE)  # already done
seurat_OE_Am <- MapQuery(anchorset = OE_Am.anchors, reference = reference.integrated, query = seurat_OE_Am, 
                           refdata = list(preAnnotation = "preAnnotation"), reference.reduction = "pca", reduction.model = "umap")


p1 <- DimPlot(reference.integrated, reduction = "umap", group.by = "preAnnotation", label = TRUE, label.size = 3, 
              repel = TRUE) + NoLegend() + ggtitle("DE Clusters in Seurat")+ scale_x_reverse() +scale_y_reverse() 
p2 <- DimPlot(seurat_OE_Am, reduction = "ref.umap", group.by = "predicted.preAnnotation", label = TRUE, 
              label.size = 3, repel = TRUE)  + ggtitle("Projection of inVtro Diff")+ scale_x_reverse() +scale_y_reverse() #+ NoLegend()
p1 + p2

p3 <- DimPlot(seurat_OE_Am, reduction = "ref.umap", group.by = "preAnnotation", label = TRUE, 
              label.size = 3, repel = TRUE)  + ggtitle("Projection of inVtro Diff")+ scale_x_reverse() +scale_y_reverse() #+ NoLegend()
p4 <-FeaturePlot(seurat_OE_Am,reduction = "ref.umap",  features = "AMBN", label = F, 
            label.size = 3, repel = TRUE, order= T,pt.size=1)  + ggtitle("AMBN Expression")+ scale_x_reverse() +scale_y_reverse() #+ NoLegend()

p1 + p2+ p3 + p4



#merge reference and query
reference.integrated$id <- 'reference'
seurat_OE_Am$id <- 'query'
refquery <- merge(reference.integrated, seurat_OE_Am)
refquery[["umap"]] <- merge(reference.integrated[["umap"]], seurat_OE_Am[["ref.umap"]])
refquery <- RunUMAP(refquery, reduction = 'umap', dims = 1:30)
DimPlot(refquery, group.by = 'id', shuffle = TRUE)


# To copy dataset data from mmonocle cds to seurat dataset
reference.import <- reference.integrated

Amelo_cds_MTi

gene_loading <- Amelo_cds_MTi@preprocess_aux@listData[["gene_loadings"]]

gene_id <-  rowData(Amelo_cds_MTi)[,1:2]
gene_id <- gene_id[row.names(gene_loading),]

row.names(gene_loading) <- gene_id$gene_short_name

gene_loading <- gene_loading[row.names(reference.import@reductions[["pca"]]@feature.loadings),]

reference.import@reductions[["pca"]]@feature.loadings <- gene_loading

reference.import@reductions[["pca"]]@cell.embeddings <- Amelo_cds_MTi@int_colData@listData[["reducedDims"]]@listData[["PCA"]][row.names(reference.import@reductions[["pca"]]@cell.embeddings),]

reference.import@reductions[["umap"]]@cell.embeddings <- Amelo_cds_MTi@int_colData@listData[["reducedDims"]]@listData[["UMAP"]][row.names(reference.import@reductions[["umap"]]@cell.embeddings),]

reference.import@reductions[["pca"]]@cell.embeddings <- Amelo_cds_MTi@int_colData@listData[["reducedDims"]]@listData[["Aligned"]][row.names(reference.import@reductions[["pca"]]@cell.embeddings),]

colnames(x = reference.import[["umap"]]@cell.embeddings) <- paste0("UMAP_", 1:2)



####################################

OE_Am.anchors <- FindTransferAnchors(reference = reference.import, query = seurat_OE_Am, 
                                     dims = 1:50, reference.reduction = "pca",scale = F)
predictions <- TransferData(anchorset = OE_Am.anchors, refdata = reference.integrated$preAnnotation, 
                            dims = 1:15)
seurat_OE_Am <- AddMetaData(seurat_OE_Am, metadata = predictions)


reference.integrated <- RunUMAP(reference.integrated, dims = 1:30, reduction = "pca", return.model = TRUE)  # already done
seurat_OE_Am <- MapQuery(anchorset = OE_Am.anchors, reference = reference.import, query = seurat_OE_Am, 
                         refdata = list(preAnnotation = "preAnnotation"), reference.reduction = "pca", reduction.model = "umap")


# Set it globally:
options(ggrepel.max.overlaps = Inf)

p1 <- DimPlot(reference.import, reduction = "umap", group.by = "preAnnotation", label = TRUE, label.size = 8, 
              repel = TRUE, pt.size= 1) + NoLegend() + ggtitle("DE Clusters in Seurat")#+ scale_x_reverse() +scale_y_reverse() 
p2 <- DimPlot(seurat_OE_Am, reduction = "ref.umap", group.by = "predicted.preAnnotation", label = TRUE, 
              label.size = 3, repel = TRUE)  + ggtitle("Projection of inVtro Diff")#+ scale_x_reverse() +scale_y_reverse() #+ NoLegend()
p1 + p2

p3 <- DimPlot(seurat_OE_Am, reduction = "ref.umap", group.by = "preAnnotation", label = TRUE, 
              label.size = 3, repel = TRUE)  + ggtitle("Projection of inVtro Diff")#+ scale_x_reverse() +scale_y_reverse()+ ggrepel::geom_text_repel(box.padding = 0.5, max.overlaps = Inf) #+ NoLegend()
p4 <-FeaturePlot(seurat_OE_Am,reduction = "ref.umap",  features = "AMBN", label = F, 
                 label.size = 3, repel = TRUE, order= T,pt.size=1)  + ggtitle("AMBN Expression")#+ scale_x_reverse() +scale_y_reverse() #+ NoLegend()

patchwork_obj <- p1 + p2#+ p3 + p4 

# function to make the patchwork   object has the same scale for each plot
p_ranges_x <- c(ggplot_build(patchwork_obj[[1]])$layout$panel_scales_x[[1]]$range$range,
                ggplot_build(patchwork_obj[[2]])$layout$panel_scales_x[[1]]$range$range)

p_ranges_y <- c(ggplot_build(patchwork_obj[[1]])$layout$panel_scales_y[[1]]$range$range,
                ggplot_build(patchwork_obj[[2]])$layout$panel_scales_y[[1]]$range$range)


patchwork_obj & 
  xlim(min(p_ranges_x), max(p_ranges_x)) & 
  ylim(min(p_ranges_y), max(p_ranges_y))



#merge reference and query
reference.import$id <- 'Fetal'
reference.import$projection <- 'fetal'
seurat_OE_Am$id <- 'In Vitro'
seurat_OE_Am$projection <- colData(OE_Am_cds[,colnames(seurat_OE_Am)])$diff_time_point 
refquery <- merge(reference.import, seurat_OE_Am)
refquery[["umap"]] <- merge(reference.import[["umap"]], seurat_OE_Am[["ref.umap"]])
#refquery <- RunUMAP(refquery, reduction = 'umap', dims = 1:30)
DimPlot(refquery, group.by = 'id', shuffle = TRUE, cols = c("grey","red"), pt.size= 1)
pd<- DimPlot(refquery, group.by = 'projection', shuffle = F, cols = c("blue","red","grey"), pt.size= 1,combine = F)


plot_seurat_density(refquery,groups_column = 'id' , nrow = NULL,
                                ncol = NULL, x_reverse = FALSE, y_reverse =FALSE)

saveRDS(refquery, "refquery_projection_combined.rds")


seurat_OE_Am <- readRDS("seurat_OE_Am_projected.rds")
reference.import<-readRDS("reference.import.rds")




seurat_OE_Am_sub <- readRDS("seurat_OE_Am_sub_projected_by_imported_PCA.rds")
#######################################33

{OE_Am.anchors <- FindTransferAnchors(reference = reference.import, query = seurat_OE_Am_sub , 
                                     dims = 1:10, reference.reduction = "pca", scale = F, k.filter = NA, k.anchor = 200, n.trees = 100, nn.method = "annoy", reduction = "pcaproject" )#,l2.norm = T,project.query = T, features = c("AMBN","SP6","ENAM","SHH"))


#reference.integrated <- RunUMAP(reference.integrated, dims = 1:30, reduction = "pca", return.model = TRUE)  # already done
seurat_OE_Am_sub  <- MapQuery(anchorset = OE_Am.anchors, reference = reference.import, query = seurat_OE_Am_sub , 
                         refdata = list(preAnnotation = "preAnnotation"), reference.reduction = "pca", reduction.model = "umap")


# Set it globally:
options(ggrepel.max.overlaps = Inf)

p1 <- DimPlot(reference.import, reduction = "umap", group.by = "preAnnotation", label = TRUE, label.size = 8, 
              repel = TRUE, pt.size= 1) + NoLegend() + ggtitle("DE Clusters in Seurat")#+ scale_x_reverse() +scale_y_reverse() 
p2 <- DimPlot(seurat_OE_Am_sub , reduction = "ref.umap", group.by = "predicted.preAnnotation", label = TRUE, 
              label.size = 3, repel = TRUE)  + ggtitle("Projection of inVtro Diff")#+ scale_x_reverse() +scale_y_reverse() #+ NoLegend()
p2 <- DimPlot(seurat_OE_Am_sub , reduction = "umap", group.by = "predicted.preAnnotation", label = TRUE, 
              label.size = 3, repel = TRUE)  + ggtitle("Projection of inVtro Diff")#+ scale_x_reverse() +scale_y_reverse() #+ NoLegend()

patchwork_obj <- p1 + p2

p3 <- DimPlot(seurat_OE_Am_sub , reduction = "umap", group.by = "preAnnotation", label = TRUE,
              label.size = 3, repel = TRUE)  + ggtitle("Projection of inVtro Diff")#+ scale_x_reverse() +scale_y_reverse()+ ggrepel::geom_text_repel(box.padding = 0.5, max.overlaps = Inf) #+ NoLegend()
p4 <-FeaturePlot(seurat_OE_Am_sub ,reduction = "umap",  features = c("AMBN"), label = F, slot = "data",
                 label.size = 3, repel = TRUE, order= T,pt.size=1, min.cutoff = NA, blend = F, blend.threshold = 0.5)  + ggtitle("AMBN Expression")#+ scale_x_reverse() +scale_y_reverse() #+ NoLegend()

patchwork_obj <- p1 + p2+ p3 + p4

# function to make the patchwork   object has the same scale for each plot
p_ranges_x <- c(ggplot_build(patchwork_obj[[1]])$layout$panel_scales_x[[1]]$range$range,
                ggplot_build(patchwork_obj[[2]])$layout$panel_scales_x[[1]]$range$range)

p_ranges_y <- c(ggplot_build(patchwork_obj[[1]])$layout$panel_scales_y[[1]]$range$range,
                ggplot_build(patchwork_obj[[2]])$layout$panel_scales_y[[1]]$range$range)


patchwork_obj & 
  xlim(min(p_ranges_x), max(p_ranges_x)) & 
  ylim(min(p_ranges_y), max(p_ranges_y))

}

#merge reference and query
reference.import$id <- 'Fetal'
reference.import$projection <- 'fetal'
seurat_OE_Am_sub$id <- 'In Vitro'
seurat_OE_Am_sub$projection <- colData(OE_Am_cds_sub[,colnames(seurat_OE_Am_sub )])$diff_time_point 
refquery <- merge(reference.import, seurat_OE_Am_sub )
refquery[["umap"]] <- merge(reference.import[["umap"]], seurat_OE_Am_sub[["ref.umap"]])
#refquery <- RunUMAP(refquery, reduction = 'umap', dims = 1:30)
DimPlot(refquery, group.by = 'id', shuffle = TRUE, cols = c("grey","red"), pt.size= 1)
pd<- DimPlot(refquery, group.by = 'projection', shuffle = F, cols = c("blue","red","grey"), pt.size= 1,combine = F)
pd
DimPlot(refquery, group.by = 'projection', shuffle = F, cols = c("red","grey"), pt.size= 1,combine = F)[[1]]+ coord_fixed()
pd
plot_seurat_density(refquery,groups_column = 'preAnnotation' , nrow = NULL,
                    ncol = 5, x_reverse = FALSE, y_reverse =FALSE, drop = c("SI","IK/EK","SR"), cluster_order = c("OE","DE-Prog","OEE","PA","AM", "diff_9","diff_7","diff_6","diff_10","diff_5","diff_8"))

saveRDS(refquery, "refquery_projection_PA_strgiht_day16only_pscaleT.rds")



refquery <- readRDS("refquery_projection_PA_strgiht_day16only_pscaleF.rds")

refquery <- readRDS("refquery_projection_combined.rds")

refquery <- subset(x = refquery, subset = preAnnotation == c("diff_11", "diff_12","diff_13","diff_14","diff_15"),invert = TRUE)

seurat_OE_Am_sub <- SplitObject(refquery,split.by = "id")$`In Vitro`
reference.import <- SplitObject(refquery,split.by = "id")$`Fetal`



refquery <- merge(reference.import, seurat_OE_Am_sub )
refquery[["umap"]] <- merge(reference.import[["umap"]], seurat_OE_Am_sub[["umap"]])

refquery@meta.data
seurat_OE_Am_sub <- SplitObject(refquery,split.by = "id")$`In Vitro`



# Run the standard workflow for visualization and clustering
seurat_OE_Am <- NormalizeData(seurat_OE_Am, verbose = FALSE)
seurat_OE_Am <- FindVariableFeatures(seurat_OE_Am, selection.method = "vst", nfeatures = 10000)
seurat_OE_Am <- ScaleData(seurat_OE_Am, verbose = FALSE)
seurat_OE_Am <- RunPCA(seurat_OE_Am, npcs = 40, verbose = FALSE,features= reference.import@assays[["integrated"]]@var.features)
seurat_OE_Am <- RunUMAP(seurat_OE_Am, reduction = "pca", dims = 1:40, verbose = FALSE,n.neighbors = 100L, min.dist= 0.1, spread = 5)
seurat_OE_Am <- FindNeighbors(object = seurat_OE_Am)
seurat_OE_Am <- FindClusters(object = seurat_OE_Am)
p1 <- DimPlot(seurat_OE_Am, reduction = "umap", group.by = "ident")+ scale_x_reverse() +scale_y_reverse() 
p2 <- DimPlot(seurat_OE_Am, reduction = "umap", group.by = "preAnnotation", label = TRUE, repel = TRUE) + 
  scale_x_reverse() +scale_y_reverse() #NoLegend()+
p2
p1 + p2 

FeaturePlot(seurat_OE_Am ,reduction = "umap",  features = c("AMBN","SP6", "PITX2","SOX2"), label = F, slot = "data",
            label.size = 3, repel = TRUE, order= T,pt.size=1, min.cutoff = NA, blend = F, blend.threshold = 0.5) # + ggtitle("AMBN Expression")#+ scale_x_reverse() +scale_y_reverse() #+ NoLegend()




reference_diff <- readRDS("reference.integrated_Oe_AM_new.rds")

p1 <- DimPlot(reference_diff, reduction = "umap", group.by = "sample")+ scale_x_reverse() +scale_y_reverse() 
p2 <- DimPlot(reference_diff, reduction = "umap", group.by = "preAnnotation", label = TRUE, repel = TRUE) + 
  scale_x_reverse() +scale_y_reverse() #NoLegend()+

p1 + p2 
FeaturePlot(reference_diff ,reduction = "umap",  features = c("AMBN","SP6", "PITX2","SOX2"), label = F, slot = "data",
            label.size = 3, repel = TRUE, order= T,pt.size=1, min.cutoff = NA, blend = F, blend.threshold = 0.5) # + ggtitle("AMBN Expression")#+ scale_x_reverse() +scale_y_reverse() #+ NoLegend()






#Function to assign closer cluster
#first finding the mid point for each base clsuter  preAnnotation 

cluster_centers <- data.frame(row.names = unique(reference.import$full_annotation ), UMAP_1 = rep(0,length(unique(reference.import$full_annotation ))),  UMAP_2 = rep(0,length(unique(reference.import$full_annotation ))))
for (clu in  unique(reference.import$full_annotation )){
  temp_obj <- subset(x = reference.import, subset = full_annotation  == clu)
  
  cluster_centers[clu,1] <- mean(temp_obj@reductions$umap@cell.embeddings[,1])
  cluster_centers[clu,2]<- mean(temp_obj@reductions$umap@cell.embeddings[,2])
  
  
}


closest_cluster <- data.frame(row.names = row.names(seurat_OE_Am_sub@reductions$umap@cell.embeddings), nearest = rep(0,nrow(seurat_OE_Am_sub@reductions$umap@cell.embeddings)))
for (cell in  row.names(closest_cluster)){

  
  temp_distances <- DataFrame(row.names = row.names(cluster_centers), Distances = rep(0,nrow(cluster_centers)))
  for (clu in  row.names(cluster_centers)){
    temp_cord <- DataFrame(UMAP_1 = numeric(0), UMAP_2 = numeric(0))
    temp_cord <- rbind(temp_cord, cluster_centers[clu,])
    temp_cord <- rbind(temp_cord, as.data.frame(t(seurat_OE_Am_sub@reductions$umap@cell.embeddings[cell,])))
    temp_distances[clu,1] <-  dist(temp_cord)
  }
  closest <- row.names(temp_distances)[which.min(temp_distances$Distances)]
  closest_cluster[cell,1] <- closest

  
  
}

seurat_OE_Am_sub$new_projection <- closest_cluster
seurat_OE_Am_sub$new_projection <- paste0(seurat_OE_Am_sub$new_projection,"-like")
refquery <- merge(reference.import, seurat_OE_Am_sub )
refquery[["umap"]] <- merge(reference.import[["umap"]], seurat_OE_Am_sub[["umap"]])
p1 <- DimPlot(refquery, group.by = 'new_projection', shuffle = F, pt.size= 1,combine = F,  na.value = "grey")[[1]]+ NoLegend() + coord_fixed()
p2 <- DimPlot(seurat_OE_Am_sub, group.by = 'new_projection', shuffle = F, pt.size= 1,combine = F,  na.value = "grey")[[1]]+ coord_fixed()



patchwork_obj <- p1 + p2
# function to make the patchwork   object has the same scale for each plot
p_ranges_x <- c(ggplot_build(patchwork_obj[[1]])$layout$panel_scales_x[[1]]$range$range,
                ggplot_build(patchwork_obj[[2]])$layout$panel_scales_x[[1]]$range$range)

p_ranges_y <- c(ggplot_build(patchwork_obj[[1]])$layout$panel_scales_y[[1]]$range$range,
                ggplot_build(patchwork_obj[[2]])$layout$panel_scales_y[[1]]$range$range)

patchwork_obj & 
  xlim(min(p_ranges_x), max(p_ranges_x)) & 
  ylim(min(p_ranges_y), max(p_ranges_y))


count(seurat_OE_Am_sub@meta.data[,c("preAnnotation","new_projection")])

data_summary <- seurat_OE_Am_sub@meta.data[,c("preAnnotation","new_projection")] %>% group_by(preAnnotation,new_projection ) %>% count()

saveRDS(data_summary, "assigned_table_for_orginal_projection_all_diff_data_Full_annot.rds")


refquery_day_16 <- subset(x = refquery, subset = projection == c("day_16", NA))

p2 <- DimPlot(refquery_day_16, group.by = 'new_projection', shuffle = F, pt.size= 1,combine = F,  na.value = "grey")[[1]]+ coord_fixed()



reference.import$full_annotation <- colData(Amelo_cds_MTi)$assigned_cell_type.4




                                                    

