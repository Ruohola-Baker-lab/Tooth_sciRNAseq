
library(tidyr)
library(monocle3)
library(ggplot2)
library(dplyr)
library(Matrix.utils)


# library(devtools)

# cds <- Amelo_cds_MTi
# 
# clustering_column = "assigned_cell_type.4"
# age_clustering_column = "age_group"
# selected_markers = "auto"
# number_of_genes2goterm = 50
# number_of_markers = 3; number_of_genes2goterm = 50; gene_font_size= 10; keywords_font_size= 10; number_of_keywords = 7; keywords_blacklist= c()
 
custom_Heatmap <- function(cds, clustering_column, age_clustering_column, selected_markers, number_of_markers = 3, number_of_genes2goterm = 50, gene_font_size= 10, keywords_font_size= 10, number_of_keywords = 7,keywords_blacklist= c(),load_previous = T, auto_sort_columns = T, normalize_by_cluster = F, disable_legends =F){

  start_time <- Sys.time()
  

  
  if (packageVersion("dplyr") != "0.8.2" | packageVersion("dbplyr") != "1.3.0"  ) {
    # #detach(dplyr, pos = 2L, unload = FALSE, character.only = FALSE,
    #        force = T)
    # #detach(dbplyr, pos = 2L, unload = FALSE, character.only = FALSE,
    #        force = T)

    remove.packages(c("dplyr","dbplyr"))

    remotes::install_github("hadley/dplyr@v0.8.2", dependencies = FALSE , upgrade = "never")
    #remotes::install_github("r-lib/gert", dependencies = FALSE, upgrade = "never")
    # #install.packages("devtools")
    # #devtools::install_github("hadley/devtools", dependencies = FALSE)
    devtools::install_url("https://cran.r-project.org/src/contrib/Archive/dbplyr/dbplyr_1.3.0.tar.gz",  dependencies = FALSE,upgrade = "never")
    #BiocManager::install("AnnotationHub", version = "3.10", dependencies = TRUE, ask=FALSE)
  }

assertthat::assert_that(!is.null(cds), msg = "cds object doesn't exist. Check the spelling!")
assertthat::validate_that(!is.null(clustering_column),msg = "please add clustering_column name. Otherwise, the default clustering will be loaded")
if (clustering_column == "cluster"){
  colData(cds)$clustering <- clusters(cds)
  clustering_column<- "clustering"
  
}
assertthat::assert_that(!is.null(age_clustering_column),msg = "please add age_clustering_column name")
assertthat::validate_that(!is.null(colData(cds)[[clustering_column]]), msg = "clustering_column doesn't exist in cds. Check the spelling!")
assertthat::assert_that(!is.null(colData(cds)[[age_clustering_column]]), msg = "age_clustering_column doesn't exist in cds. Check the spelling!")

cds_name <-   deparse(match.call()$cds) # deparse get the charcter nameof an object, while gsub here remove . in the name #gsub("\\.", "_", .deparse(substitute(cds)) doesn't work in function
data.dir <- cds_name  # Create a directory to save data
dir.create("custom_Heatmap_output", showWarnings = F)
data_file <- paste("custom_Heatmap_output","/",cds_name, "_heatmapData.Rdata", sep="")


  if (load_previous == T & file.exists(data_file)){
    load(data_file)
    loaded <- T
  
}else {
  loaded <- F
}
  




# if (!dir.exists(data.dir)){
#   
#   loaded = F
#   #setwd(data.dir)
# }




if (loaded == F) {
  

# remove rowsum == 0 (aka genes with zero expression in all cells)
print("removing genes with zero expression in all cells")
cds <- cds[rowSums(counts(cds)) !=0 ,]
colData(cds)[[age_clustering_column]] <- factor(colData(cds)[[age_clustering_column]], levels = gtools::mixedsort(unique(colData(cds)[[age_clustering_column]])))
levels(colData(cds)[[age_clustering_column]])[levels(colData(cds)[[age_clustering_column]])=="9_11w"] <- "09_11w"

# to make sure to drop unused levels
colData(cds)[[clustering_column]] <- droplevels(colData(cds)[[clustering_column]])
colData(cds)[[age_clustering_column]] <- droplevels(colData(cds)[[age_clustering_column]])

# extract expression matrix per group
print("extracting genes expression matrix per cluster")
cell_group_df <- data.frame(cell=colnames(cds) , group = colData(cds)[[clustering_column]] ) #clusters(cds)
mat_all <- aggregate_gene_expression(cds,cell_group_df = cell_group_df, norm_method ="size_only" )
#colnames(mat_all) <- seq(1,ncol(mat_all))
if (normalize_by_cluster){
  mat_all <- as.data.frame(as.matrix(mat_all))/table(colData(cds)[[clustering_column]])
}



# extract expression matrix per age group
print("extracting genes expression matrix per age_group")
cell_group_df <- data.frame(cell=colnames(cds) , group = colData(cds)[[age_clustering_column]] ) #clusters(cds)
mat_age <- aggregate_gene_expression(cds,cell_group_df = cell_group_df, norm_method ="size_only" )
mat_age <- mat_age[,levels(colData(cds)[[age_clustering_column]])] #c("09_11w","12_13w","14_16w","17_19w","20_22w")

# reorder the columns if necessary
# mat_all<- mat_all[,c("OE","DE-Prog","OEE","SI-1","SI-2","SI-3","EK","SR-1","SR-2","PA-1","PA-2","AM-1","AM-2")]

# extract a matrix of percentage of cells expressing each gene per cluster
print("calculating the percentage of cells expressing each gene per cluster")
percentage_mat <-Matrix::Matrix(0, nrow = nrow(cds), ncol  = ncol(mat_all) )
colnames(percentage_mat) <- colnames(mat_all)
row.names(percentage_mat) <- row.names(mat_all)

for (i in colnames(mat_all)){
  non_zero<- Matrix::rowSums(counts(cds[, colData(cds)[[clustering_column]] == i])!=0)
  all_sum <- Matrix::rowSums(counts(cds[, colData(cds)[[clustering_column]] == i])!=0 | counts(cds[, colData(cds)[[clustering_column]] == i]) ==0)
  percentage_mat[,which(i == colnames(mat_all))] <-   (non_zero/all_sum)*100
}



# Scale by row
print("scaling the expression matrix per row")
mat_all <- t(scale(t(mat_all)))  #scale = apply(t(mat_all), 2, sd, na.rm = TRUE)
mat_all[is.nan(mat_all)] <- 0

mat_age <- t(scale(t(mat_age)))  #scale = apply(t(mat_age), 2, sd, na.rm = TRUE)
mat_age[is.nan(mat_age)] <- 0

# determine max
print("selectig the max expression per gene")
max_rows <-apply(X=mat_all, MARGIN=1, FUN=max)
#mat_all[mat_all ==max_rows] <-  0



# select the max genes and order them diagonal 
print("ordering the max expression genes diagonally")
df <- c()

for (c in seq(ncol(mat_all))){
  df <- c(df,sort(mat_all[mat_all[,c] >= max_rows, c], decreasing = T))
}


# to keep track of which cluster has the heighest expression per particular gene
print("matching max genes to their clusters")
highest <-c()

for (c in seq(ncol(mat_all))){
  highest <- c(highest,rep(colnames(mat_all)[c], length(mat_all[mat_all[,c] >= max_rows, c])))
}


row.names(mat_all) = rowData(cds[row.names(mat_all),])$gene_short_name

percentage_mat <-  percentage_mat[names(df),]

mat_age <- mat_age[names(df),]

ensemble_df <-as.data.frame(names(df))
names(df) = rowData(cds[names(df),])$gene_short_name
row.names(mat_age) = rowData(cds[row.names(mat_age),])$gene_short_name




# create dataframe of the ordered genes
print("building the annotation dataframe")
df <- data.frame(id= names(df), var = df)

# label rows by serial numbers for easier slicing
row.names(df) <- 1:nrow(df)
names(highest) <- 1:nrow(df)
row.names(percentage_mat) <- 1:nrow(df)
row.names(mat_age) <- 1:nrow(df)


print("remove duplicate genes while keeping max")
## using unique  to exclude duplicate maxima
df <- unique(subset(df, var==ave(var, id, FUN=max)) )


# slice out duplicate from highest list 
highest <- highest[row.names(df)]
ensemble_df <- ensemble_df[row.names(df),]

# slice out duplicate from percentage_mat 
percentage_mat <- percentage_mat[row.names(df),]
row.names(percentage_mat) = df$id

mat_age <- mat_age[row.names(df),]
row.names(mat_age) = df$id

# add it to the dataframe
df <- cbind(df,highest )
df <- cbind(df,ensemble_df )

# rename rownames back
row.names(df) <- df$id
#mat_all <- mat_all[df$id,]

cumulative_index <- cumsum(table(df$highest)[colnames(mat_all)])

# matching cell percentage to highest expressing celltype per gene
print("matching cell percentage to highest expressing celltype per geney")
cell_percentage <- c(percentage_mat[1:cumulative_index[1],1])
for (c in seq(2,ncol(mat_all))){
  cell_percentage <- c(cell_percentage,percentage_mat[(cumulative_index[c-1]+1):cumulative_index[c],c])
}


# add it to the dataframe
df <- cbind(df, cell_percentage)

#df$major_age <- factor(colnames(mat_age)[apply(mat_age,1,which.max)],levels = levels(colData(cds)[[age_clustering_column]]))





# important
df_sub <- as.data.frame(colnames(mat_all))
colnames(df_sub) <- "Clusters"
row.names(df_sub) <- df_sub$Clusters

# cds <- Amelo_cds_MTi
# clustering_column <- "assigned_cell_type.4"
# age_clustering_column <- "age_group"
# add annotation for age group
print("calculating the age_score and the age_group annotation")
age_group_df <- table(colData(cds)[,c(age_clustering_column,clustering_column)])
age_group_sum <- rowSums(age_group_df)
age_group_df <- (age_group_df/age_group_sum)*100
df_sub$major_age <-row.names(age_group_df)[apply(age_group_df,2,which.max)]
df_sub$major_age <- factor(df_sub$major_age, levels = levels(colData(cds)[[age_clustering_column]]))



age_group_df <- t(t(age_group_df)/colSums(age_group_df))
df_sub$age_score <- (scales::rescale(colSums((age_group_df) * c(1:nrow(age_group_df))), to = c(0,100)))

# auto order
if (auto_sort_columns == T){
age_group_df <- age_group_df[,order(df_sub$age_score) ]
mat_all <- mat_all[, order(df_sub$age_score) ]
df_sub <- df_sub[ order(df_sub$age_score),]
}

mat_all <- mat_all[row.names(df),]


#library(ComplexHeatmap)

# print("plotting heatmap to pdf")
# pdf("custom_Heatmap_TESTING.pdf")
# ch <- ComplexHeatmap::Heatmap(
#   name = "expression",
#   mat_all[row.names(df), ] ,
#   cluster_columns = F,
#   cluster_rows = F,
#   show_row_name = F,
#   show_column_names = T,
#   column_names_side = "top",
#   column_names_rot = 45,
#   bottom_annotation  = HeatmapAnnotation(age_score =df_sub[, 3] , age_group=t(age_group_df), annotation_name_side="right"), #major_age= df_sub[, 2],
#   left_annotation  = rowAnnotation(foo =anno_mark(at = which(row.names(mat_all) %in% selected_markers ), labels =row.names(df)[row.names(df) %in% selected_markers ], side = "left", labels_gp = gpar(fontsize  = 5))),
#   row_split  = factor(df$highest, levels = colnames(mat_all)),
#   row_title = " ",
#   use_raster = T,
#   show_heatmap_legend = T
# )
# draw(ch)
# 
# dev.off()

 

library(ViSEAGO)

# load genes background
library(org.Hs.eg.db)
background=keys(org.Hs.eg.db, keytype ='ENSEMBL')
background=intersect(background,row.names(cds))
background <- background[background %in% row.names(cds)]


# connect to Ensembl . #biomart 2.42.1
# print("loading GO annotation database")
# Ensembl<-ViSEAGO::Ensembl2GO(biomart = "ensembl", host = "www.ensembl.org",
#                              version = NULL)
# # Display table of available organisms with Ensembl
# #ViSEAGO::available_organisms(Ensembl)
# # load GO annotations from Ensembl
# myGENE2GO<-ViSEAGO::annotate(
#   "hsapiens_gene_ensembl",
#   Ensembl
# )

# load gene selection
myGENE2GO <- readRDS("/home/ubuntu/data_sci_RNA_seq/myGENE2GO.rds")
  
go_df <- mat_all[1:200,]
row.names(go_df) <- 1:200

go_term_df <- mat_all[1:200,]
row.names(go_term_df) <- 1:200

gene_info <- read.delim("/home/ubuntu/data_sci_RNA_seq/Homo_sapiens.gene_info")
gene_info <- gene_info[gene_info$type_of_gene == "protein-coding",] # remove pesudogenes
gene_info <- gene_info[stringr::str_detect(gene_info$description, "ribosomal", negate = T), ]  #remove ribosomal genes
df_coding <- df[df$id %in% gene_info$Symbol,]
print("Finding top marker genes")
top <- monocle3::top_markers(cds[rowData(cds)$id %in% df_coding$ensemble_df,],group_cells_by= clustering_column , genes_to_test_per_group = number_of_genes2goterm +50)
top <- top[order(top$marker_score, decreasing = T),]
top <- top[!duplicated(top$gene_short_name),]
row.names(top) <- top$gene_short_name

df_coding$marker_score <- 0
df_coding[row.names(top),]$marker_score <- top$marker_score

# m = "AM-1"

print(paste("Calculate the top 200  goterm enrichment from the top ",number_of_genes2goterm," genes per cluster",sep = ""))
for (m in colnames(mat_all)){
  selection<- df_coding[df_coding$highest== m,]
  extra_genes <- selected_markers[selected_markers %in% row.names(selection)]
  extra_genes <- as.character(df_coding[extra_genes,]$ensemble_df) 

  selection <- selection %>% arrange(desc(marker_score),desc(cell_percentage)) %>%
    mutate(rank=cumsum(!duplicated(.))) %>% filter(rank<=max(number_of_genes2goterm,10))
  #selection<- dplyr::top_n(selection,max(number_of_genes2goterm,10),marker_score)
  
  selection<- background[background %in% as.character(selection$ensemble_df)]
  
  selection <- c(selection , extra_genes)
  
  BP <- ViSEAGO::create_topGOdata(
    geneSel = selection,
    allGenes = background,
    gene2GO = myGENE2GO,
    ont = "BP",
    nodeSize = 5
  )
  # perform TopGO test using clasic algorithm
  classic <- topGO::runTest(BP,
                            algorithm = "classic",
                            statistic = "fisher")
  # merge results from topGO
  BP_sResults <- ViSEAGO::merge_enrich_terms(Input = list(condition = c("BP", "classic")), envir = environment())
  go_df[,m] <-  as.data.frame(BP_sResults@data[order(condition.pvalue),])[,1][1:200]
  go_term_df[,m] <-  as.data.frame(BP_sResults@data[order(condition.pvalue),])[,2][1:200]
  write.csv(as.data.frame(BP_sResults@data), file = paste("custom_Heatmap_output","/",stringr::str_replace(m,"/","_") ,"_goterm_results.csv", sep=""))
  
}











print("save meta variables for quick rerun")
save(list= c("background","selection", "df","df_sub" ,"mat_all","myGENE2GO","BP","go_df","go_term_df","age_group_df", "df_coding"), file = paste("custom_Heatmap_output","/",cds_name, "_heatmapData.Rdata", sep=""))
print("write annotations into csv files for reference")
write.csv(df, file = paste("custom_Heatmap_output","/",cds_name, "_gene_Annoation.csv", sep=""))
write.csv(df_sub,file = paste("custom_Heatmap_output","/",cds_name, "_cluste_Annotation.csv", sep=""))
write.csv(go_term_df, file = paste("custom_Heatmap_output","/",cds_name, "_top_200_goterm.csv", sep=""))
write.csv(as.data.frame(BP_sResults@data), file = paste("custom_Heatmap_output","/",cds_name, "_goterm_results.csv", sep=""))
saveRDS(selected_markers, file = paste("custom_Heatmap_output","/",cds_name, "_selected_markers.rds", sep=""))
saveRDS(df_coding, file = paste("custom_Heatmap_output","/",cds_name, "_filtered_gene_annot.rds", sep=""))
}
# end of load skip
# 


#devtools::install_github("jokergoo/simplifyEnrichment", upgrade = "never")
library(simplifyEnrichment)
library(stringr)

print("exctract keywords from goterms")
list_word_df <- list()
for (s in  colnames(mat_all)){
  words <-count_word_from_GO(
    go_df[, which(colnames(mat_all) == s)],
    exclude_words = c(
      "system",
      "biosynthetic",
      "rna",
      "muscle",
      "neuron",
      "transcription",
      "compound",
      "tissue",
      "cell",
      "process",
      "multicellular",
      "nucleic",
      "tube",
      "macromolecule",
      "lung",
      "cardiac",
      "anatomical",
      "organismal",
      "nervous",
      "dnatemplated",
      "limb",
      "polymerase",
      "acidtemplated",
      "nucleobasecontaining",
      "nitrogen",
      "skeletal",
      "animal",
      "cellular",
      "pathway",
      "involved",
      "component",
      "receptor",
      "expression",
      "gene",
      "biological",
      "action",
      "potential",
      "immune",
      "neutrophil",
      "primirna",
      "contraction",
      "leukocyte",
      "myeloid",
      "histone",
      "viral",
      "ventricular",
      "myocyte",
      "organisms",
      "immunity",
      "enzyme",
      "part",
      "heart",
      "brain",
      "circulatory",
      "central",
      "supramolecule",
      "supramolecular",
      "organ",
      "linked",
      "cytosolic",
      "development",
      "negative",
      "positive",
      "organism",
      "diadenosine",
      "cotranslation",
      "cotranslational",
      "targeting",
      "endoderm",
      "translation",
      "regulation",
      "developmental",
      "localization",
      "stem",
      "symbiotic",
      "heterocycle",
      "aromatic",
      "nuclear",
      "translational",
      "establishment",
      "neurogenesis",
      "catabolic",
      "organic",
      "cyclic",
      "mediated",
      "organonitrogen",
      "decay",
      "interspecies",
      keywords_blacklist
      
      
    )
  )
  #word <- words$word[charmatch(words$word, go_term_df[,which(colnames(mat_all) == s)])]
  ranked_go <- str_extract(go_term_df[,which(colnames(mat_all) == s)], paste(words$word, collapse = "|"))
  ranked_go <- unique(ranked_go[!is.na(ranked_go)])
  word <- ranked_go[ranked_go != "NA"]
  
  # word <-ranked_go[,which(colnames(mat_all) == s)][charmatch(go_term_df[,which(colnames(mat_all) == s)],words$word)]
  # word <- word[!is.na(word)]
  # word <- word[word != "NA"]
  
  freq <- seq(length(word),1)
  words <- data.frame(word = word, freq = freq)
  list_word_df[[s]] <- words
}





# a function to handle graph unit division 
`/.unit` <- function (x, y) {
  x <- convertUnit(x, "pt")
  x <- as.numeric(x)
  unit(x / y, "pt")
}

# load gene selection



if ("auto"  %in% selected_markers) {
  print(paste("automatic selection of top ",number_of_markers," marker genes per group",sep = ""))
  selected_markers <- selected_markers[!selected_markers %in% "auto"  ]
  for (m in colnames(mat_all)){
    selection<- df_coding[df_coding$highest== m,]
    selection <- selection %>% arrange(desc(marker_score),desc(cell_percentage)) %>%
      mutate(rank=cumsum(!duplicated(.))) %>% filter(rank<=max(number_of_markers,1))
    #selection<- dplyr::top_n(selection,max(number_of_markers,1),marker_score)
    selected_markers <- c(selected_markers,as.character(selection$id))
    saveRDS(selected_markers, file = paste("custom_Heatmap_output","/",cds_name, "_selected_markers.rds", sep=""))
    # selection<- top[top$cell_group== m,]
    # selection<- top_n(selection,max(number_of_markers,1),marker_score)
    # selected_markers <- c(selected_markers,as.character(selection$gene_short_name))
    # 
  } 
  
}
#mat_all <- mat_all[row.names(df),]

df_sub$age_score <- sort(df_sub$age_score)
# rewrite colnames
#colnames(mat_all) <- levels(colData(cds)[[clustering_column]])[colnames(mat_all)]

print("construct the layout of keywords annotations")
splitting_factor <- factor(df$highest, levels = colnames(mat_all))
align_to = split(seq_len(nrow(mat_all)),splitting_factor )
align_to = align_to[names(align_to) != "0"]
align_to = align_to[names(align_to) %in% names(list_word_df)]

fontsize_range = c(1, keywords_font_size)
gbl = lapply(names(align_to), function(nm) {
  kw = list_word_df[[nm]][, 1][1:number_of_keywords]
  freq = list_word_df[[nm]][, 2][1:number_of_keywords]
  fontsize = scale_fontsize(freq, rg = c(1, max(10, freq)), fs = fontsize_range)
  
  word_cloud_grob(text = kw, fontsize = fontsize, max_width = unit(105, "mm"))
})
names(gbl) = names(align_to)

margin = unit(6, "pt")
gbl_h = lapply(gbl, function(x) convertHeight(grobHeight(x), "cm") + margin)
gbl_h = do.call(unit.c, gbl_h)

gbl_w = lapply(gbl, function(x) convertWidth(grobWidth(x), "cm"))
gbl_w = do.call(unit.c, gbl_w)
gbl_w = max(gbl_w) + margin


panel_fun = function(index, nm) {
  # background
  grid.rect(gp = gpar(fill = "#DDDDDD", col = NA))
  # border
  grid.lines(c(0, 1, 1, 0), c(0, 0, 1, 1), gp = gpar(col = "#AAAAAA"), 
             default.units = "npc")
  gb = gbl[[nm]]
  # a viewport within the margins
  pushViewport(viewport(x = margin/2, y = margin/2, 
                        width = grobWidth(gb), height = grobHeight(gb),
                        just = c("left", "bottom")))
  grid.draw(gb)
  popViewport()
}
library(ComplexHeatmap)
library(circlize)
print("plotting heatmap to pdf")
pdf(paste("custom_Heatmap_output","/",cds_name,"_custom_Heatmap.pdf", sep = ""), 10,7)
ht <- ComplexHeatmap::Heatmap(
  name = "expression",
  mat_all ,
  cluster_columns = F,
  cluster_rows = F,
  show_row_name = F,
  show_column_names = T,
  column_names_side = "top",
  column_names_rot = 45,
  bottom_annotation  = HeatmapAnnotation(age_score =df_sub[, 3] , age_group=t(age_group_df), annotation_name_side="right", show_legend = F, col = list(age_score =  colorRamp2(c(0,100 ), c("#FFEFDA", "orange")),age_group= colorRamp2(c(0,1 ), c("white", "darkgreen")))), #major_age= df_sub[, 2],
  left_annotation  = rowAnnotation(foo =anno_mark(at = which(row.names(mat_all) %in% selected_markers ), labels =row.names(df)[row.names(df) %in% selected_markers ], side = "left", labels_gp = gpar(fontsize  = gene_font_size))),
  row_split  = splitting_factor,
  row_title = " ",
  use_raster = T,
  show_heatmap_legend = F
) + rowAnnotation(keywords = anno_link(align_to = align_to, 
                                       which = "row", panel_fun = panel_fun, 
                                       size = gbl_h, gap = unit(2, "mm"), 
                                       width = gbl_w + unit(5, "mm"), # 5mm for the link
                                       link_gp = gpar(fill = "#DDDDDD", col = "#AAAAAA"), 
                                       internal_line = FALSE), show_legend = F)# you can set it to TRUE to see what happens
draw(ht, ht_gap = unit(2, "pt"))  

# col_fun = colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))
# lgd = Legend(col_fun = col_fun, title = "expression", direction = "vertical")
# draw(lgd, x = unit(0.995, "npc"), y = unit(0.88, "npc"), just = c("right", "top"))

if (disable_legends == F){
col_fun = colorRamp2(c(0,1 ), c("white", "darkgreen"))
lgd = Legend(col_fun = col_fun, title = "age_group",at = c(0, 0.25,0.50, 0.75,  1.0),labels =c("0", "", "50", "", "100") , direction = "vertical")
draw(lgd, x = unit(0.9, "npc"), y = unit(0.177, "npc"), just = c("left", "top"))

col_fun = colorRamp2(c(0,100 ), c("#FFEFDA", "orange"))
lgd = Legend(col_fun = col_fun, title = "age_score",at = c(0, 25,50,75, 100), labels =c("0", "", "50", "", "100"), direction = "vertical")
draw(lgd, x = unit(0.8, "npc"), y = unit(0.177, "npc"), just = c("left", "top"))

col_fun = colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))
lgd = Legend(col_fun = col_fun, title = "expression", direction = "vertical")
draw(lgd, x = unit(0.7, "npc"), y = unit(0.177, "npc"), just = c("left", "top"))
}

dev.off()
print("Done! results are saved in custom_Heatmap_output folder")
end_time <- Sys.time()

end_time - start_time

}

