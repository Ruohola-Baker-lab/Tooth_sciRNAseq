top_pathway <- function(scMlnet_results, deg, receiver_cell = receiver_cell,  TF_fun = "sum", main_fun = "sum", path_fun = "sum", lr_glom_normal =lr_glom_normal, show.sub.pathway = F, method = "KL_rec",proximity_order= NULL, palette, mirror.y = F, output_dir){
  
  
  suppressWarnings({
    
    ### normalize the database
    tf_ori_data <- read.delim("./database/TFTargetGene.txt", header = T)
    tf_size_factor <- table(tf_ori_data$TF)
    rec_ori_data <- read.delim("./database/RecTF.txt", header = T)
    rec_size_factor <- table(rec_ori_data$Receptor)
    lig_ori_data <- read.delim("./database/LigRec.txt", header = T)
    lig_size_factor <- table(lig_ori_data$Ligand)
    expres <- (deg$foldChange)
    names(expres) <- deg$Gene
    expres[expres == Inf] <- max(expres[is.finite(expres)])
    expres[expres == -Inf] <- min(expres[is.finite(expres) & expres != 0])
    
    expres[expres < 1] <- (1/expres[expres < 1])
    
    expres[expres == Inf] <- max(expres[is.finite(expres)])
    expres[expres == -Inf] <- min(expres[is.finite(expres) & expres != 0])
    
    rec_exp <- setNames((lr_glom_normal[!duplicated(lr_glom_normal$Receptor.ApprovedSymbol),][[paste("receptor_",receiver_cell,sep="")]]), unique(lr_glom_normal$Receptor.ApprovedSymbol))

    
    lig_col <- 17:(17+ncol(lig_mat)-1)
    lig_exp <-unique((lr_glom_normal[,lig_col]))
    rownames(lig_exp) <- unique(lr_glom_normal$Ligand.ApprovedSymbol)
    
    
    if (is.null(proximity_order)){
      proximity_order <- colnames(lig_mat)
    }
    
    colnames(lig_exp) <- colnames(lig_mat)
    lig_exp <- lig_exp[,proximity_order]
    lig_mat <- as.matrix(lig_mat)
    row.names(lig_mat) <- lr_glom_normal$Ligand.ApprovedSymbol
    lig_mat <- lig_mat[!duplicated(row.names(lig_mat)),]
    lig_exp <- lig_mat[,proximity_order]
    
    
    KL_DF <- lr_glom_normal[,c("Ligand.ApprovedSymbol","Receptor.ApprovedSymbol","KL")]
    
    
    scMlnet_results$LigRec <- stringr::str_split(scMlnet_results$LigRec, pattern = "_", simplify = T)
    scMlnet_results$RecTF <- stringr::str_split(scMlnet_results$RecTF, pattern = "_", simplify = T)
    scMlnet_results$TFTar <- stringr::str_split(scMlnet_results$TFTar, pattern = "_", simplify = T)
    
    
    scMlnet_results$TFTar <- as.data.frame(scMlnet_results$TFTar)
    scMlnet_results$TFTar$V3 <- expres[scMlnet_results$TFTar$V2]
    
    tf_target <- scMlnet_results$TFTar %>% group_by(V1) %>% summarise(get(TF_fun)(V3),.groups = 'drop')
    colnames(tf_target)[2] <- "aggregate"
    tf_n <- setNames(tf_target$aggregate, as.character(tf_target$V1) )
    
    #tf_n <- table(scMlnet_results$TFTar[,1])
    tf_n <- tf_n/tf_size_factor[names(tf_n)]
    
    scMlnet_results$RecTF <- as.data.frame(scMlnet_results$RecTF)
    scMlnet_results$RecTF$tf_n <- 0
    
    names_tf <- names(tf_n)
    for (tf in seq(tf_n)){
      scMlnet_results$RecTF[scMlnet_results$RecTF$V2 %in% names_tf[tf],]$tf_n <- tf_n[tf]
      
    }
    
    rec_n <- scMlnet_results$RecTF %>% group_by(V1) %>% summarise(get(main_fun)(tf_n),.groups = 'drop')
    colnames(rec_n)[2] <- "aggregate"
    rec_n <- setNames(rec_n$aggregate, as.character(rec_n$V1) )
    
    #multilyby rec expression
    if(method == "KL_rec"){
      rec_exp <- setNames((lr_glom_normal[!duplicated(lr_glom_normal$Receptor.ApprovedSymbol),]$KL), unique(lr_glom_normal$Receptor.ApprovedSymbol))
      rec_exp <- scales::rescale(rec_exp, from= min(rec_n), to = max(rec_n))
      rec_n <- rec_n + rec_exp[names(rec_n)]
    }
    
    if(method == "min_max"){
      rec_exp <- scales::rescale(rec_exp, from= min(rec_n), to = max(rec_n))
      rec_n <- rec_n + rec_exp[names(rec_n)]
    }
    
    if(method == "scale_no_center"){
      rec_n<- scale(rec_n,center= F, scale= TRUE)+scale(rec_exp[names(rec_n)],center= F, scale= TRUE)
      rec_n <- setNames(rec_n[,1], row.names(rec_n))
    }
    
    if(method == "scale_sd"){
      rec_n<- scale(rec_n,center= F, scale= sd(rec_n, na.rm = TRUE))+scale(rec_exp[names(rec_n)],center= F, scale= sd(rec_exp[names(rec_n)] , na.rm = TRUE))
      rec_n <- setNames(rec_n[,1], row.names(rec_n))
    }
    
    if(method == "rescale"){
      rec_n<- scales::rescale(as.vector(rec_n))+scales::rescale(rec_exp[names(rec_n)])
    }
    
    
    
    rec_n <- rec_n/rec_size_factor[names(rec_n)]
    
    scMlnet_results$LigRec <- as.data.frame(scMlnet_results$LigRec)
    scMlnet_results$LigRec$rec_n <- 0

    
    names_rec <- names(rec_n)
    for (rec in seq(rec_n)){
      scMlnet_results$LigRec[scMlnet_results$LigRec$V2 %in% names_rec[rec],]$rec_n <- rec_n[rec]
      
    }
    
    # add KL score
    
    if (method == "KL_rec_lig"){
      merged <- merge(scMlnet_results$LigRec,KL_DF, by.x = c("V1","V2"), by.y= c("Ligand.ApprovedSymbol","Receptor.ApprovedSymbol"))
      merged$KL <- scales::rescale(merged$KL, from= min(scMlnet_results$LigRec$rec_n), to = max(scMlnet_results$LigRec$rec_n))
      scMlnet_results$LigRec$rec_n <- scMlnet_results$LigRec$rec_n+ merged$KL
    }
    
    lig_n <- scMlnet_results$LigRec %>% group_by(V1) %>% summarise(get(main_fun)(rec_n),.groups = 'drop')
    colnames(lig_n) <- c("lig","score")
    lig_n <- setNames(lig_n$score, as.character(lig_n$lig) )
    
    if(method == "min_max"){
      # consider by lig expression
      #lig_exp <- scales::rescale(lig_exp, from= min(lig_n), to = max(lig_n))
      lig_exp <- scales::rescale(as.matrix((lig_exp)), from= min(rec_n), to = max(rec_n))
      weight_v <-  seq(from = 1, to = 0.05, length.out = ncol(lig_exp))
      lig_exp <- lig_exp * weight_v
      lig_exp <- rowMeans(lig_exp)
      
      lig_n <- (lig_n) + (lig_exp[names(lig_n)])
    }
    
    #lig_n <- lig_n/lig_size_factor[names(lig_n)]
    
    lig_rank_all <- lig_n[order(lig_n, decreasing = T)]
    
    
    
    
    ## aggregte pathways
    lig_rank_all <- as.data.frame(lig_rank_all)
    lig_rank_all$pathway <- ""
    
    #ligand_name <- lig_rank_all$Var1 
    for (p in rownames(lig_rank_all)){
      lig_rank_all[ p,]$pathway <- CellChatDB[CellChatDB$ligand == p,]$pathway_name[1]
      
    }
    
    pathway_n <- lig_rank_all %>% group_by(pathway) %>% summarise(get(path_fun)(lig_rank_all ),.groups = 'drop')
    colnames(pathway_n)[2] <- "score"
    
    pathway_n <- pathway_n[order(pathway_n$score, decreasing = T),]
    pathway_n$score <- (pathway_n$score/sum(pathway_n$score))*100
    pathway_n$score.perc <- scales::percent((pathway_n$score/sum(pathway_n$score)),
                                            accuracy= 0.1)
    

    
    # Package
    library(treemap)
    #library(randomcoloR)
    n <- length(unique(lig_rank_all$pathway))
    # # palette <- distinctColorPalette(n)
    # # 
    # # palette <- setNames(distinctColorPalette(length(path_included)), path_included)
    # # dput(palette)
    # palette <-c(TGFb = "#DBC6DE", NRG = "#E2E5D5", BMP10 = "#DA3BE0", GDF = "#679C8F", 
    #             GDNF = "#E2E2A1", NODAL = "#E24C66", ACTIVIN = "#C2E579", WNT = "#DA77E8", 
    #             ncWNT = "#64B1DF", EGF = "#66E1E1", BMP = "#E88E50", FGF = "#8580DF", 
    #             PDGF = "#62E4B7", VEGF = "#AAD8E1", IGF = "#D5AB90", INSULIN = "#7945DF", 
    #             HH = "#6EE07F", EDA = "#E2A1E0", NGF = "#8D984C", NT = "#DE8D9F", 
    #             FLT3 = "#70E041", HGF = "#D9E83D", NOTCH = "#8D92BF", NRXN = "#CC4E99", 
    #             OCLN = "#E1C959", ROBO = "#B0E7BE")
    
    pathway_n$pathway <- as.factor(pathway_n$pathway)
    palette <- palette[levels(pathway_n$pathway)]
    #palette <- as.factor(palette)
    
    pdf(paste(output_dir,"/","Tree_path_",receiver_cell,".pdf", sep = ""), 6,7.5) # 6,7.5
    
    if (show.sub.pathway ==F){
      treemap::treemap(pathway_n,
                       
                       # data
                       index=c("pathway", "score.perc"),
                       vSize="score",
                       type="categorical",
                       vColor = "pathway",
                       
                       
                       # Main
                       title="",
                       palette= palette ,
                       
                       # Borders:
                       border.col=c("black", "black"),             
                       border.lwds= c(3,1), 
                       mirror.y = mirror.y,
                       
                       # Labels
                       fontsize.labels=c(0.6,0.6),      #0.5,
                       fontcolor.labels=c("black","black"),
                       fontface.labels=c(2,1),            
                       bg.labels=c("transparent"),              
                       align.labels= list(c("left", "top"),c("right","bottom")),                                  
                       overlap.labels=0,
                       inflate.labels=T ,                       # If true, labels are bigger when rectangle is bigger.
                       position.legend =  "none"
                       
      )
      dev.off()
      show.sub.pathway = T
    }
    
    if (show.sub.pathway ==T){
      pdf(paste(output_dir,"/","Tree_path_",receiver_cell,"_sub.pdf", sep = ""), 6,7.5)
      lig_rank_all$lig <- row.names(lig_rank_all)
      lig_rank_all <- merge(lig_rank_all, pathway_n, by.x ="pathway", by.y= "pathway")
      lig_rank_all$pathway <- as.factor(lig_rank_all$pathway)
      treemap::treemap(lig_rank_all,
                       
                       # data
                       index=c("pathway","score.perc", "lig"),
                       vSize="lig_rank_all",
                       type="categorical",
                       vColor = "pathway",
                       
                       
                       # Main
                       title="",
                       palette= palette,
                       
                       # Borders:
                       border.col=c("black", "black", "white"),             
                       border.lwds= c(3,1,3),  
                       mirror.y = mirror.y,
                       
                       # Labels
                       fontsize.labels=c(0.6,0,0.6),      #0.5,
                       fontcolor.labels=c("black","black","white"),
                       fontface.labels=2,            
                       bg.labels=c("transparent"),              
                       align.labels= list(c("left", "top"),c("right", "top"),c("center","bottom")),                                  
                       overlap.labels=0,
                       inflate.labels=T ,                       # If true, labels are bigger when rectangle is bigger.
                       position.legend =  "none"
                       
      )
    }
    
    dev.off()
  })
}
