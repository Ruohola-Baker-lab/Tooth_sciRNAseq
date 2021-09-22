
# load the customized heatmap function.
source("R/custom_Heatmap.R")


# load your cell dataset object (cds)
amelo_odonto_cds <- readRDS("sample_datasets/dental_epithelial_lineage_cds.rds")


# custom_Heatmap: only supply your cds,  a clustering_column & age_group_column
# you can specify the markers genes by selected_markers = c("DSPP", "AMBN")
# or leave it to automaticly find markers by selected_markers = "auto", then set the number_of_markers
# or do both selected_markers = c("auto", "DSPP","NES","DMP1","FGF3") to get top genes + yoursupplied genes. Note supllied genes will be included in goterm analysis
# number_of_genes2goterm define the number of top genes to  use for goterm analysis. You have to find the sweet spot number. not too many and you dliute the result or too little and you miss important things
# gene_font_size: size of the labeled genes
# keywords_font_size: size of goterm key words
# number_of_keywords: number of top keywords to display.
# keywords_blacklist: this manually remove redundant, meaningless or irrleveant words from the goterm keywords
# load_previous is to load all the saved metadata if present from previous run on same cds, so the second run on the same dataset can be faster. if you made alot of changes then select FALSE and calculate everyting again. Some variables won't do any effect if you change them unless you do a complete rerun
# to be more spcific changes to number_of_genes2goterm variable will update only if you not load prevous metadata
# all results are saved in custome_Heatmap_output folder
# it will tell you how many minutes it took to run at the end. This data take about 10min in first run


###Example
custom_Heatmap(
  cds = amelo_odonto_cds,
  clustering_column = "assigned_cell_type",
  age_clustering_column = "age_group",
  selected_markers = c("auto"),
  number_of_markers = 3, 
  number_of_genes2goterm = 50,
  gene_font_size= 10,
  keywords_font_size= 13, 
  number_of_keywords = 7,
  keywords_blacklist = c(
    "symbiont",
    "head",
    "cardiovascular",
    "transcript",
    "regulated", 
    "cartilage",
    "generation",
    "response",
    "chondrocyte",
    "platelet",
    "sensory",
    "neurons",
    "osteoblast",
    "bone",
    "ureteric",
    "urogenital",
    "renal",
    "relaxation",
    "gland",
    "blood",
    "quality",
    "smooth",
    "ossification",
    "odontogenesis",
    "kidney",
    "head",
    "septum",
    "aortic",
    "taxis",
    "endogenous",
    "valve",
    "locomotion",
    "subcellular",
    "movement",
    "axon",
    "luteolysis",
    "chemical",
    "stimulus",
    "mesonephric",
    "substance",
    "forebrain",
    "homotypic",
    "embryonic",
    "endodermal",
    "drug",
    "nephron",
    "renal",
    "tubule",
    "metanephros",
    "free",
    "ubiquitin",
    "chain",
    "poly",
    "polymerization",
    "synaptic",
    "transmission",
    "spontaneous",
    "neurotransmitter",
    "synapse",
    "vessel",
    "lineage",
    "unfolded",
    "inorganic",
    "estradiol",
    "ketone",
    "lipid",
    "endochondral",
    "topologically",
    "posttranscriptional",
    "cortisol",
    "reverse",
    "signal",
    "cholesterol",
    "interleukin",
    "endoplasmic",
    "incorrect",
    "axonogenesis",
    "sterol",
    "motor",
    "telomere",
    "peptide",
    "folding",
    "anion",
    "reticulum",
    "recognition",
    "programmed",
    "death",
    "guidance",
    "hydroxylation",
    "plasma",
    "transduction",
    "cycle",
    "innervation",
    "organelle",
    "presynaptic",
    "vasculature",
    "heterophilic",
    "chemotaxis",
    "tripeptide",
    "assembly",
    "odontoblast",
    "telencephalon",
    "amyloid",
    "fibril",
    "structure",
    "oligopeptide",
    "artery",
    "pallium",
    "aorta",
    "cerebral",
    "cortex",
    "axonal",
    "fasciculation"
    
    
  ),
  auto_sort_columns = F,
  load_previous = T,
  normalize_by_cluster = F

)
