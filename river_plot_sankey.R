# Libraries
library(tidyverse)
library(viridis)
library(patchwork)
library(hrbrthemes)
library(circlize)

# Load dataset from github
data <- read.table("https://raw.githubusercontent.com/holtzy/data_to_viz/master/Example_dataset/13_AdjacencyDirectedWeighted.csv", header=TRUE)
# Package
library(networkD3)

# I need a long format
data_long <- data %>%
  rownames_to_column %>%
  gather(key = 'key', value = 'value', -rowname) %>%
  filter(value > 0)
colnames(data_long) <- c("source", "target", "value")
data_long$target <- paste(data_long$target, " ", sep="")

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(name=c(as.character(data_long$source), as.character(data_long$target)) %>% unique())

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
data_long$IDsource=match(data_long$source, nodes$name)-1 
data_long$IDtarget=match(data_long$target, nodes$name)-1

# prepare colour scale
ColourScal ='d3.scaleOrdinal() .range(["#FDE725FF","#B4DE2CFF","#6DCD59FF","#35B779FF","#1F9E89FF","#26828EFF","#31688EFF","#3E4A89FF","#482878FF","#440154FF"])'

# Make the Network
sankeyNetwork(Links = data_long, Nodes = nodes,
              Source = "IDsource", Target = "IDtarget",
              Value = "value", NodeID = "name", 
              sinksRight=FALSE, colourScale=ColourScal, nodeWidth=40, fontSize=13, nodePadding=20)



##### my data
refquery_day_16 <- subset(x = refquery, subset = projection == c("day_16", NA))

data_summary <- seurat_OE_Am_sub@meta.data[,c("preAnnotation","new_projection")] %>% group_by(preAnnotation,new_projection ) %>% count() %>% mutate(freq = n / sum(n))

data_summary <- seurat_OE_Am_sub@meta.data[,c("preAnnotation","new_projection")] %>% count(preAnnotation, new_projection) %>% mutate(freq = n / sum(n))

data_summary <- seurat_OE_Am_sub@meta.data[,c("projection","new_projection")] %>% count(projection, new_projection) %>% mutate(freq = n / sum(n))

seurat_OE_Am_sub_day_16 <- subset(x = seurat_OE_Am_sub, subset = projection == c("day_16", NA))

data_summary <- seurat_OE_Am_sub_day_16@meta.data[,c("preAnnotation","new_projection")] %>% count(preAnnotation, new_projection) %>% mutate(freq = n / sum(n))


data_summary <- seurat_OE_Am_sub_day_16@meta.data[,c("clustering","new_projection")] %>% count(clustering, new_projection) %>% mutate(freq = n / sum(n))
data_summary <- seurat_OE_Am_sub_day_16@meta.data[,c("projection","new_projection")] %>% count(projection, new_projection) %>% mutate(freq = n / sum(n))




# filter the data if needed

data_summary<- data_summary[data_summary$freq >0.01,]
data_summary<- data_summary[data_summary$n > 25,]


# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(name=c(as.character(data_summary$clustering), as.character(data_summary$new_projection)) %>% unique())

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
data_summary$IDsource=match(data_summary$clustering, nodes$name)-1 
data_summary$IDtarget=match(data_summary$new_projection, nodes$name)-1


# prepare colour scale
ColourScal ='d3.scaleOrdinal() .range(["#FDE725FF","#B4DE2CFF","#6DCD59FF","#35B779FF","#1F9E89FF","#26828EFF","#31688EFF","#F8766D","#B79F00","#00BA38","#00BFC4","#619CFF","#F564E3"])'

# Make the Network
sankeyNetwork(Links = data_summary, Nodes = nodes,
              Source = "IDsource", Target = "IDtarget",
              Value = "n", NodeID = "name", 
              sinksRight=FALSE, colourScale=ColourScal, nodeWidth=40, fontSize=13, nodePadding=20)


table(data_summary$new_projection)
