#### Load necessary packages ####

# If you do not already have these packages
# you may need to install each of them with the command "install.packages()": 
install.packages('dplyr')
library('dplyr')
library(Seurat)
library(ggplot2)
library(patchwork)
library("berryFunctions")
library(collections)
library(tibble)
library(janitor)
library(Metrics)
library(reshape)
library(scales)
library('tidyr')
library(HSMMSingleCell)
library(SingleCellExperiment)
library('scater')
library(destiny)
library(BisqueRNA)
library(plotly)
library("scatterplot3d")
library(monocle)
library(ComplexHeatmap)
library(stringr)
library(reshape)
library(data.table)
library(CATALYST)

#### Importing Seurat object and viewing UMAP ####

# importing Sujit's Seurat object

# For Qingkang: change to your local file location of CD4T_CAVA.rds
CD4_cava <- readRDS("/Users/jmakings/Desktop/Antoine_stuff/AntoinePsuedotime/CD4T_CAVA.rds")

# This will give you an overview of the Seurat object
# Should show that RNA and ADT assays are present, and
# 2 dimensionality reductions have been performed: pca, and umap
CD4_cava

# Because umap has already been preformed, we can plot it with this command
UMAPPlot(CD4_cava, group.by = 'Cluster')

#### Destiny Diffusion Map analysis and plotting ####

# set seed so that pseudotime graphs always appear in same direction
set.seed(2)

##### Preprocessing and computations to get Diffusion Map #####

# We will follow the workflow from the Destiny documentation:  
# https://theislab.github.io/destiny/reference/DPT.html

# converting to Single Cell Experiment (necessary format to work with Destiny)
ssCD4_cava <- as.SingleCellExperiment(CD4_cava)

# Create diffusion map object from logcounts data, with 50 principal components
# The DiffusionMap function may take a few minutes to run, be patient
matrix <- as.data.frame(assay(ssCD4_cava, i='logcounts'))
dm2 <- DiffusionMap(t(matrix), n_pcs = 50)

# Add eigenvectors of diffusion map to Single Cell Experiment object
reducedDim(ssCD4_cava, type = 'DC') <- dm2@eigenvectors

# This will give a plot of DC1 and DC2, with cells labeled by cluster
dcPlot1_2 <- plotReducedDim(ssCD4_cava, dimred = 'DC', ncomponents = 1:2,colour_by = 'Cluster')
dcPlot1_2

# This creates a pseudo time object from the diffusion map's transition probabilities
dpt <- DPT(dm2)

# this will show specific branches and a path between branches 
branchPlot <- plot(dpt, root = 2, paths_to = c(1,3), col_by = 'branch')
branchPlot

# loading top DCs into single cell experiment
ssCD4_cava$dc1 <- dm2$DC1
ssCD4_cava$dc2 <- dm2$DC2
ssCD4_cava$dc3 <- dm2$DC3

##### Filtering and plotting diffusion plot with different methods #####

# To filter Single Cell Experiment by clusters 1,2,7,14,and 17, which contained 
# almost all Treg and exTreg cells
reduced <- filterSCE(ssCD4_cava, ssCD4_cava$Cluster=='CD4T_7' | 
                       ssCD4_cava$Cluster=='CD4T_17' | ssCD4_cava$Cluster=='CD4T_1' | 
                       ssCD4_cava$Cluster=='CD4T_2' | ssCD4_cava$Cluster=='CD4T_14' )
plotReducedDim(reduced, dimred = 'DC', ncomponents = 1:2,colour_by = 'Cluster')

# To show just clusters 7 and 17, we can use a similar method
seven17 <- filterSCE(ssCD4_cava, ssCD4_cava$Cluster=='CD4T_7' | ssCD4_cava$Cluster=='CD4T_17' )
plotReducedDim(seven17, dimred = 'DC', ncomponents = 1:2,colour_by = 'Cluster')

# Different color palletes for creating graphs
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
palette16 <- c("#000000","#999999","#004949","#009292","#ff6db6","#ffb6db",
               "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
               "#920000","#924900","#db6d00","#24ff24","#ffff6d")
grayPalette <- c("gray","gray","gray","gray","gray","gray","gray","gray","blue","gray","gray","gray","gray","gray","red", "gray")

fiveClusters <- ggplot(as.data.frame(colData(reduced)),
                       aes(x = dc1,
                           y = dc2, color = Cluster)) + geom_point() + theme_classic() +
  xlab("DC1") + ylab("DC2") + scale_colour_manual(values=cbPalette) + 
  ggtitle("Destiny DC Plot by Cluster")

allClusters <- ggplot(as.data.frame(colData(ssCD4_cava)),
                      aes(x = dc1,
                          y = dc2, color = Cluster)) + geom_point() + theme_classic() +
  xlab("DC1") + ylab("DC2") + scale_colour_manual(values=palette16) + 
  ggtitle("Destiny DC Plot by Cluster")

# different colorings of the DC plots, originally created for colorblind viewing
fiveClusters
allClusters

alpha_vec <- vector()
for (x in ssCD4_cava$Cluster) {
  if (x == "CD4T_7" || x == "CD4T_17") {
    alpha_vec <- append(alpha_vec, 1)
  }
  else 
    alpha_vec <- append(alpha_vec, 0.4)
}

grayPlot <- ggplot(as.data.frame(colData(ssCD4_cava)),
                   aes(x = dc1,
                       y = dc2, color = Cluster)) + geom_point(alpha = alpha_vec) + theme_classic() +
  xlab("DC1") + ylab("DC2") + scale_colour_manual(values=grayPalette) + 
  scale_fill_manual(values = alpha(alpha_vec)) + 
  ggtitle("Destiny DC Plot by Cluster")

# This plot colors only cluster 7 and 17 and grays out the rest of the clusters
grayPlot

##### Adding Bin Lines to split cells by DC1 value ##### 

# function to get bin lines, include the plot to put over, the number of bins, and single cell experiment object
plotBinLines <- function(plot, bins, ssObj) {
  binCD4_ss <- cut(ssObj$dc1, breaks=bins)
  binChange <- str_replace_all(binCD4_ss, pattern = ".*,([^-]*)]*", replacement = "\\1")
  binChange2 <- gsub(']','',binChange)
  binChange2 <- sort(as.numeric(unique(binChange2)))
  dif <- binChange2[2] -binChange2[1]
  
  bin20plot <- plot + geom_vline(xintercept = c(unique(binChange2),binChange2[1]-dif))
  return(bin20plot)
}

# plot of bin lines
binPlot <- plotBinLines(grayPlot,10, ssCD4_cava)
binPlot

#### Differential Expression Analysis #### 

# get bins then add to Single Cell Experiment object
binCD4_ss <- cut(ssCD4_cava$dc1, breaks=10, labels=c(1:10))
ssCD4_cava$Bins <- binCD4_ss

newSrt <- as.Seurat(ssCD4_cava)

# setting ident to cluster, but will later set it to bin and cluster together 
newSrt <- SetIdent(newSrt, value=newSrt@meta.data$Cluster)
newSrt <- SetIdent(newSrt, value = newSrt@meta.data$Bins)

# scaling marker and trancript data 
markers <- rownames(newSrt)
newSrt <- ScaleData(newSrt, features = markers)
transcripts <- rownames(newSrt@assays$RNA@counts)
newSrt <- ScaleData(newSrt, assay = 'RNA', features = transcripts)

exTregBinsDE <- function(Srt) {
  l <- list()
  for (i in 1:7) {
    print(i)
    l[[i]] <- FindMarkers(subset(x = Srt, Cluster == 'CD4T_7'), ident.1 =i,indent.2=c(8,9,10), min.pct=0, logfc.threshold = 0)
  }
  return(l)
}

#Cluster 7 DE of each bin vs exTreg cluster on the right (Bins 8,9,10)
binsList <- exTregBinsDE(newSrt)

#Cluster 7 DE of each bin vs all other bins
allBins <- FindAllMarkers(subset(x = newSrt, Cluster == 'CD4T_7'), min.pct=0, logfc.threshold = 0)

newSrt <- SetIdent(newSrt, value = newSrt@meta.data$Bins)

#Cluster 7 DE Bins 1-7 vs Bins 8-10
leftvsRight <- FindMarkers(subset(x = newSrt, Cluster == 'CD4T_7'), ident.1 =1:7,indent.2=c(8,9,10), min.pct=0, logfc.threshold = 0)

# add bin and cluster together
newSrt$binCluster <- str_c(newSrt$Cluster, "_Bin_",newSrt$Bins)
counts <- data.frame(t(newSrt@assays$RNA@counts), newSrt$Cell_Index, newSrt$binCluster)

# set it as ident 
newSrt <- SetIdent(newSrt, value = newSrt@meta.data$binCluster)

# finding DEs for each bin/cluster pairing (For 170 bin/cluster pairs this will take some time)
binClustervsRest <- FindAllMarkers(newSrt, min.pct=0, logfc.threshold = 0,min.cells.feature = 0,min.cells.group = 0)

#### Heatmaps from Marker Differential Expression Data ####

# subset for just clusters 7 and 17
newSrt <- SetIdent(newSrt, value = newSrt@meta.data$Cluster)
Cluster.7.17 <- subset(x = newSrt, idents = c('CD4T_7','CD4T_17' ))
newSrt <- SetIdent(newSrt, value = newSrt@meta.data$binCluster)

# order the clusters accordingly
clusterLevels <- c('CD4T_7_Bin_1','CD4T_17_Bin_1','CD4T_7_Bin_2',
                   'CD4T_17_Bin_2','CD4T_7_Bin_3','CD4T_17_Bin_3',
                   'CD4T_7_Bin_4', 'CD4T_17_Bin_4', 'CD4T_7_Bin_5',
                   'CD4T_7_Bin_6', 'CD4T_7_Bin_7','CD4T_17_Bin_7',
                   'CD4T_7_Bin_8', 'CD4T_17_Bin_8','CD4T_7_Bin_9', 
                   'CD4T_17_Bin_9', 'CD4T_7_Bin_10')

# change active idents and assay
Cluster.7.17 <- SetIdent(Cluster.7.17, value = Cluster.7.17@meta.data$binCluster)
Cluster.7.17@active.assay <- 'ADT'

# depending on which direction psuedotime goes (either right or left) in the DC 
# plot, the bin/cluster pairs may differ. The 'whichClusters' function determines 
# which bin/cluster pairs are present and orders them accordingly
binclusters.7.17 <- unique(Idents(Cluster.7.17 ))

whichClusters <- function(binclusters) {
  if ('CD4T_17_Bin_1' %in% binclusters == TRUE)
    clusterLevels <- c('CD4T_7_Bin_1','CD4T_17_Bin_1','CD4T_7_Bin_2',
                       'CD4T_17_Bin_2','CD4T_7_Bin_3','CD4T_17_Bin_3',
                       'CD4T_7_Bin_4', 'CD4T_17_Bin_4', 'CD4T_7_Bin_5',
                       'CD4T_7_Bin_6', 'CD4T_7_Bin_7','CD4T_17_Bin_7',
                       'CD4T_7_Bin_8', 'CD4T_17_Bin_8','CD4T_7_Bin_9', 
                       'CD4T_17_Bin_9', 'CD4T_7_Bin_10')
  else 
    clusterLevels <- c('CD4T_7_Bin_1','CD4T_7_Bin_2','CD4T_17_Bin_2',
                       'CD4T_7_Bin_3','CD4T_17_Bin_3',
                       'CD4T_7_Bin_4', 'CD4T_17_Bin_4', 'CD4T_7_Bin_5',
                       'CD4T_7_Bin_6', 'CD4T_7_Bin_7','CD4T_17_Bin_7',
                       'CD4T_7_Bin_8', 'CD4T_17_Bin_8','CD4T_7_Bin_9', 
                       'CD4T_17_Bin_9', 'CD4T_7_Bin_10','CD4T_17_Bin_10')
  return(clusterLevels)
}

bincluster.7.17order <- whichClusters(binclusters.7.17)

# ordering the heatmap so bins are in order
Cluster.7.17@active.ident <- factor(x=Cluster.7.17@active.ident, levels = bincluster.7.17order)
Cluster.7.17@active.assay <- 'ADT'

# plotting top 10 log2FC Differential Expression for each cluster
binClustervsRest %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC) -> top10

# finally, a marker heatmap for all bins in cluster 7 and 17 
MarkerHeatmapAllBins <- DoHeatmap(Cluster.7.17, features = top10$gene, angle = 90, size = 1.5)
MarkerHeatmapAllBins

TregIdents = c('CD4T_17_Bin_1','CD4T_17_Bin_2',
               'CD4T_17_Bin_3',  'CD4T_17_Bin_4', 'CD4T_17_Bin_5',
               'CD4T_17_Bin_6',
               'CD4T_17_Bin_7', 'CD4T_17_Bin_8', 'CD4T_17_Bin_9','CD4T_17_Bin_10')
exTregIdents = c('CD4T_7_Bin_1','CD4T_7_Bin_2', 
                 'CD4T_7_Bin_3','CD4T_7_Bin_4'
                 ,'CD4T_7_Bin_5','CD4T_7_Bin_6', 
                 'CD4T_7_Bin_7', 'CD4T_7_Bin_8', 
                 'CD4T_7_Bin_9', 'CD4T_7_Bin_10')

CustomHeatmap <- function(Srt, idents, features, binclusters, size) {
  new_idents = c()
  for (i in idents) {
    if (i %in% binclusters)
      new_idents <- append(new_idents, i)
  }
  subset_data <- subset(x = Srt, idents = new_idents)
  heatmap1 <- DoHeatmap(subset_data, features = features, angle = 90, size = size)
  return(heatmap1)
}

# Marker heatmap just for cluster 17
tregMarkerHeatmap <- CustomHeatmap(Cluster.7.17, TregIdents, top10$gene, binclusters.7.17, 2)
tregMarkerHeatmap

# Marker heatmap just for cluster 7
extregMarkerHeatmap <- CustomHeatmap(Cluster.7.17, exTregIdents, top10$gene, binclusters.7.17, 2)
extregMarkerHeatmap 

#### Heatmaps from Transcript Differential Expression Data ####

# First, get differential expression for transcripts (done by specifying 'assay="RNA" ')
# because there are many more transcripts than markers, this computation takes especially long
binClusterGenes <- FindAllMarkers(newSrt, assay = "RNA",min.pct=0, logfc.threshold = 0,min.cells.feature = 0,min.cells.group = 0)

binClusterDF <- as.data.frame(binClusterGenes)

#just cluster7 genes as a function of bin, for Klaus (8/30/22) 
cluster7genes <- binClusterDF[(binClusterDF$cluster=='CD4T_7_Bin_1' |
                                 binClusterDF$cluster=='CD4T_7_Bin_2' |
                                 binClusterDF$cluster=='CD4T_7_Bin_3' |
                                 binClusterDF$cluster=='CD4T_7_Bin_4' |
                                 binClusterDF$cluster=='CD4T_7_Bin_5' |
                                 binClusterDF$cluster=='CD4T_7_Bin_6' |
                                 binClusterDF$cluster=='CD4T_7_Bin_7' |
                                 binClusterDF$cluster=='CD4T_7_Bin_8' |
                                 binClusterDF$cluster=='CD4T_7_Bin_9' |
                                 binClusterDF$cluster=='CD4T_7_Bin_10'),]

#just cluster 17 genes as a function of bin
cluster17genes <- binClusterDF[(binClusterDF$cluster=='CD4T_17_Bin_1' |
                                 binClusterDF$cluster=='CD4T_17_Bin_2' |
                                 binClusterDF$cluster=='CD4T_17_Bin_3' |
                                 binClusterDF$cluster=='CD4T_17_Bin_4' |
                                 binClusterDF$cluster=='CD4T_17_Bin_5' |
                                 binClusterDF$cluster=='CD4T_17_Bin_6' |
                                 binClusterDF$cluster=='CD4T_17_Bin_7' |
                                 binClusterDF$cluster=='CD4T_17_Bin_8' |
                                 binClusterDF$cluster=='CD4T_17_Bin_9' |
                                 binClusterDF$cluster=='CD4T_17_Bin_10'),]

# cluster 7 and 17 genes
filtered <- binClusterDF[(binClusterDF$cluster=='CD4T_7_Bin_1' |
                            binClusterDF$cluster=='CD4T_7_Bin_2' |
                            binClusterDF$cluster=='CD4T_7_Bin_3' |
                            binClusterDF$cluster=='CD4T_7_Bin_4' |
                            binClusterDF$cluster=='CD4T_7_Bin_5' |
                            binClusterDF$cluster=='CD4T_7_Bin_6' |
                            binClusterDF$cluster=='CD4T_7_Bin_7' |
                            binClusterDF$cluster=='CD4T_7_Bin_8' |
                            binClusterDF$cluster=='CD4T_7_Bin_9' |
                            binClusterDF$cluster=='CD4T_7_Bin_10' |
                            binClusterDF$cluster=='CD4T_17_Bin_1' |
                            binClusterDF$cluster=='CD4T_17_Bin_2' |
                            binClusterDF$cluster=='CD4T_17_Bin_3' |
                            binClusterDF$cluster=='CD4T_17_Bin_4' |
                            binClusterDF$cluster=='CD4T_17_Bin_5' |
                            binClusterDF$cluster=='CD4T_17_Bin_6' |
                            binClusterDF$cluster=='CD4T_17_Bin_7' |
                            binClusterDF$cluster=='CD4T_17_Bin_8' |
                            binClusterDF$cluster=='CD4T_17_Bin_10' |
                            binClusterDF$cluster=='CD4T_17_Bin_9'),]

# filter for top differentiall expressed genes
filtered %>% 
  group_by(cluster) %>% 
  top_n(n = 15, wt = avg_log2FC) -> topgenes

# Heatmap showing all cluster 7 and 17 transcripts across bins 
TranscriptHeatmapAllBins <- DoHeatmap(Cluster.7.17, assay = 'RNA', features = topgenes$gene, angle = 90, size = 2) +
  theme(text=element_text(size = 4), legend.text = element_text(size = 5))
TranscriptHeatmapAllBins

# subsetting for Cluster 7 and Cluster 17 individually
Cluster.7.17 <- SetIdent(Cluster.7.17, value = Cluster.7.17@meta.data$Cluster)
Cluster.7 <- subset(x = Cluster.7.17, idents = c('CD4T_7'))
Cluster.17 <- subset(x = Cluster.7.17, idents = c('CD4T_17'))
Cluster.7 <- SetIdent(Cluster.7, value = Cluster.7@meta.data$binCluster)
Cluster.17 <- SetIdent(Cluster.17, value = Cluster.17@meta.data$binCluster)

Cluster.7.17 <- SetIdent(Cluster.7.17, value = Cluster.7.17@meta.data$binCluster)

# For cluster 7 bins only
TranscriptHeatmapC7 <- DoHeatmap(Cluster.7, assay = 'RNA', features = topgenes$gene, angle = 90, size = 2) +
  theme(text=element_text(size = 4), legend.text = element_text(size = 5))
TranscriptHeatmapC7

# For cluster 17 bins only
TranscriptHeatmapC17 <- DoHeatmap(Cluster.17, assay = 'RNA', features = topgenes$gene, angle = 90, size = 2) +
  theme(text=element_text(size = 4), legend.text = element_text(size = 5))
TranscriptHeatmapC17

##### Create transcript Pseudobulk heatmap from Sujit's way (standard Ley workflow)##### 

# This code chunk for preprocessing taken from Sujit
gene_df <- as.data.frame(newSrt@assays$RNA@data) %>% t()
gene_meta <- newSrt@meta.data %>% select(binCluster)
gene_merge <- merge(gene_df, gene_meta, by=0)
gene_avg <- gene_merge %>% group_by(binCluster) %>% summarise_all(funs(mean)) %>% data.table()
gene_avg$Row.names <- NULL
gene_avg$binCluster <- gsub('CD4T','C',gene_avg$binCluster)
gene_avg <- gene_avg %>% t()
gene_avg <- gene_avg %>% row_to_names(row_number = 1)

# This code chunk for preprocessing taken from Sujit
geneScaled <- CreateSeuratObject(gene_avg)
genedata <- as.data.frame(t(as.matrix(GetAssayData(geneScaled))))
scale_limits <- c(-2,2)
scaled_gene <- data.frame(lapply(genedata,function(x)rescale(x,to=scale_limits)))
rownames(scaled_gene) <-  rownames(genedata)
scaled_gene <- t(scaled_gene)
geneScaled@assays$RNA@scale.data <- scaled_gene
geneScaled = as.matrix(geneScaled@assays$RNA@scale.data)

pal_cols <- colorRampPalette(c("blue","yellow","red"))

geneScaled <- as.data.frame(geneScaled)

# co_1 and co_2 are the two options for bin/cluster pairing
# based on the orientation of the diffusion map
co_1 <- c('C_7_Bin_1','C_7_Bin_2', 
          'C_7_Bin_3','C_7_Bin_4'
          ,'C_7_Bin_5','C_7_Bin_6', 
          'C_7_Bin_7', 'C_7_Bin_8', 
          'C_7_Bin_9', 'C_17_Bin_1', 'C_17_Bin_2',
          'C_17_Bin_3',  'C_17_Bin_4', 
          'C_17_Bin_7', 'C_17_Bin_8', 'C_17_Bin_9', 
          'C_7_Bin_10')

co_2 <- c('C_7_Bin_1','C_7_Bin_2',
          'C_17_Bin_2','C_7_Bin_3','C_17_Bin_3',
          'C_7_Bin_4', 'C_17_Bin_4', 'C_7_Bin_5',
          'C_7_Bin_6', 'C_7_Bin_7','C_17_Bin_7',
          'C_7_Bin_8', 'C_17_Bin_8','C_7_Bin_9', 
          'C_17_Bin_9', 'C_7_Bin_10', 'C_17_Bin_10')

# subset scaled.7.17 with co_1 or co_2 based on which bin/cluster pairs are present
clusters <- co_2
#clusters <- co_1

# if the below line throws an error 'undefined columns selected', then uncomment: 
# clusters <- co_1
scaled.7.17 <- geneScaled[clusters]

newnames_row1 <- lapply(
  rownames(scaled.7.17) ,
  function(x) bquote(bold(.(x))))

newnames_cols1 <- lapply(
  colnames(scaled.7.17),
  function(x) bquote(bold(.(x))))

# general use of pheatmap for all genes, which clusters based on bin/cluster pairing
pheatmap::pheatmap(scaled.7.17, color = pal_cols(100), cluster_rows = FALSE, 
                   main = "Scaled Heatmap for Cluster 7 and 17 Bulk Transcriptome, by Bin", 
                   fontsize = 10, fontsize_col = 10, fontsize_row = 2, fontsize_number=2, 
                   show_rownames = T, 
                   labels_row = as.expression(newnames_row1),
                   labels_col = as.expression(newnames_cols1))


# heatmap for every gene, in bin order
gene_pseudo_heatmap <- Heatmap(scaled.7.17, column_order = clusters, cluster_rows=T, show_row_names = F,
        row_names_gp = gpar(fontsize=6))
gene_pseudo_heatmap

# now we will look at a heatmaps with only the top differentially expressed genes
scaledDown <- scaled.7.17[row.names(scaled.7.17) %in% topgenes$gene,]

newnames_row <- lapply(
  rownames(scaledDown) ,
  function(x) bquote(bold(.(x))))

newnames_cols <- lapply(
  colnames(scaledDown),
  function(x) bquote(bold(.(x))))

# pheatmap with only the top differentially expressed genes, with bin/cluster pairs clustered
pheatmap::pheatmap(scaledDown, color = pal_cols(100), cluster_rows = FALSE, 
                   main = "Scaled Heatmap for Cluster 7 and 17 Bulk Transcriptome, by Bin", 
                   fontsize = 10, fontsize_col = 10, fontsize_row = 5.5, fontsize_number=2, 
                   labels_row = as.expression(newnames_row),
                   labels_col = as.expression(newnames_cols))

scaledDown <- as.matrix(scaledDown)

# This heatmap has top DE genes, with bin/cluster pairs clustered
Heatmap(scaledDown, cluster_columns=T, cluster_rows=T, 
        row_names_gp = gpar(fontsize=6))

# This heatmap has top DE genes, with bin/cluster pairs ordered from 1 to 10
Heatmap(scaledDown, column_order=clusters, cluster_rows=T, 
        row_names_gp = gpar(fontsize=6))

#### Create Protein/Marker Pseudobulk heatmap from Sujit's way (standard Ley workflow) ####

# This code chunk for preprocessing taken from Sujit
norm_df <- as.data.frame(newSrt@assays$ADT@data) %>% t()
norm_meta <- newSrt@meta.data %>% select(binCluster)
norm_merge <- merge(norm_df, norm_meta, by=0)
norm_avg <- norm_merge %>% group_by(binCluster) %>% summarise_all(funs(mean)) %>% data.table()
norm_avg$Row.names <- NULL
norm_avg$binCluster <- gsub('CD4T','C',norm_avg$binCluster)
norm_avg <- norm_avg %>% t()
norm_avg <- norm_avg %>% row_to_names(row_number = 1)

# This code chunk for preprocessing taken from Sujit
normScaled <- CreateSeuratObject(norm_avg)
data <- as.data.frame(t(as.matrix(GetAssayData(normScaled))))
scale_limits <- c(-2,2)
scaled_data <- data.frame(lapply(data,function(x)rescale(x,to=scale_limits)))
rownames(scaled_data) <-  rownames(data)
scaled_data <- t(scaled_data)
normScaled@assays$RNA@scale.data <- scaled_data
normScaled = as.matrix(normScaled@assays$RNA@scale.data)

pal_cols <- colorRampPalette(c("blue","yellow","red"))

normScaled <- as.data.frame(normScaled)

scaled.7.17 <- normScaled[clusters]

newnames_row <- lapply(
  rownames(scaled.7.17) ,
  function(x) bquote(bold(.(x))))

newnames_cols <- lapply(
  colnames(scaled.7.17),
  function(x) bquote(bold(.(x))))

# Marker heatmap clustered by bin/cluster pairing
pheatmap::pheatmap(scaled.7.17, color = pal_cols(100), cluster_rows = FALSE,
                   main = "Scaled Heatmap for Cluster 7 and 17 Bulk Marker Expression, by Bin", fontsize = 10, fontsize_col = 10, fontsize_row = 10, fontsize_number=10, 
                   labels_row = as.expression(newnames_row),
                   labels_col = as.expression(newnames_cols))

scaled.7.17 <- as.matrix(scaled.7.17)

# Heatmap ordered by bin/cluster pairing
Heatmap(scaled.7.17, column_order = clusters, cluster_rows=T)

# Heatmap with clustering of bin/cluster pair
Heatmap(scaled.7.17, cluster_columns = T, cluster_rows=T)










