
library(Seurat)
library(Matrix)
library(patchwork)
library(dplyr)
library(magrittr)
library(ggplot2)
#
# #For each seurat object i use this function to store them into a list. So the orig.ident is the folder name
# ReadSeurat <- function(x){
#   seurat <- readRDS(sprintf("%s/seurat_%s.rds",x, x))
#   seurat$subID <- seurat$orig.ident
#   seurat$orig.ident <- sprintf('%s', x)
#   return(seurat)
# }
#
PlotGeneList <- read.delim2('Scripts_and_Data/Gene_Marker_list.txt',sep = ',')
seurat_combi <- readRDS('seurat_untouched_0718.rds')
# Seurat_List <- readRDS('seurat_list.rds')
# seurat_combi <- merge(x = Seurat_List[[1]], y = Seurat_List[-1], project = 'CombinedAll')

pdf(file="Integration_Visualization_0721.pdf",
    width = 24,
    height = 24)

Vplot1  <- VlnPlot(seurat_combi, features = c("nFeature_RNA"), pt.size = 0, group.by = 'orig.ident')
Vplot2  <- VlnPlot(seurat_combi, features = c("nCount_RNA"), pt.size = 0, group.by = 'orig.ident')
Vplot <- Vplot1 | Vplot2
print(Vplot)

seurat_combi <- seurat_combi %>%
    NormalizeData() %>%
    FindVariableFeatures(nfeatures = 3000) %>%
    ScaleData() %>%
    RunPCA(npcs = 50)
#
#
  pcs <- 35
# #
# # #First time visualization of the dataset in terms of different sources
plot1_Unint <- DimPlot(seurat_combi, group.by = "orig.ident") + plot_annotation(title="Visualization without integration")
plot2_Unint <- FeaturePlot(seurat_combi, features=colnames(PlotGeneList), ncol=6, pt.size = 0.5) + plot_annotation(title="Visualization without integration")

print(plot1_Unint)
print(plot2_Unint)
print("Finish plotting unintegrated")
#

library(harmony)
#should we use the group.by.vars??
seurat_combi <- RunHarmony(seurat_combi, group.by.vars = "orig.ident", dims.use = 1:pcs, max.iter.harmony = 50)
seurat_combi <- RunUMAP(seurat_combi, reduction = "harmony", dims = 1:pcs, reduction.name = 'har_umap',reduction.key="UMAPHAR_")
seurat_combi <- FindNeighbors(seurat_combi, reduction = "harmony", dims = 1:pcs, graph.name =  c('har_nn', 'har_snn')) %>%
  FindClusters(graph.name = 'har_snn', resolution = 0.6)
seurat_combi$har_clu <- seurat_combi$seurat_clusters
#
# #visualization of the results
plot1_Har <- DimPlot(seurat_combi, group.by="orig.ident", reduction = 'har_umap', pt.size = 0.5)
plot2_Har <- DimPlot(seurat_combi, label = T, reduction = 'har_umap', pt.size = 0.5)
plot3_Har <- FeaturePlot(seurat_combi, features=colnames(PlotGeneList),reduction = 'har_umap', ncol=6, pt.size = 0.5) + plot_annotation(title="harmony Integration")
#plot_Har <- plot1_Har + plot2_Har + plot_annotation(title = "harmony integration")
print(plot1_Har)
print(plot2_Har)
print(plot3_Har)
#
# print("Finish Harmony")

# ------------- CSS integration -------------
library(simspec)

#Should we use the label tag? Cuz for each dataset, there're multiple possible reasons for doing this
seurat_combi <- cluster_sim_spectrum(seurat_combi, label_tag = "orig.ident", cluster_resolution = 0.3)
#seurat_combi[['css']] <- CreateDimReducObject(Embeddings(seurat_combi, "css")[colnames(seurat_combi),], key="CSS_")
seurat_combi <- RunUMAP(seurat_combi, reduction="css", dims = 1:ncol(Embeddings(seurat_combi, "css")), reduction.name = 'css_umap',reduction.key="UMAPCSS_")
seurat_combi <- FindNeighbors(seurat_combi, reduction = "css", dims = 1:ncol(Embeddings(seurat_combi, "css")), graph.name = c('css_nn', 'css_snn')) %>%
  FindClusters(graph.name = 'css_snn', resolution = 0.6)
seurat_combi$css_clu <- seurat_combi$seurat_clusters

# #Visualization
# plot1_CSS <- DimPlot(seurat_combi, group.by="orig.ident", reduction = 'css_umap',pt.size = 0.5)
# plot2_CSS <- DimPlot(seurat_combi, label = T, reduction = 'css_umap', pt.size = 0.5)
# plot3_CSS <- FeaturePlot(seurat_combi, features=colnames(PlotGeneList), reduction = 'css_umap', ncol=6, pt.size = 0.5) + plot_annotation(title="CSS Integration")
# #plot_CSS <- plot1_CSS + plot2_CSS + plot_annotation(title = "CSS integration")
# print(plot1_CSS)
# print(plot2_CSS)
# print(plot3_CSS)
#dev.off()
print("Finish CSS")


# ------------- MNN integration -------------
library(SeuratWrappers)

seurat_samples <- SplitObject(seurat_combi, "orig.ident")
seurat_mnn <- RunFastMNN(seurat_samples)
seurat_combi[['mnn']] <- CreateDimReducObject(Embeddings(seurat_mnn, "mnn")[colnames(seurat_combi),], key="MNN_")
seurat_combi <- RunUMAP(seurat_combi, dims = 1:pcs, reduction = "mnn", reduction.name = 'mnn_umap',reduction.key="UMAPMNN_" )
seurat_combi <- FindNeighbors(seurat_combi, reduction = "mnn", dims = 1:pcs, graph.name = c('mnn_nn', 'mnn_snn') ) %>%
  FindClusters(graph.name = 'mnn_snn', resolution = 0.6)
seurat_combi$mnn_clu <- seurat_combi$seurat_clusters

#Visualization
plot1_MNN <- DimPlot(seurat_combi, group.by="orig.ident", reduction = 'mnn_umap', pt.size = 0.5)
plot2_MNN <- DimPlot(seurat_combi, label = T,reduction = 'mnn_umap',pt.size = 0.5)

#Change the width and height of pdf files.
plot3_MNN <- FeaturePlot(seurat_combi, reduction = 'mnn_umap', features=colnames(PlotGeneList), ncol=6, pt.size = 0.5) + plot_annotation(title="MNN Integration")
#plot_MNN <- plot1_MNN + plot2_MNN + plot_annotation(title = "MNN integration")
print(plot1_MNN)
print(plot2_MNN)
print(plot3_MNN)
print("Finish MNN")
#Scaronma Integration

dev.off()
saveRDS(seurat_combi, 'seurat_0718.rds')
