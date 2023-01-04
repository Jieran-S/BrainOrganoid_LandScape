# Added the filter for each seurat object to leave only the intersectional genes before merging
# removal of the artificially introduced batch effects

# ------------- Environment setup -------------
library(Seurat)
library(Matrix)
library(patchwork)
library(dplyr)
library(magrittr)
library(ggplot2)

PlotGeneList <- read.delim2('Scripts_and_Data/Gene_Marker_list.txt',sep = ',')

# ------------- Seurat Object Handling (Batch effect removal) -------------
seurat.list <- readRDS('seurat_list.rds')

# Taking the intersect of the features while keeping the annotation genes attached
list.genes <- lapply(seurat.list, function(x){rownames(x)})
intersect.genes <- Reduce(intersect, list.genes) %>% 
  union(colnames(PlotGeneList))

# list.intersect <- list()
# for (i in 1:length(list.genes)){
#   list.intersect[[i]] <- Reduce(intersect, list.genes[-i]) 
# }
# intersect.16 <- unique(unlist(list.intersect)) %>% union(colnames(PlotGeneList))


# Trimming seurat objects 
seurat.list <- lapply(seurat.list, function(x){x[intersect.genes, ]})
seurat_combi <- merge(x = seurat.list[[1]], y = seurat.list[-1], project = 'CombinedAll')
saveRDS(seurat_combi, 'seurat_0721_filteredGenes.rds')

# ------------- General Protocol before Integration -------------
pdf(file="Visualization/Integration_FilterGenes.pdf",
    width = 24,
    height = 24)

seurat_combi <- seurat_combi %>%
  NormalizeData() %>%
  FindVariableFeatures(nfeatures = 3000) %>%
  ScaleData() %>%
  RunPCA(npcs = 50)

Vplot1  <- VlnPlot(seurat_combi, features = c("nFeature_RNA"), pt.size = 0, group.by = 'orig.ident')
Vplot2  <- VlnPlot(seurat_combi, features = c("nCount_RNA"), pt.size = 0, group.by = 'orig.ident')
Vplot <- Vplot1 | Vplot2
print(Vplot)

#Identify no. of dimension to use for integration
pct <- seurat_combi[["pca"]]@stdev / sum(seurat_combi[["pca"]]@stdev) * 100
choice1 <- which(cumsum(pct) > 80 & pct < 5)[1]
choice2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1 & cumsum(pct) > 60 ), decreasing = T)[1] + 1
pcs <- min(choice1, choice2, 40)
print(sprintf('pct: %s', pcs))

seurat_combi <- RunUMAP(seurat_combi, dims= 1:pcs)
plot1_Unint <- DimPlot(seurat_combi, reduction = 'umap', group.by = "orig.ident", raster=FALSE) + plot_annotation(title="Visualization without integration")
plot2_Unint <- FeaturePlot(seurat_combi, features=colnames(PlotGeneList), ncol=6, pt.size = 0.5) + plot_annotation(title="Visualization without integration")

print(plot1_Unint)
print(plot2_Unint)
print("Finish plotting unintegrated")

# ------------- Harmony Integration -------------
library(harmony)

seurat_combi <- RunHarmony(seurat_combi, group.by.vars = "orig.ident", dims.use = 1:pcs, max.iter.harmony = 50)
seurat_combi <- RunUMAP(seurat_combi, reduction = "harmony", dims = 1:pcs, reduction.name = 'har_umap',reduction.key="UMAPHAR_")
seurat_combi <- FindNeighbors(seurat_combi, reduction = "harmony", dims = 1:pcs, graph.name =  c('har_nn', 'har_snn')) %>%
  FindClusters(graph.name = 'har_snn', resolution = 0.6)
seurat_combi$har_clu <- seurat_combi$seurat_clusters

# #visualization of the results
plot1_Har <- DimPlot(seurat_combi, group.by="orig.ident", reduction = 'har_umap', pt.size = 0.5, raster=FALSE)+ plot_annotation(title="harmony Integration")
plot2_Har <- DimPlot(seurat_combi, label = T, reduction = 'har_umap', pt.size = 0.5, raster=FALSE)+ plot_annotation(title="harmony Integration")
plot3_Har <- FeaturePlot(seurat_combi, features=colnames(PlotGeneList),reduction = 'har_umap', ncol=6, pt.size = 0.5) + plot_annotation(title="harmony Integration")

print(plot1_Har)
print(plot2_Har)
print(plot3_Har)

print("Finish Harmony")

# ------------- CSS integration -------------

library(simspec)

seurat_combi <- cluster_sim_spectrum(seurat_combi, label_tag = "orig.ident", cluster_resolution = 0.3)
seurat_combi <- RunUMAP(seurat_combi, reduction="css", dims = 1:ncol(Embeddings(seurat_combi, "css")), reduction.name = 'css_umap',reduction.key="UMAPCSS_")
seurat_combi <- FindNeighbors(seurat_combi, reduction = "css", dims = 1:ncol(Embeddings(seurat_combi, "css")), graph.name = c('css_nn', 'css_snn')) %>%
  FindClusters(graph.name = 'css_snn', resolution = 0.6)
seurat_combi$css_clu <- seurat_combi$seurat_clusters

#Visualization
plot1_CSS <- DimPlot(seurat_combi, group.by="orig.ident", reduction = 'css_umap',pt.size = 0.5,  raster=FALSE)+ plot_annotation(title="CSS Integration")
plot2_CSS <- DimPlot(seurat_combi, label = T, reduction = 'css_umap', pt.size = 0.5,  raster=FALSE) + plot_annotation(title="CSS Integration")
plot3_CSS <- FeaturePlot(seurat_combi, features=colnames(PlotGeneList), reduction = 'css_umap', ncol=6, pt.size = 0.5) + plot_annotation(title="CSS Integration")

print(plot1_CSS)
print(plot2_CSS)
print(plot3_CSS)

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
plot1_MNN <- DimPlot(seurat_combi, group.by="orig.ident", reduction = 'mnn_umap', pt.size = 0.5, raster=FALSE) + plot_annotation(title="MNN Integration")
plot2_MNN <- DimPlot(seurat_combi, label = T,reduction = 'mnn_umap',pt.size = 0.5, raster=FALSE) + plot_annotation(title="MNN Integration")
plot3_MNN <- FeaturePlot(seurat_combi, reduction = 'mnn_umap', features=colnames(PlotGeneList), ncol=6, pt.size = 0.5) + plot_annotation(title="MNN Integration")

print(plot1_MNN)
print(plot2_MNN)
print(plot3_MNN)
print("Finish MNN")

dev.off()
saveRDS(seurat_combi, 'seurat_0721_RemoveBatch.rds')
