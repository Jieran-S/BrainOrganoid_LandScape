library(SeuratWrappers)

seurat_samples <- SplitObject(seurat_combi, "orig.ident")
seurat_mnn <- RunFastMNN(seurat_samples)
seurat_combi[['mnn']] <- CreateDimReducObject(Embeddings(seurat_mnn, "mnn")[colnames(seurat_combi),], key="MNN_")
seurat_combi <- RunUMAP(seurat_combi, dims = 1:pcs, reduction = "mnn", reduction.name = 'mnn_umap',reduction.key="UMAPMNN_" )
seurat_combi <- FindNeighbors(seurat_combi, reduction = "mnn", dims = 1:pcs, graph.name = 'mnn_snn' ) %>%
  FindClusters(graph.name = 'mnn_snn' , resolution = 0.6)
seurat_combi$mnn_clu <- seurat_combi$seurat_clusters

#Visualization
plot1_MNN <- UMAPPlot(seurat_combi, group.by="orig.ident", reduction = 'mnn_umap')
plot2_MNN <- UMAPPlot(seurat_combi, label = T,reduction = 'mnn_umap')

#Change the width and height of pdf files.
plot3_MNN <- FeaturePlot(seurat_combi, reduction = 'mnn_umap', features=colnames(PlotGeneList), ncol=6) + plot_annotation(title="MNN Integration")
#plot_MNN <- plot1_MNN + plot2_MNN + plot_annotation(title = "MNN integration")
print(plot1_MNN)
print(plot2_MNN)
print(plot3_MNN)
print("Finish MNN")
