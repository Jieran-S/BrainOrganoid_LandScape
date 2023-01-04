library(Seurat)
library(Matrix)
library(patchwork)
library(dplyr)
library(magrittr)
library(ggplot2)

# 
pdf(file="Visualization/Integration_CCA_0721_filteredGene.pdf",
         width = 24,
         height = 24)
# 
# # ------------- import and merge all seurat objects (Only need to run once) ------------ 
#remove <- c("2019_Yoon_BD_New", "2021_Bowles_Drop_New", "Scripts_and_Data", 'Integration_Visualization_All_3.pdf', 'seurat_all_unint.rds', 'seurat_all.rds', 'Integration_Visualization_All_4.pdf', "scPipeline")
#Seurat_List <- lapply(setdiff(list.files(), remove), ReadSeurat)
#seurat_combi <- merge(x = Seurat_List[[1]], y = Seurat_List[-1], project = 'CombinedAll')
PlotGeneList <- read.delim2('Scripts_and_Data/Gene_Marker_list.txt',sep = ',')

#  ------------- CCA integration  -------------
#We identify anchors using the FindIntegrationAnchors function, which takes a list of Seurat objects as input.s
seurat_combi <- readRDS('seurat_0721_filteredGenes.rds')
Seurat_List <- readRDS('seurat_list.rds') %>%
  lapply(function(x){
    x[rownames(seurat_combi),] %>%
      NormalizeData() %>%
      FindVariableFeatures(nfeatures = 3000) %>%
      ScaleData() %>%
      RunPCA(npcs = 50)
  })

integ_features <- SelectIntegrationFeatures(object.list = Seurat_List, 
                                            nfeatures = 3000) 
split_seurat <- PrepSCTIntegration(object.list = Seurat_List, 
                                   anchor.features = integ_features)

seurat_integrate <- Seurat_List %>%
  FindIntegrationAnchors(dims = 1:50, reduction = 'cca', anchor.features = integ_features, normalization.method = "SCT") %>%
  IntegrateData(dims = 1:50, new.assay.name = 'CCA', normalization.method = "SCT")

#transfer new assay into the old seurat subject
seurat_combi[['CCA']] <- seurat_integrate[['CCA']]
DefaultAssay(seurat_combi) <- 'CCA'

#Normal pipeline conduction: Before trying this maybe add a copy just in case 
seurat_combi <- seurat_combi %>%
 ScaleData() %>%
 RunPCA(npcs = 50)
  
pct <- seurat_combi[["pca"]]@stdev / sum(seurat_combi[["pca"]]@stdev) * 100
choice1 <- which(cumsum(pct) > 80 & pct < 5)[1]
choice2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1 & cumsum(pct) > 60 ), decreasing = T)[1] + 1
pcs <- min(choice1, choice2, 40)

seurat_combi <- seurat_combi %>% 
  RunUMAP(dims = 1:pcs,reduction.name = 'CCA_umap', reduction.key = 'UMAPCCA_') %>%
  FindNeighbors(dims = 1:pcs, graph.name = c('CCA_nn', 'CCA_snn')) %>%
  FindClusters(graph.name = 'CCA_snn', resolution = 0.6)

# -------------  Visualization -------------
plot1_CCA <- DimPlot(seurat_combi, group.by="orig.ident", reduction = 'CCA_umap')+ 
  plot_annotation(title="CCA Integration")
plot2_CCA <- DimPlot(seurat_combi, label = T,reduction = 'CCA_umap')+ 
  plot_annotation(title="CCA Integration")
plot3_CCA <- FeaturePlot(seurat_combi, reduction = 'CCA_umap', features=colnames(PlotGeneList), ncol=6) + 
  plot_annotation(title="CCA Integration")
seurat_combi$cca_clu <- seurat_combi$seurat_clusters
DefaultAssay(seurat_combi) <- 'RNA'
plot4_CCA <- FeaturePlot(seurat_combi, reduction = 'CCA_umap', features=colnames(PlotGeneList), ncol=6) + 
  plot_annotation(title="CCA Integration")

print(plot1_CCA)
print(plot2_CCA)
print(plot3_CCA)
print(plot4_CCA)
print("Finish CCA")

dev.off()
saveRDS(seurat_combi, 'seurat_CCA.rds')