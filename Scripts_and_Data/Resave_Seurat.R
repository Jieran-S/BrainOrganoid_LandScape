#Save all the previous data that hasn't been saved

library(Seurat)
library(patchwork)
library(dplyr)
library(magrittr)
library(Matrix)

remove = c("2019_Yoon_BD_New", "2021_Bowles_Drop_New")

for (file in setdiff(list.files(), remove)){
  seurat_combi <- readRDS(sprintf("%s/seurat_%s.rds",file, file)) 
  
  if(any(names(seurat_combi)== 'mnn')){

    seurat_combi <- seurat_combi %>%
      NormalizeData() %>%
      FindVariableFeatures(nfeatures = 3000) %>%
      ScaleData() %>%
      RunPCA(npcs = 50)
    
    pct <- seurat_combi[["pca"]]@stdev / sum(seurat_combi[["pca"]]@stdev) * 100
    choice1 <- which(cumsum(pct) > 80 & pct < 5)[1]
    choice2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1 & cumsum(pct) > 60 ), decreasing = T)[1] + 1
    pcs <- min(choice1, choice2, 40)
    
    library(harmony)
    seurat_combi <- RunHarmony(seurat_combi, group.by.vars = "orig.ident", dims.use = 1:pcs, max.iter.harmony = 50)
    
    library(simspec)
    seurat_combi <- cluster_sim_spectrum(seurat_combi, label_tag = "orig.ident", cluster_resolution = 0.3)
    
    library(SeuratWrappers)
    seurat_samples <- SplitObject(seurat_combi, "orig.ident")
    seurat_mnn <- RunFastMNN(seurat_samples)
    seurat_combi[['mnn']] <- CreateDimReducObject(Embeddings(seurat_mnn, "mnn")[colnames(seurat_combi),], key="MNN_")
    seurat_combi <- RunUMAP(seurat_combi, dims = 1:pcs, reduction = "mnn", reduction.name = 'mnn_umap',reduction.key="UMAPMNN_" )
    seurat_combi <- FindNeighbors(seurat_combi, reduction = "mnn", dims = 1:pcs, graph.name = 'mnn_snn' ) %>%
      FindClusters(resolution = 0.6)
    
  }else{
    
    pct <- seurat_combi[["pca"]]@stdev / sum(seurat_combi[["pca"]]@stdev) * 100
    choice1 <- which(cumsum(pct) > 80 & pct < 5)[1]
    choice2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1 & cumsum(pct) > 60 ), decreasing = T)[1] + 1
    pcs <- min(choice1, choice2, 40)
  }
  
  seurat_combi <- RunUMAP(seurat_combi, dims= 1:pcs)
  seurat_combi <- RunUMAP(seurat_combi, reduction = "harmony", dims = 1:pcs, reduction.name = 'har_umap',reduction.key="UMAPHAR_")
  seurat_combi <- FindNeighbors(seurat_combi, reduction = "harmony", dims = 1:pcs, graph.name = 'har_snn') %>% FindClusters(resolution = 0.6)
  seurat_combi <- RunUMAP(seurat_combi, reduction="css", dims = 1:ncol(Embeddings(seurat_combi, "css")), reduction.name = 'css_umap',reduction.key="UMAPCSS_")
  seurat_combi <- FindNeighbors(seurat_combi, reduction = "css", dims = 1:ncol(Embeddings(seurat_combi, "css")), graph.name = 'css_snn') %>%
    FindClusters(resolution = 0.6)
  
  saveRDS(seurat_combi, sprintf("%s/seurat_%s.rds",file, file))
  print(sprintf("%s finished",file))
}

