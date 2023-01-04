# ------------ Setting up the environment ------------
# devtools::install_github("cellgeni/sceasy")

library(Seurat)
library(Matrix)
library(patchwork)
library(dplyr)
library(magrittr)
library(ggplot2)
library(sceasy)


# #Set default into conda environment scvi-env
# # in terminal:
# #    conda activate scvi-env
# #    conda info --envs
# #    export LD_LIBRARY_PATH=/home/jiesun/packages/anaconda3/envs/scvi-env/lib:$LD_LIBRARY_PATH:
# 
# # Sys.setenv(RETICULATE_PYTHON = "~/packages/anaconda3/envs/scvi-env/bin/python")

# # py_config()


library(reticulate)
use_condaenv("scvi-env", required=TRUE, conda = '~/packages/anaconda3/condabin/conda')
# py_config()

sc <- import("scanpy", convert = FALSE)
scvi <- import("scvi", convert = FALSE)

# ------------ Convert Seurat Object ------------ 
seurat_combi <- readRDS('seurat_0721_filteredGenes.rds') %>%
  NormalizeData() %>%
  FindVariableFeatures(nfeatures = 2000) 

top2000 <- head(VariableFeatures(seurat_combi), 2000)
seurat_raw <- seurat_combi
seurat_combi <- seurat_combi[top2000]

# adata <- sc$AnnData(
#   X   = t(as.matrix(GetAssayData(seurat_combi,slot='counts'))), #scVI requires raw counts
#   obs = seurat_combi[[]],
#   var = GetAssay(seurat_combi)[[]]
# )

# Also try this: 
adata <- convertFormat(seurat_combi, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)
# identical(adata, adata1)

# ------------  ScVI model training ------------ 
# run seteup_anndata
scvi$model$SCVI$setup_anndata(adata)

# create the model
model = scvi$model$SCVI(adata,n_layers=2, n_latent=30, gene_likelihood="nb")

# train the model
model$train()

# to specify the number of epochs when training:
# model$train(n_epochs = as.integer(400))

# ------------ Convert back to Seurat and Clusteirng ------------ 
# get the latent represenation
latent = model$get_latent_representation()

# put it back in our original Seurat object
latent <- as.matrix(latent)
rownames(latent) = colnames(seurat_combi)
seurat_combi[['scvi']] <- CreateDimReducObject(embeddings = latent, key = "scvi_", assay = DefaultAssay(seurat_combi))

#Clustering
seurat_combi <- RunUMAP(seurat_combi, reduction="scvi", dims = 1:35, reduction.name = 'scvi_umap',reduction.key="UMAPscVI_")
seurat_combi <- FindNeighbors(seurat_combi, reduction = "scvi", dims = 1:35, graph.name =  c('scvi_nn', 'scvi_snn')) %>%
  FindClusters(graph.name = 'scvi_snn', resolution = 0.6)
seurat_combi$csvi_clu <- seurat_combi$seurat_clusters

# ------------ Visualization ------------ 
plot1_scvi <- DimPlot(seurat_combi, group.by="orig.ident", reduction = 'scvi_umap', pt.size = 0.5, raster=FALSE)
plot2_scvi <- DimPlot(seurat_combi, label = T,reduction = 'scvi_umap',pt.size = 0.5, raster=FALSE)
plot3_scvi <- FeaturePlot(seurat_combi, reduction = 'scvi_umap', features=colnames(PlotGeneList), ncol=6, pt.size = 0.5) + plot_annotation(title="scvi Integration")

pdf(file="Integration_scVI.pdf",
    width = 24,
    height = 24)

print(plot1_scvi)
print(plot2_scvi)
print(plot3_scvi)
print("Finish scvi")

dev.off()
saveRDS(seurat_combi, 'seurat_scvi.rds')
