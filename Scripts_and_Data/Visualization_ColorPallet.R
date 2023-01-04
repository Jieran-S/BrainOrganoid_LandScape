library(Seurat)
library(Matrix)
library(patchwork)
library(dplyr)
library(magrittr)
library(ggplot2)
library(RColorBrewer)

seurat_combi <- readRDS('seurat_0721_int_RemoveBatch.rds')

integraiton.list <- c('har', 'css', 'mnn')

pdf(file="Visualization/Integration_ColorPallet.pdf",
    width = 24,
    height = 24)

colorVec <- DiscretePalette(n, palette = 'set3')

for (name in integration.list){
  plot1 <- DimPlot(seurat_combi, group.by="orig.ident", reduction = sprintf('%s_umap', name), pt.size = 0.5, raster=FALSE) + plot_annotation(title=sprintf('%s integration', name))
  plot2 <- DimPlot(seurat_combi, label = T,reduction = sprintf('%s_umap', name),pt.size = 0.5, raster=FALSE) + plot_annotation(title=sprintf('%s integration', name))
  print(plot1)
  pirnt(plot2)
  }

