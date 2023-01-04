#Set default into conda environment
# in terminal: export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/jiesun/packages/anaconda3/lib/
library(Seurat)
library(Matrix)
library(patchwork)
library(dplyr)
library(magrittr)
library(ggplot2)

Sys.setenv(RETICULATE_PYTHON = "~/packages/anaconda3/bin/python")
library(reticulate)
scanorama <- import("scanorama")

seurat_combi <- readRDS('seurat_0721_filteredGenes.rds') 

Seurat_List <- readRDS('seurat_list.rds') %>%
    lapply(function(x){
      x[rownames(seurat_combi),] %>%
      NormalizeData() %>%
      FindVariableFeatures(nfeatures = 3000) %>%
      ScaleData() %>%
      RunPCA(npcs = 50)
    })

pdf(file="Integration_Visualization_Scan.pdf", 
    width = 24,
    height = 24)

#Some issue here: What kind of genes are we looking for in the integration? The selection
# should be based on all data or just a few data in individual datasets?
hvgs <- lapply(Seurat_List, function(x){x@assays$RNA@var.features[1:3000]})%>%
  unlist()

#???? hvg.selected <- unique(hvgs[which(table(hvgs)>1)])

#hvgs <- lapply(seq(hvgs1),function(i) Reduce(intersect,lapply(hvgs1,"[[",i)))

#Set up cell matrix and gene lists
assaylist <- list()
genelist <- list()
for (i in 1:length(Seurat_List)) {
  #print(i)
  assaylist[[i]] <- t(as.matrix( GetAssayData(Seurat_List[[i]], "data")[Reduce(intersect, list(hvgs, rownames(Seurat_List[[i]]))), ] ))
  genelist[[i]] <- Reduce(intersect, list(hvgs, rownames(Seurat_List[[i]])))
}

#Integration
integrated.data <- scanorama$integrate(datasets_full = assaylist, genes_list = genelist)
intdimred <- do.call(rbind, integrated.data[[1]])
colnames(intdimred) <- paste0("PC_", 1:100)
rownames(intdimred) <- colnames(seurat_combi)
stdevs <- apply(intdimred, MARGIN = 2, FUN = sd)

#Introducing the reduction objects into the system
seurat_combi[["scan"]] <- CreateDimReducObject(embeddings = intdimred, stdev = stdevs,
                                               key = "scan_", assay = "RNA")
#New elbow plot to choose how much to use for UMAP 
scant <- seurat_combi[["scan"]]@stdev / sum(seurat_combi[["scan"]]@stdev) * 100
choice1 <- which(cumsum(scant) > 80 & scant < 5)[1]
choice2 <- sort(which((scant[1:length(scant) - 1] - scant[2:length(scant)]) > 0.1 & cumsum(scant) > 60 ), decreasing = T)[1] + 1
scans <- min(choice1, choice2, 80)
if (is.na(scans)){scans <- 80}

seurat_combi <- RunUMAP(seurat_combi, dims = 1:scans, reduction = "scan", reduction.name = 'scan_umap',reduction.key="UMAPscan_" )
seurat_combi <- FindNeighbors(seurat_combi, dims = 1:scans, reduction = "scan", graph.name =  c('scan_nn', 'scan_snn') ) %>%
  FindClusters(graph.name= 'scan_snn', resolution = 0.6)
seurat_combi$scan_clu <- seurat_combi$seurat_clusters

#Visualization
plot1_scan <- DimPlot(seurat_combi, group.by="orig.ident", reduction = 'scan_umap', raster = F)
plot2_scan <- DimPlot(seurat_combi, label = T,reduction = 'scan_umap', raster = F)

#Change the width and height of pdf files.
plot3_scan <- FeaturePlot(seurat_combi, reduction = 'scan_umap', features=colnames(PlotGeneList), ncol=6, raster = F) + plot_annotation(title="scanorama Integration")

print(plot1_scan)
print(plot2_scan)
print(plot3_scan)
print("Finish scanorama")

# ------------- Finishing up -------------  
dev.off()
saveRDS(seurat_combi,'seurat_integrated_scan.rds')