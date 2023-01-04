# A script for integration based on previous trials. Trying CSS, MNN, Harmony and SEVI

# A further improvement may be based on the https://doi.org/10.1038/s41592-022-01480-9 
# to evaluate the integration efficiency and check out the possible best method for the 
# integration

# ------------- Environment ------------- 
library(Seurat)
library(Matrix)
library(patchwork)
library(dplyr)
library(magrittr)
library(ggplot2)
ReadSeurat <- function(x){
  seurat <- readRDS(sprintf("%s/seurat_%s.rds",x, x))
  seurat$subID <- seurat$orig.ident
  seurat$orig.ident <- sprintf('%s', x)
  return(seurat)
}

# ------------- import and merge all seurat objects (Only need to run once) ------------ 
#remove <- c("2019_Yoon_BD_New", "2021_Bowles_Drop_New", "Scripts_and_Data", 'Integration_Visualization_All_3.pdf', 'seurat_all_unint.rds', 'seurat_all.rds', 'Integration_Visualization_All_4.pdf', "scPipeline")
#Seurat_List <- lapply(setdiff(list.files(), remove), ReadSeurat)
#seurat_combi <- merge(x = Seurat_List[[1]], y = Seurat_List[-1], project = 'CombinedAll')
PlotGeneList <- read.delim2('Scripts_and_Data/Gene_Marker_list.txt',sep = ',')

#General processing for integration
#seurat_combi <- seurat_combi %>%
#  NormalizeData() %>%
#  FindVariableFeatures(nfeatures = 3000) %>%
#  ScaleData() %>%
#  RunPCA(npcs = 50)

#saveRDS(Seurat_List, 'all_seurat_list.rds')
#saveRDS(seurat_combi, 'seurat_all.rds')

# ------------- Start integration ------------
#seurat_combi <- readRDS('seurat_all.rds')
Seurat_List <- readRDS('seurat_list.rds') 
seurat_combi <- merge(x = Seurat_List[[1]], y = Seurat_List[-1], project = 'CombinedAll')

pdf(file="Integration_Visualization_0708.pdf", 
    width = 24,
    height = 24)

#Violin plot without integration
Vplot1  <- VlnPlot(seurat_combi, features = c("nFeature_RNA"), pt.size = 0, group.by = 'orig.ident')
Vplot2  <- VlnPlot(seurat_combi, features = c("nCount_RNA"), pt.size = 0, group.by = 'orig.ident')
Vplot <- Vplot1 | Vplot2
print(Vplot)

pct <- seurat_combi[["pca"]]@stdev / sum(seurat_combi[["pca"]]@stdev) * 100
#seurat_combi[['pca']] <- CreateDimReducObject(Embeddings(seurat_combi, "pca")[colnames(seurat_combi),], key="PCA_")
choice1 <- which(cumsum(pct) > 80 & pct < 5)[1]
choice2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1 & cumsum(pct) > 60 ), decreasing = T)[1] + 1
pcs <- min(choice1, choice2, 40)

seurat_combi <- RunUMAP(seurat_combi, dims= 1:pcs)

plot_elbow <- ElbowPlot(seurat_combi, ndims= 50) + 
  geom_point(x = pcs, y = seurat_combi[["pca"]]@stdev[pcs] , colour = "red") +
  geom_label(
    label=sprintf("PC selected: %1.0f",pcs), 
    x=pcs,
    y=10,
    label.padding = unit(0.55, "lines"), # Rectangle size around label
    label.size = 1,
    color = "black",
    fill="#69b3a2"
  )
print(plot_elbow)

#First time visualization of the dataset in terms of different sources
plot1_Unint <- DimPlot(seurat_combi, reduction = 'umap', group.by = "orig.ident") + plot_annotation(title="Visualization without integration")
plot2_Unint <- FeaturePlot(seurat_combi, features=colnames(PlotGeneList), ncol=6, pt.size = 0.5) + plot_annotation(title="Visualization without integration")

#plot3_Unint <- plot1_Unint + plot2_Unint + plot_layout(widths = c(1.5, 2)) + plot_annotation(title="Visualization without integration")
print(plot1_Unint)
print(plot2_Unint)
print("Finish plotting unintegrated")

# ------------- harmony integration  ------------- 
library(harmony)

seurat_combi <- RunHarmony(seurat_combi, group.by.vars = "orig.ident", dims.use = 1:pcs, max.iter.harmony = 50)
seurat_combi <- RunUMAP(seurat_combi, reduction = "harmony", dims = 1:pcs, reduction.name = 'har_umap',reduction.key="UMAPHAR_")
seurat_combi <- FindNeighbors(seurat_combi, reduction = "harmony", dims = 1:pcs, graph.name = 'har_snn') %>% 
  FindClusters(graph.name = 'har_snn', resolution = 0.6)
seurat_combi$har_clu <- seurat_combi$seurat_clusters

#visualization of the results
plot1_Har <- UMAPPlot(seurat_combi, group.by="orig.ident", reduction = 'har_umap', pt.size = 0.5)
plot2_Har <- UMAPPlot(seurat_combi, label = T, reduction = 'har_umap', pt.size = 0.5)
plot3_Har <- FeaturePlot(seurat_combi, features=colnames(PlotGeneList),reduction = 'har_umap', ncol=6, pt.size = 0.5) + plot_annotation(title="harmony Integration")
#plot_Har <- plot1_Har + plot2_Har + plot_annotation(title = "harmony integration")
print(plot1_Har)
print(plot2_Har)
print(plot3_Har)

print("Finish Harmony")

# ------------- CSS integration ------------- 
library(simspec)

seurat_combi <- cluster_sim_spectrum(seurat_combi, label_tag = "orig.ident", cluster_resolution = 0.3)
plot1_Har <- UMAPPlot(seurat_combi, group.by="orig.ident", pt.size = 0.5)
#seurat_combi[['css']] <- CreateDimReducObject(Embeddings(seurat_combi, "css")[colnames(seurat_combi),], key="CSS_")
seurat_combi <- RunUMAP(seurat_combi, reduction="css", dims = 1:ncol(Embeddings(seurat_combi, "css")), reduction.name = 'css_umap',reduction.key="UMAPCSS_")
seurat_combi <- FindNeighbors(seurat_combi, reduction = "css", dims = 1:ncol(Embeddings(seurat_combi, "css")), graph.name = 'css_snn') %>%
  FindClusters(graph.name = 'css_snn', resolution = 0.6)
seurat_combi$css_clu <- seurat_combi$seurat_clusters

#Visualization
plot1_CSS <- UMAPPlot(seurat_combi, group.by="orig.ident", reduction = 'css_umap',pt.size = 0.5)
plot2_CSS <- UMAPPlot(seurat_combi, label = T, reduction = 'css_umap', pt.size = 0.5)
plot3_CSS <- FeaturePlot(seurat_combi, features=colnames(PlotGeneList), reduction = 'css_umap', ncol=6, pt.size = 0.5) + plot_annotation(title="CSS Integration")
#plot_CSS <- plot1_CSS + plot2_CSS + plot_annotation(title = "CSS integration")
print(plot1_CSS)
print(plot2_CSS)
print(plot3_CSS)
#dev.off()
print("Finish CSS")


# ------------- MNN integration -------------
library(SeuratWrappers)

seurat_samples <- SplitObject(seurat_combi, "orig.ident")
seurat_mnn <- RunFastMNN(seurat_samples)
seurat_combi[['mnn']] <- CreateDimReducObject(Embeddings(seurat_mnn, "mnn")[colnames(seurat_combi),], key="MNN_")
seurat_combi <- RunUMAP(seurat_combi, dims = 1:pcs, reduction = "mnn", reduction.name = 'mnn_umap',reduction.key="UMAPMNN_" )
seurat_combi <- FindNeighbors(seurat_combi, reduction = "mnn", dims = 1:pcs, graph.name = 'mnn_snn' ) %>%
  FindClusters(graph.name = 'mnn_snn', resolution = 0.6)
seurat_combi$mnn_clu <- seurat_combi$seurat_clusters

#Visualization
plot1_MNN <- UMAPPlot(seurat_combi, group.by="orig.ident", reduction = 'mnn_umap', pt.size = 0.5)
plot2_MNN <- UMAPPlot(seurat_combi, label = T,reduction = 'mnn_umap',pt.size = 0.5)

#Change the width and height of pdf files.
plot3_MNN <- FeaturePlot(seurat_combi, reduction = 'mnn_umap', features=colnames(PlotGeneList), ncol=6, pt.size = 0.5) + plot_annotation(title="MNN Integration")
#plot_MNN <- plot1_MNN + plot2_MNN + plot_annotation(title = "MNN integration")
print(plot1_MNN)
print(plot2_MNN)
print(plot3_MNN)
print("Finish MNN")

#  ------------- CCA integration  ------------- 
#We identify anchors using the FindIntegrationAnchors function, which takes a list of Seurat objects as input.
seurat_integrate <- Seurat_List %>%
  FindIntegrationAnchors(dims = 1:50, reduction = 'cca') %>%
  IntegrateData(dims = 1:50, new.assay.name = 'CCA')

#transfer new assay into the old seurat subject
seurat_combi[['CCA']] <- seurat_integrate[['CCA']]
DefaultAssay(seurat_combi) <- 'CCA'

#Normal pipeline conduction: Before trying this maybe add a copy just in case 
seruat_combi <- seurat_combi %>% 
  NormalizeData() %>%
  ScaleData() %>%
  RunPCA(npcs = 50) %>%

pct <- seurat_combi[["pca"]]@stdev / sum(seurat_combi[["pca"]]@stdev) * 100
choice1 <- which(cumsum(pct) > 80 & pct < 5)[1]
choice2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1 & cumsum(pct) > 60 ), decreasing = T)[1] + 1
pcs <- min(choice1, choice2, 40)

seruat_combi <- seurat_combi %>% 
  RunUMAP(dims = 1:pcs,reduction.name = 'CCA_umap', reduction.key = 'UMAPCCA_') %>%
  FindNeighbors(dims = 1:pcs, graph.name = 'CCA_snn') %>%
  FindClusters(graph.name = 'CCA_snn', resolution = 0.6)

#Visualization
plot1_CCA <- UMAPPlot(seurat_combi, group.by="orig.ident", reduction = 'CCA_umap')
plot2_CCA <- UMAPPlot(seurat_combi, label = T,reduction = 'CCA_umap')
plot3_CCA <- FeaturePlot(seurat_combi, reduction = 'CCA_umap', features=colnames(PlotGeneList), ncol=6) + 
  plot_annotation(title="CCA Integration")
seurat_combi$cca_clu <- seurat_combi$seurat_clusters
DefaultAssay(seurat_combi) <- 'RNA'

print(plot1_CCA)
print(plot2_CCA)
print(plot3_CCA)
print("Finish CCA")

# ------------- scanorama integration ------------- 
#Set default into conda environment
# in terminal: export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/jiesun/packages/anaconda3/lib/
Sys.setenv(RETICULATE_PYTHON = "~/packages/anaconda3/bin/python")
library(reticulate)
scanorama <- import("scanorama")

hvgs <- lapply(Seurat_List, function(x){x@assays$RNA@var.features}) %>%
  unlist() %>%
  unique()
#hvgs <- unique(unlist(hvgs_per_dataset))

#Set up cell matrix and gene lists
assaylist <- list()
genelist <- list()
for (i in 1:length(alldata.list)) {
  assaylist[[i]] <- t(as.matrix(GetAssayData(Seurat_List[[i]], "data")[hvgs, ]))
  genelist[[i]] <- hvgs
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
seurat_combi <- FindNeighbors(seurat_combi, dims = 1:scans, reduction = "scan", graph.name = 'scan_snn' ) %>%
  FindClusters(graph.name= 'scan_snn', resolution = 0.6)
seurat_combi$scan_clu <- seurat_combi$seurat_clusters

#Visualization
plot1_scan <- UMAPPlot(seurat_combi, group.by="orig.ident", reduction = 'scan_umap')
plot2_scan <- UMAPPlot(seurat_combi, label = T,reduction = 'scan_umap')

#Change the width and height of pdf files.
plot3_scan <- FeaturePlot(seurat_combi, reduction = 'scan_umap', features=colnames(PlotGeneList), ncol=6) + plot_annotation(title="scanorama Integration")

print(plot1_scan)
print(plot2_scan)
print(plot3_scan)
print("Finish scanorama")

# ------------- Finishing up -------------  
dev.off()
saveRDS(seurat_combi,'seurat_integrated.rds')

