# ------------- Environment setup -------------
library(Seurat)
library(Matrix)
library(patchwork)
library(dplyr)
library(magrittr)

# for plotting and color tuning
library(ggplot2)
library(scales)
library(RColorBrewer)


# ------------- Plotting attempts for seurat integraiton without information -------------

seurat.plot <- readRDS('seurat_0721_int_RemoveBatch.rds')
PlotGeneList <- read.delim2('Scripts_and_Data/Gene_Marker_list.txt',sep = ',')

pdf(file="Visualization/Integration_ColorToned.pdf",
    width = 24,
    height = 24)

# ------------- Assigning new sequencing categories into seurat objects  -------------
Seq_tech <- c(
  "2017_Birey_BD"             = "BD Resolve system",             
  "2017_Sloan_NextSeq_New"    = "Smart-Seq2",
  "2019_Kanton_10x"           = "10x",        
  "2019_Marton_SmartSeq2_New" = "Smart-Seq2",
  "2019_Trujillo_10x"         = "10x",            
  "2019_Velasco_10x"          = "10x",
  "2019_Xiang_10x"            = "10x",
  "2020_Andersen_10x_New"     = "10x",
  "2020_Miura_10x_New"        = "10x",
  "2020_Qian_SPLiT_New"       = "SPLit-Seq",
  "2020_Sawada_BD"            = "BD Resolve system",             
  "2020_Smit_Drop"            = 'Drop-seq',
  "2020_Themasap_10x_New"     = "10x",
  "2021_Agoglia_SmartSeq2_New"= "Smart-Seq2",
  "2021_Fiorenzano_10x"       = "10x", 
  "2021_Huang_Drop"           = 'Drop-seq',
  "2022_Kelava_10x_New"       = "10x"
)

Cultech_List <- c(
  "2017_Birey_BD"             = "Pasca, et al., 2015",             
  "2017_Sloan_NextSeq_New"    = "Pasca, et al., 2015",
  "2019_Kanton_10x"           = "Pasca, et al., 2015",        
  "2019_Marton_SmartSeq2_New" = "Pasca, et al., 2015",
  "2019_Trujillo_10x"         = "Pasca, et al., 2015",            
  "2019_Velasco_10x"          = "Lancaster, et al., 2013",
  "2019_Xiang_10x"            = "Pasca, et al., 2015",
  "2020_Andersen_10x_New"     = "Pasca, et al., 2015",
  "2020_Miura_10x_New"        = "Pasca, et al., 2015",
  "2020_Qian_SPLiT_New"       = "Lancaster, et al., 2013",
  "2020_Sawada_BD"            = "Lancaster, et al., 2013",             
  "2020_Smit_Drop"            = 'Lancaster, et al., 2013',
  "2020_Themasap_10x_New"     = "Pasca, et al., 2015",
  "2021_Agoglia_SmartSeq2_New"= "Pasca, et al., 2015",
  "2021_Fiorenzano_10x"       = "Lancaster, et al., 2013", 
  "2021_Huang_Drop"           = 'Lancaster, et al., 2013',
  "2022_Kelava_10x_New"       = "Lancaster, et al., 2013"
)

seurat.plot$Culture.Protocol <- Cultech_List[seurat.plot$orig.ident]
seurat.plot$Sequence.Technology <- Seq_tech[seurat.plot$orig.ident]

# ------------- Changing color tone for visualization -------------

cols <- rev(brewer.pal(8, 'RdBu'))
color.cul <- c("#8DD3C7", "#FFA68E")
color.tech <- c("#7bcec6","#D87F81","#2aadf3","#fcd200","#2eae8f")
color.features <- brewer.pal(8, 'YlGnBu')
# Harmony: 49 clusters, CSS: 30, MNN: 28
color.cluster <- c("#6bd6db","#486492","#ff9b6b","#fff814","#f44655")
pal <- colorRampPalette(color.cluster)

# ------------- plotting graphs -------------
seurat_combi <- seurat.plot
Vplot1  <- VlnPlot(seurat_combi, features = c("nFeature_RNA"), pt.size = 0, group.by = 'Culture.Protocol')
Vplot2  <- VlnPlot(seurat_combi, features = c("nCount_RNA"), pt.size = 0, group.by = 'Culture.Protocol')
Vplot <- Vplot1 | Vplot2
print(Vplot)

Vplot1  <- VlnPlot(seurat_combi, features = c("nFeature_RNA"), pt.size = 0, group.by = 'Sequence.Technology')
Vplot2  <- VlnPlot(seurat_combi, features = c("nCount_RNA"), pt.size = 0, group.by = 'Sequence.Technology')
Vplot <- Vplot1 | Vplot2
print(Vplot)


plot1_Unint <- DimPlot(seurat_combi, reduction = 'umap', group.by = 'Culture.Protocol', raster=FALSE, cols = color.cul) + plot_annotation(title="Visualization without integration") + NoAxes()
plot2_Unint <- DimPlot(seurat_combi, reduction = 'umap', group.by = 'Sequence.Technology', raster=FALSE, cols = color.tech) + plot_annotation(title="Visualization without integration") + NoAxes()
plot3_Unint <- FeaturePlot(seurat_combi, features=colnames(PlotGeneList), ncol=6, pt.size = 0.5, cols = color.features) & NoAxes() # + plot_annotation(title="Visualization without integration")

print(plot1_Unint)
print(plot2_Unint)
print(plot3_Unint)
print("Finish plotting unintegrated")

# ------------- Harmony Integration -------------
# library(harmony)

seurat_combi$seurat_clusters <- seurat_combi$har_clu
# visualization of the results
plot1_Har <- DimPlot(seurat_combi, group.by='Culture.Protocol', reduction = 'har_umap', pt.size = 0.5, raster=FALSE, cols = color.cul)+ plot_annotation(title="harmony Integration") + NoAxes()
plot2_Har <- DimPlot(seurat_combi, group.by='Sequence.Technology', reduction = 'har_umap', pt.size = 0.5, raster=FALSE, cols = color.tech)+ plot_annotation(title="harmony Integration") + NoAxes()
plot3_Har <- DimPlot(seurat_combi, group.by='har_clu', label = T, reduction = 'har_umap', pt.size = 0.5, raster=FALSE, cols = pal(length(unique(seurat_combi$har_clu)))) + 
              plot_annotation(title="harmony Integration") + 
              NoAxes()
plot4_Har <- FeaturePlot(seurat_combi, features=colnames(PlotGeneList),reduction = 'har_umap', ncol=6, pt.size = 0.5, cols = color.features) & NoAxes() # + plot_annotation(title="harmony Integration") 

print(plot1_Har)
print(plot2_Har)
print(plot3_Har)
print(plot4_Har)

print("Finish Harmony")

# ------------- CSS integration -------------
# library(simspec)

seurat_combi$seurat_clusters <- seurat_combi$css_clu 

#Visualization
plot1_CSS <- DimPlot(seurat_combi, group.by='Culture.Protocol', reduction = 'css_umap',pt.size = 0.5,  raster=FALSE, cols = color.cul)+ plot_annotation(title="CSS Integration") + NoAxes()
plot2_CSS <- DimPlot(seurat_combi, group.by='Sequence.Technology', reduction = 'css_umap',pt.size = 0.5,  raster=FALSE, cols = color.tech)+ plot_annotation(title="CSS Integration") + NoAxes()
plot3_CSS <- DimPlot(seurat_combi, group.by='css_clu', label = T, reduction = 'css_umap', pt.size = 0.5,  raster=FALSE, cols = pal(length(unique(seurat_combi$css_clu)))) + plot_annotation(title="CSS Integration") + NoAxes()
plot4_CSS <- FeaturePlot(seurat_combi, features=colnames(PlotGeneList), reduction = 'css_umap', ncol=6, pt.size = 0.5, cols = color.features) & NoAxes()

print(plot1_CSS)
print(plot2_CSS)
print(plot3_CSS)
print(plot4_CSS)

print("Finish CSS")

# ------------- MNN integration -------------
# library(SeuratWrappers)

seurat_combi$seurat_clusters <- seurat_combi$mnn_clu 

#Visualization
plot1_MNN <- DimPlot(seurat_combi, group.by='Culture.Protocol', reduction = 'mnn_umap', pt.size = 0.5, raster=FALSE, cols = color.cul) + plot_annotation(title="MNN Integration") + NoAxes()
plot2_MNN <- DimPlot(seurat_combi, group.by='Sequence.Technology', reduction = 'mnn_umap', pt.size = 0.5, raster=FALSE,  cols = color.tech) + plot_annotation(title="MNN Integration") + NoAxes()
plot3_MNN <- DimPlot(seurat_combi, group.by='mnn_clu',label = T,reduction = 'mnn_umap',pt.size = 0.5, raster=FALSE, cols = pal(length(unique(seurat_combi$mnn_clu)))) + plot_annotation(title="MNN Integration") + NoAxes()
plot4_MNN <- FeaturePlot(seurat_combi, reduction = 'mnn_umap', features=colnames(PlotGeneList), ncol=6, pt.size = 0.5, cols = color.features) & NoAxes() # + plot_annotation(title="MNN Integration") 

print(plot1_MNN)
print(plot2_MNN)
print(plot3_MNN)
print(plot4_MNN)
print("Finish MNN")

dev.off()
# saveRDS(seurat_combi, 'seurat_1115_ColorToned.rds')
