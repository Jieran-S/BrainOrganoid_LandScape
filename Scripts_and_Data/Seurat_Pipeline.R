#This is the pipeline used for treating individual single cell sequences from the collected published papers listed in
# the folder. The pipeline follows the instruction from the single cell pipeline instruction.
# The scipt will take the directory of folder which contains the seurat object file as input and use it to extract and 
# analyze the data.
#The output will explore MNN, CSS and Harmony integration and listed them all 
#in separate pages of pdf for the final results. 

library(Seurat)
library(patchwork)
library(dplyr)
library(magrittr)
library(Matrix)
library(ggplot2)

Seurat_Pipeline <- function(x){

  #Import Seurat object (combined)
  seurat_combi <- readRDS(x) 
  seurat_combi[["percent.mt"]] <- PercentageFeatureSet(seurat_combi, pattern = "^MT[-\\.]")
  seurat_combi$log10GenesPerUMI <- log10(seurat_combi$nFeature_RNA)/ log10(seurat_combi$nCount_RNA)
  
  #Start PDF
  pdf(file="Integration_Visualization.pdf", 
      width = 15,
      height = 7.5)
  
  seurat_combi <- subset(seurat_combi, ((nCount_RNA >0)&(nFeature_RNA >0)&(log10GenesPerUMI>0)))
  
  #Firstly the violin plot of the cell counts and feature counts
  Vplot1  <- VlnPlot(seurat_combi, features = c("nFeature_RNA"), pt.size = 0) +
                    geom_hline(yintercept = 10^( mean(log10(seurat_combi$nFeature_RNA)) - 2*sd(log10(seurat_combi$nFeature_RNA)))) +
                    geom_hline(yintercept = 10^( mean(log10(seurat_combi$nFeature_RNA)) + 2*sd(log10(seurat_combi$nFeature_RNA))))
  
  Vplot2  <- VlnPlot(seurat_combi, features = c("nCount_RNA"), pt.size = 0) +
                    geom_hline(yintercept = 10^( mean(log10(seurat_combi$nCount_RNA)) - 2*sd(log10(seurat_combi$nCount_RNA)))) +
                    geom_hline(yintercept = 10^( mean(log10(seurat_combi$nCount_RNA)) + 2*sd(log10(seurat_combi$nCount_RNA))))
  
  if (all(seurat_combi$percent.mt == seurat_combi$percent.mt[1], na.rm = T)){
      Vplot <- Vplot1 + Vplot2 + plot_layout( ncol=2 ) + plot_annotation(title = basename(getwd()))
  }else{
    Vplot3 <- VlnPlot(seurat_combi, features = c("percent.mt"), pt.size = 0) +
      geom_hline(yintercept =  mean(seurat_combi$percent.mt) + 2*sd(seurat_combi$percent.mt))
    
    Vplot <- Vplot1 + Vplot2 + Vplot3 + plot_layout( ncol=3 ) + plot_annotation(title = basename(getwd()))
  }
  
  print(Vplot)
  
  #cell-wise QC to eliminate cells with too much or too little UMI, genes, mitoRNA and cell complexity value 
  #Doing QC: Log10 transform the data then assuming it Gaussian Distribution, filtering out all data that lies outside 2 std of the nCount_RNA
  if (all(seurat_combi$percent.mt == seurat_combi$percent.mt[1], na.rm = T)){
    if (is.na( mean(log10(seurat_combi$nCount_RNA)+2*sd(log10(seurat_combi$nCount_RNA)))) 
        | is.na(mean(log10(seurat_combi$nFeature_RNA)) + 2*sd(log10(seurat_combi$nFeature_RNA)))){
          seurat_combi <- subset(x=seurat_combi,
                                 subset = ((nCount_RNA > 100) & 
                                           (nCount_RNA < 30000) & 
                                           (nFeature_RNA > 100 )& 
                                           (nFeature_RNA < 15000) &
                                           (log10GenesPerUMI > 0.2)))
    }else{          
          seurat_combi <- subset(x=seurat_combi,
                                 subset = ((log10(seurat_combi$nCount_RNA) > mean(log10(seurat_combi$nCount_RNA)) - 2*sd(log10(seurat_combi$nCount_RNA))) & 
                                             (log10(seurat_combi$nCount_RNA) < mean(log10(seurat_combi$nCount_RNA)) + 2*sd(log10(seurat_combi$nCount_RNA))) & 
                                             (log10(seurat_combi$nFeature_RNA) > mean(log10(seurat_combi$nFeature_RNA)) - 2*sd(log10(seurat_combi$nFeature_RNA))) & 
                                             (log10(seurat_combi$nFeature_RNA) < mean(log10(seurat_combi$nFeature_RNA)) + 2*sd(log10(seurat_combi$nFeature_RNA))) &
                                             (log10GenesPerUMI > mean(log10(seurat_combi$log10GenesPerUMI)) - 2*sd(log10(seurat_combi$log10GenesPerUMI))))
          )
    }
  }else if (is.na( mean(log10(seurat_combi$nCount_RNA)+2*sd(log10(seurat_combi$nCount_RNA)))) | is.na(mean(log10(seurat_combi$nFeature_RNA)) + 2*sd(log10(seurat_combi$nFeature_RNA)))){
    seurat_combi <- subset(x=seurat_combi,
                          subset = ((nCount_RNA > 100) & 
                                    (nCount_RNA < 30000) & #What should the threhold set to be?
                                    (nFeature_RNA > 100 )& 
                                    (nFeature_RNA < 15000) &
                                    (percent.mt < 20)& 
                                    (log10GenesPerUMI > 0.2)))
                  
  }else{
    seurat_combi <- subset(x=seurat_combi,
                           subset = ((log10(seurat_combi$nCount_RNA) > mean(log10(seurat_combi$nCount_RNA)) - 2*sd(log10(seurat_combi$nCount_RNA))) & 
                                       (log10(seurat_combi$nCount_RNA) < mean(log10(seurat_combi$nCount_RNA)) + 2*sd(log10(seurat_combi$nCount_RNA))) & 
                                       (log10(seurat_combi$nFeature_RNA) > mean(log10(seurat_combi$nFeature_RNA)) - 2*sd(log10(seurat_combi$nFeature_RNA))) & 
                                       (log10(seurat_combi$nFeature_RNA) < mean(log10(seurat_combi$nFeature_RNA)) + 2*sd(log10(seurat_combi$nFeature_RNA))) &
                                       (percent.mt < min(mean(percent.mt) + 2*sd(percent.mt), 20)) & 
                                       (log10GenesPerUMI > mean(log10(seurat_combi$log10GenesPerUMI)) - 2*sd(log10(seurat_combi$log10GenesPerUMI)))))
  }

  
  #Gene filtering: Leaving only features with at least 10 true values detected among all cells
  counts <- GetAssayData(object = seurat_combi, slot = "counts")
  nonzero <- counts > 0
  #leaving only genes with expression positive in at least 10 cells
  keep_genes <- Matrix::rowSums(nonzero) >= 10
  #reshape the count matrix
  filtered_counts <- counts[keep_genes, ]
  #assign the new matrix into the seurat object
  seurat_combi <- CreateSeuratObject(filtered_counts, meta.data = seurat_combi@meta.data)
  
  seurat_combi <- seurat_combi %>%
    NormalizeData() %>%
    FindVariableFeatures(nfeatures = 3000) %>%
    ScaleData() %>%
    RunPCA(npcs = 50)
  #Scaled, normalized and PCA the gene data 
  
  # PC selection: Choosing the combo either the PC contributes less than 5% whereas cummulative sum over 80%
  # or that the the percent change in variation between the consecutive PCs is less than 0.1%.
  # or that the number of PC selected is over 25 (half of the 50)
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
                  label.size = 0.35,
                  color = "black",
                  fill="#69b3a2"
                )
  print(plot_elbow)
  
  #First time visualization of the dataset in terms of different sources
  plot1_Unint <- DimPlot(seurat_combi, group.by = "orig.ident")
  plot2_Unint <- FeaturePlot(seurat_combi, c("MKI67","NES","DCX","FOXG1","DLX2","EMX1","OTX2","LHX9","nFeature_RNA"), ncol=3, pt.size = 0.5)
  plot3_Unint <- plot1_Unint + plot2_Unint + plot_layout(widths = c(1.5, 2)) + plot_annotation(title="Visualization without integration")
  print(plot3_Unint)
  print("Finish plotting unintegrated")
  
  
  #Harmony Integration
  library(harmony)
  
  seurat_combi <- RunHarmony(seurat_combi, group.by.vars = "orig.ident", dims.use = 1:pcs, max.iter.harmony = 50)
  #seurat_combi[['harmony']] <- CreateDimReducObject(Embeddings(seurat_combi, "harmony")[colnames(seurat_combi),], key="HAR_")
  seurat_combi <- RunUMAP(seurat_combi, reduction = "harmony", dims = 1:pcs, reduction.name = 'har_umap',reduction.key="UMAPHAR_")
  seurat_combi <- FindNeighbors(seurat_combi, reduction = "harmony", dims = 1:pcs, graph.name = 'har_snn') %>% FindClusters(resolution = 0.6)
  
  #visualization of the results
  plot1_Har <- UMAPPlot(seurat_combi, group.by="orig.ident", pt.size = 0.5)
  plot2_Har <- UMAPPlot(seurat_combi, label = T, pt.size = 0.5)
  plot3_Har <- FeaturePlot(seurat_combi, c("MKI67","NES","DCX","FOXG1","DLX2","EMX1","OTX2","LHX9","nFeature_RNA"), ncol=3, pt.size = 0.5)
  plot_Har <- ((plot1_Har / plot2_Har) | plot3_Har) + plot_layout(width = c(1,2)) + plot_annotation(title = "harmony integration")
  print(plot_Har)
  
  print("Finish Harmony")
  
  #CSS integration
  library(simspec)
  
  seurat_combi <- cluster_sim_spectrum(seurat_combi, label_tag = "orig.ident", cluster_resolution = 0.3)
  #seurat_combi[['css']] <- CreateDimReducObject(Embeddings(seurat_combi, "css")[colnames(seurat_combi),], key="CSS_")
  seurat_combi <- RunUMAP(seurat_combi, reduction="css", dims = 1:ncol(Embeddings(seurat_combi, "css")), reduction.name = 'css_umap',reduction.key="UMAPCSS_")
  seurat_combi <- FindNeighbors(seurat_combi, reduction = "css", dims = 1:ncol(Embeddings(seurat_combi, "css")), graph.name = 'css_snn') %>%
    FindClusters(resolution = 0.6)
  
  #Visualization
  plot1_CSS <- UMAPPlot(seurat_combi, group.by="orig.ident", pt.size = 0.5)
  plot2_CSS <- UMAPPlot(seurat_combi, label = T, pt.size = 0.5)
  plot3_CSS <- FeaturePlot(seurat_combi, c("MKI67","NES","DCX","FOXG1","DLX2","EMX1","OTX2","LHX9","nFeature_RNA"), ncol=3, pt.size = 0.5)
  plot_CSS <- ((plot1_CSS / plot2_CSS) | plot3_CSS) + plot_layout(width = c(1,2)) + plot_annotation(title = "CSS integration")
  print(plot_CSS)
  #dev.off()
  print("Finish CSS")
  
  #if only one object then skip MNN
  if(all(seurat_combi$orig.ident == seurat_combi$orig.ident[1], na.rm = T)){
    dev.off()
  }else{
  #MNN integration
  library(SeuratWrappers)
  
  seurat_samples <- SplitObject(seurat_combi, "orig.ident")
  seurat_mnn <- RunFastMNN(seurat_samples)
  seurat_combi[['mnn']] <- CreateDimReducObject(Embeddings(seurat_mnn, "mnn")[colnames(seurat_combi),], key="MNN_")
  seurat_combi <- RunUMAP(seurat_combi, dims = 1:pcs, reduction = "mnn", reduction.name = 'mnn_umap',reduction.key="UMAPMNN_" )
  seurat_combi <- FindNeighbors(seurat_combi, reduction = "mnn", dims = 1:pcs, graph.name = 'mnn_snn' ) %>%
    FindClusters(resolution = 0.6)
  
  #Visualization
  plot1_MNN <- UMAPPlot(seurat_combi, group.by="orig.ident", pt.size = 0.5)
  plot2_MNN <- UMAPPlot(seurat_combi, label = T, pt.size = 0.5)
  
  #Change the width and height of pdf files.
  plot3_MNN <- FeaturePlot(seurat_combi, c("MKI67","NES","DCX","FOXG1","DLX2","EMX1","OTX2","LHX9","nFeature_RNA"), ncol=3, pt.size = 0.1)
  plot_MNN <- ((plot1_MNN / plot2_MNN) | plot3_MNN) + plot_layout(width = c(1,2)) + plot_annotation(title = "MNN integration")
  print(plot_MNN)
  print("Finish MNN")
  dev.off()
  }
  
  saveRDS(seurat_combi, sprintf("seurat_%s.rds",basename(getwd())))
  return(seurat_combi)
}