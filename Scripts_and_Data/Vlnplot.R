library(Seurat)
library(patchwork)
library(ggplot2)

remove = c("2019_Yoon_BD_New", "2021_Bowles_Drop_New")

for (file in setdiff(list.files(), remove)){
  seurat_combi <- readRDS(sprintf("%s/seurat_%s.rds",file, file)) 
  
  #Start PDF
  pdf(file=sprintf("%s/Vlnplot_%s.pdf", file, file), 
      width = 15,
      height = 30)
  
  #Firstly the violin plot of the cell counts and feature counts
  Vplot1  <- VlnPlot(seurat_combi, features = c("nFeature_RNA"),group.by="orig.ident", pt.size = 0)
  
  Vplot2  <- VlnPlot(seurat_combi, features = c("nCount_RNA"), group.by="orig.ident", pt.size = 0)
  
  if (all(seurat_combi$percent.mt == seurat_combi$percent.mt[1], na.rm = T)){
    Vplot <- Vplot1 + Vplot2 + plot_layout( nrow=2 ) + plot_annotation(title = file)
    
  }else{
    Vplot3 <- VlnPlot(seurat_combi, features = c("percent.mt"),group.by="orig.ident",  pt.size = 0) 
    Vplot <- Vplot1 + Vplot2 + Vplot3 + plot_layout( nrow=3 ) + plot_annotation(title = file)
  }
  print(Vplot)
  dev.off()
}
  