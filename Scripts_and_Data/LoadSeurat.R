library(Matrix)
library(Seurat)

ProcessData <- function(x){
  Data <- read.table(x)
  M = data.frame(Data[-1,-1], row.names=Data[,1][-1])
  colnames(M) <- Data[1,][-1]
  Mat <- Matrix(as.matrix(M),sparse= T)
  Seurat <- CreateSeuratObject(Mat, project = x)
  return(Seurat)
}

Seu1 <- ProcessData(list.files()[7])
List_Seu <- lapply(list.files()[grep("txt.gz", list.files())], ProcessData)
SeuCombi <- merge(Seu1, List_Seu[-1], project=basename(getwd()))
saveRDS(SeuCombi, file= sprintf("seurat_%s.rds",basename(getwd())))

SeuAA <- Seurat_Pipeline(sprintf("seurat_%s.rds",basename(getwd())))