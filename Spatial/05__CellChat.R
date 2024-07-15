library(ggplot2)
library(Seurat)
library(CellChat)
load('GBM4.rdata')



GBM4@meta.data$cell_type <- paste("C",GBM4$seurat_clusters,sep = "")
Idents(GBM4) <- "cell_type"    ## Previous is GBM4@meta.data$Region
p1 <- SpatialDimPlot(GBM4, label = TRUE, label.size = 5)





