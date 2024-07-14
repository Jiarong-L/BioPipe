library(ggplot2)
library(Seurat)
library(CARD)
load('GBM4.rdata')
load('sc.rdata')


scRNA <- FindNeighbors(scRNA, reduction = "CCA", dims = 1:30, graph.name = "CCA_g")
scRNA <- FindClusters(scRNA, graph.name = "CCA_g", resolution = 0.6, cluster.name = "CCA_cluster") 


DimPlot(scRNA, reduction = "umap_CCA",group.by = "CCA_cluster",label = T)




