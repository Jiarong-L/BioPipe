library(ggplot2)
library(Seurat)
library(dplyr)
load('sc.rdata')


scRNA <- FindNeighbors(scRNA, reduction = "harmony1", dims = 1:30, graph.name = "harmony1_g")
scRNA <- FindClusters(scRNA, graph.name = "harmony1_g", resolution = 0.6 ) ## default output to Idents(scRNA) / scRNA@meta.data$seurat_clusters, or you can set cluster.name = "harmony1_cluster"
p1 <- DimPlot(scRNA, reduction = "umap_harmony1",group.by = "seurat_clusters",label = T)


## Just an example,not true!! 
## Assign Cell types to clusters by the expression of known markers & check it via umap
# markers <- FindAllMarkers(scRNA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# top10 <- markers %>%
#   group_by(cluster) %>%
#   top_n(n = 10, wt = avg_log2FC)


Marker = list('A like'=c('OLIG1'),
              'B like'=c('HOPX'),
              'C like' =c('PTGDS'))

p2 <- DotPlot(scRNA,features=Marker,cols = c('gray','red')) + theme(axis.text.x = element_text(angle = 90, hjust = 1))


scRNA@meta.data$cell_type <- NA
scRNA@meta.data$cell_type[scRNA@meta.data$seurat_clusters %in% c('0','1','6','8')] <- "A like"
scRNA@meta.data$cell_type[scRNA@meta.data$seurat_clusters %in% c('2','3','5','7','9')] <- "B like"
scRNA@meta.data$cell_type[scRNA@meta.data$seurat_clusters %in% c('4','10','11')] <- "C like"
Idents(scRNA)<-scRNA$cell_type

p3 <- DimPlot(scRNA, reduction = "umap_harmony1",label = T)

ggsave('../img/04_2.png', p1 + p2 + p3, width= 15 , height= 4)


save(scRNA,file = 'scRNA.rdata')
