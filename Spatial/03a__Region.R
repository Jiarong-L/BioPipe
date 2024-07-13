
library(ggplot2)
library(Seurat)
load('GBM4.rdata')



GBM4@meta.data$Region<-NA
GBM4@meta.data$Region[GBM4@meta.data$seurat_clusters %in% c('1','5')] <- "Normal"
GBM4@meta.data$Region[GBM4@meta.data$seurat_clusters %in% c('3','4','2')] <- "Transition"
GBM4@meta.data$Region[GBM4@meta.data$seurat_clusters %in% c('0')] <- "Tumor"
 
## or:  Idents(GBM4)<-GBM4$Region  to update   GBM4@meta.data$seurat_clusters
p1 <- SpatialPlot(GBM4, label = FALSE, label.size = 5,
                  group.by = 'Region',
                  cols = c('Normal'='#4b5cc4','Transition'='#FE8D3C','Tumor'='#AA0000'))

ggsave('../img/03_1.png', p1, width= 6 , height= 3)


save(GBM4,file = 'GBM4.rdata')

