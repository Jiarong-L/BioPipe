library(ggplot2)
library(Seurat)
library(dplyr)
library(clusterProfiler)
load('GBM4.rdata')



## Find top10 markers for regions
Idents(GBM4)<-GBM4$Region
markers <- FindAllMarkers(GBM4, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top10<-markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)


## Heatmap 

p2 <- DoHeatmap(GBM4, features = top10$gene,group.colors = c('Normal'='#4b5cc4','Transition'='#FE8D3C','Tumor'='#AA0000')) + NoLegend()

ggsave('../img/03_2.png', p2, width= 6 , height= 3)


## BiocManager::install("org.Hs.eg.db")
## SYMBOL to ENTREZID
gid <- bitr(unique(markers$gene), 'SYMBOL', 'ENTREZID', OrgDb= 'org.Hs.eg.db')
markers <- full_join(markers, gid, by=c('gene' = 'SYMBOL'))

KEGG = compareCluster(ENTREZID ~ cluster, data = markers, fun='enrichKEGG')

p3 <- dotplot(KEGG, label_format=40) + theme(axis.text.x = element_text(angle=45, hjust=1)) + scale_color_gradient(high="#4b5cc4",low="#FE8D3C")
ggsave('../img/03_3.png', p3, width= 6 , height= 3)





save(GBM4,file = 'GBM4.rdata')