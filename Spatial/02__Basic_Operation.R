

library(ggplot2)
library(Seurat)


## Load Seurat Object, GBM4@ +tab to see all L2
GBM4 <-Load10X_Spatial(
       data.dir ="GBM4_spaceranger_out/", 
       filename = "filtered_feature_bc_matrix.h5", 
       slice ="GBM4_HE_img")  ## key: GBM4HEimg


## Rename cell names
GBM4$orig.ident <-"GBM4"    ## Default all "SeuratProject"


## ToDo: Filtering


## SCTransform() = NormalizeData() + FindVariableFeatures() + ScaleData()
## Then Clustering and UMAP , be aware of @active.assay
GBM4 <- SCTransform(GBM4, assay = "Spatial", verbose = FALSE)   ## --> GBM4@assays$SCT
GBM4 <- RunPCA(GBM4, assay = "SCT", verbose = FALSE)            ## --> GBM4@reductions$pca
GBM4<- FindNeighbors(GBM4, reduction = "pca", dims = 1:10)      ## construct a KNN graph based on the euclidean distance in 1:10 pca space
GBM4 <- FindClusters(GBM4, verbose = FALSE, resolution = 0.4)   ## GBM4$seurat_clusters
GBM4 <- RunUMAP(GBM4, reduction = "pca", dims = 1:10)           ## --> GBM4@reductions$umap  UMAP of 1:10 pca embeddings


## Plot:   head(GBM4@meta.data) to see current features
## if label = TRUE, it will be Idents(GBM4) or GBM4$seurat_clusters
p1 <- SpatialFeaturePlot(GBM4, features = c("SOX10","nCount_Spatial"))
p2 <- SpatialPlot(GBM4, label = TRUE, label.size = 2) 
p3 <- DimPlot(GBM4, reduction = "umap", label = TRUE)
p4 <- VlnPlot(GBM4, features = c("nCount_Spatial", "nCount_Spatial"), ncol = 2,raster=FALSE)  ## , "mt_percent"


ggsave('../img/02_1.png', p1, width= 6 , height= 3)
ggsave('../img/02_2.png', p2 + p3, width= 6 , height= 3)
ggsave('../img/02_3.png', p4, width= 6 , height= 3)

save(GBM4,file = 'GBM4.rdata')

