library(ggplot2)
library(Seurat)
library(dplyr)
library(harmony)
load('sc.rdata')


## Harmony: group.by.vars should be batch_labels (you would like to eliminate diff between these groups)
scRNA <- RunHarmony(scRNA,reduction = "pca",group.by.vars = "Source",reduction.save = "harmony1")


## Harmony  to integrate layers  (same as before?)
scRNA <- IntegrateLayers(
  object = scRNA, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony2",
  verbose = FALSE
)


## CCA  when cell_type similar across samples   slow
scRNA <- IntegrateLayers(
  object = scRNA, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "CCA",
  verbose = FALSE
)

## RPCA  when cell_type diff across samples   slowerrr
scRNA <- IntegrateLayers(
  object = scRNA, method = RPCAIntegration ,
  orig.reduction = "pca", new.reduction = "RPCA",
  verbose = FALSE
)

scRNA <- JoinLayers(scRNA)  

## Check Batch Effect before/after correction
scRNA <- RunUMAP(scRNA, reduction = "pca", dims = 1:30, verbose = F, reduction.name = "umap_pca")
scRNA <- RunUMAP(scRNA, reduction = "harmony1", dims = 1:30, verbose = F, reduction.name = "umap_harmony1")
scRNA <- RunUMAP(scRNA, reduction = "harmony2", dims = 1:30, verbose = F, reduction.name = "umap_harmony2")
scRNA <- RunUMAP(scRNA, reduction = "CCA", dims = 1:30, verbose = F, reduction.name = "umap_CCA")
scRNA <- RunUMAP(scRNA, reduction = "RPCA", dims = 1:30, verbose = F, reduction.name = "umap_RPCA")

p1 <- DimPlot(scRNA,reduction = "umap_pca", group.by = "Source")
p2 <- DimPlot(scRNA,reduction = "umap_harmony1", group.by = "Source")
p3 <- DimPlot(scRNA,reduction = "umap_harmony2", group.by = "Source")
p4 <- DimPlot(scRNA,reduction = "umap_CCA", group.by = "Source")
p5 <- DimPlot(scRNA,reduction = "umap_RPCA", group.by = "Source")

ggsave('../img/04_1.png', p1 + p2 + p3+ p4 +p5, width= 10 , height= 5)


save(scRNA,file = 'sc.rdata')