
### DecontX的步骤和原理
if(FALSE) {
"
使用贝叶斯模型, 不需要 empty droplet information, 也可拆分出单个细胞中的环境RNA分布

decontX()适用于scRNA:  
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1950-6
DecontPro()适用于CITE-seq ADT:  
https://pmc.ncbi.nlm.nih.gov/articles/PMC9979990/
"
}



library(Seurat)
library(decontX)
library(ggplot2)
source("../scProcess_SCT.R")


counts <- Read10X(data.dir = "pbmc3k/",gene.column = 2)   ## 1:ENSG00xxx   2:geneName
pbmc3k <- CreateSeuratObject(counts,  min.cells=3, min.features = 200)

## decontX + filterX
decontX_res <- decontX(pbmc3k[["RNA"]]$counts)    ## pbmc3k@assays$RNA@counts    (filtered min cells/counts)
pbmc3k$contam <- decontX_res$contamination        ## runParams/decontXcounts(矫正后)/z/estimates/contamination(污染度)
pbmc3k_filt <- pbmc3k[,pbmc3k$contam<0.2]         ## 2700 -> 2499 cells


## show diff (not obvious in this case)
pbmc3k <- scProcess_SCT(pbmc3k)
pbmc3k_filt <- scProcess_SCT(pbmc3k_filt)
p1 <- DimPlot(pbmc3k, reduction = "umap") + DimPlot(pbmc3k_filt, reduction = "umap")
ggsave('../img/02c_1.png', p1, width= 5 , height= 3)

