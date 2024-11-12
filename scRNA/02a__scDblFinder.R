
### scDblFinder的步骤和原理
if(FALSE) {
"
    也是生成人工Doublet,但是为了训练区分Doublet和Singlet的Classifier模型,于是可以使用更多feature(Figure 2)


SingleCellExperiment 结构: https://bioconductor.org/books/3.14/OSCA.intro/the-singlecellexperiment-class.html
Mannual: https://bioconductor.org/packages/release/bioc/vignettes/scDblFinder/inst/doc/scDblFinder.html
Mannual: https://github.com/plger/scDblFinder
原理: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9204188/
"
}




library(Seurat)
library(scDblFinder)
library(stringr)
library(ggplot2)
source("../scProcess_SCT.R")

counts <- Read10X(data.dir = "pbmc3k/",gene.column = 2)   ## 1:ENSG00xxx   2:geneName
pbmc3k <- CreateSeuratObject(counts,  min.cells=3, min.features = 200)
pbmc3k <- scProcess_SCT(pbmc3k)

## 多样本的情况，Mannual中搜索关键词 BPPARAM
sce <- as.SingleCellExperiment(pbmc3k, assay ="SCT" ) ## 设置DefaultAssay(pbmc)没用
sce <- scDblFinder(sce, clusters="seurat_clusters")  ## clusters=NULL即生成随机的人工Doublet；
pbmc3k <- as.Seurat(sce)                             ## dbr每1k细胞加1%，不过建议还是参考10X(见DoubletFinder.R中链接)



p2 <- DimPlot(pbmc3k, reduction = "UMAP", group.by = "scDblFinder.class") ## colnames(pbmc3k@meta.data)，且UMAP变大写了
ggsave('../img/02a_2.png', p2, width= 5 , height= 3)

table(sce$scDblFinder.class)
# singlet doublet
#    2670      27