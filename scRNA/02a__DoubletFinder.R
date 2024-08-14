
### DoubletFinder的步骤和原理
if(FALSE) {
"

1. 根据预先定义的细胞类型, 模拟生成一些Doublet
2. 混合真实-模拟数据
3. 对混合数据进行PCA, 原则上人工模拟的Doublet会与真实的Doublet距离较近, 因此邻居中人工Doublet比例(pANN)较多者更可能是Doublet
4. 根据10X建议设定Doublet的总数, 计算其中Heterotypic的期望数,(假设数据符合泊松分布),如此可知(pANN)的cut-off值
(DoubletRate随着上机细胞的浓度而改变, 数据符合Poisson分布)

## pN指定模拟双细胞的比例: 真实/模拟      一般设定  =25%         
## pK指定近邻的比例: 将最近的 ALL*pK 个细胞视为近邻    一般对多个pN-pK组合计算pANN、择优选择使pANN呈双峰的pK
## pANN: 近邻中, 模拟双细胞的比例
## Homotypic 同型来源双细胞主要来自模拟, 其比例 pHomo = 每个Cell_Type_频率的平方和
## Heterotypic 异型来源双细胞是期望去除是部分, 其比例 pHeter = 1 - pHomo 

DoubletRate Estimation: 
A:  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6126471/
B:  https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled

Mannual: https://blaserlab.r-universe.dev/DoubletFinder/doc/manual.html
"
}




library(Seurat)
library(DoubletFinder)
library(stringr)
library(ggplot2)


counts <- Read10X(data.dir = "pbmc3k/",gene.column = 1)
pbmc3k <- CreateSeuratObject(counts,  min.cells=3, min.features = 200)

pbmc3k <- SCTransform(pbmc3k, assay = "RNA", verbose = FALSE) 
pbmc3k <- RunPCA(pbmc3k, assay = "SCT", verbose = FALSE)        
pbmc3k <- RunUMAP(pbmc3k, reduction = "pca", dims = 1:30)             
pbmc3k <- RunTSNE(pbmc3k, reduction = "pca", dims = 1:30)          
pbmc3k <- FindNeighbors(pbmc3k, reduction = "pca", dims = 1:30)     
pbmc3k <- FindClusters(pbmc3k, verbose = FALSE, resolution = 0.1)  



### 计算多种pN-pK参数组合对应的双峰系数，以寻找最优pK
sweep.list <- paramSweep(pbmc3k, PCs = 1:30, sct = T)      ##（10x上限是 10k cell）pN = 0.05-0.3, pK = 0.0005-0.3
sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)      ## 获得(pN-pK-BCreal)，若提供真实的双细胞标签GT、还可以ROC
bcmvn <- find.pK(sweep.stats)                              ## 获得(pK--BCmetric)
pK_bcmvn <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))   ## 最优pK


### 假设有0.15的细胞是Doublet，其中Homotypic部分是模拟引入的，Heterotypic部分是真实待去除的
DoubletRate = 0.15    ## This value can best be estimated from cell loading densities into the 10X/Drop-Seq device
AllCellNum = nrow(pbmc3k@meta.data)
Annotations = pbmc3k$seurat_clusters

nExp_poi <- round(DoubletRate*AllCellNum)
pHomo <- modelHomotypic(Annotations)
nExp_poi.adj <- round(nExp_poi*(1-pHomo))


### DoubletFinder
sig_PC = 1:10
pbmc3k <- doubletFinder(pbmc3k, PCs = sig_PC, pN = 0.25, pK = pK_bcmvn, 
                        nExp = nExp_poi.adj, 
                        reuse.pANN = FALSE,   ## pANN_x_y_z/FALSE, stored as pbmc3k@meta.data column
                        sct = TRUE) 


### Plot Singlet/Doublet 
p1 <- DimPlot(pbmc3k, reduction = "umap", group.by = "DF.classifications_0.25_0.12_289") ## colnames(pbmc3k@meta.data)
ggsave('../img/02a_1.png', p1, width= 5 , height= 3)


sum(pbmc3k@meta.data$DF.classifications_0.25_0.12_289 == "Singlet")  ## 2700 --> 2411
