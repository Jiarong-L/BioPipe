
### SingleR自动注释细胞
if(FALSE) {
"
基于test/ref中重合的基因, 计算其中的Marker以进行注释

ref_DB格式---SummarizedExperiment   https://zhuanlan.zhihu.com/p/55001703


参考: https://cloud.tencent.com/developer/article/2194308
SingleR 自带数据库: 
人-- HumanPrimaryCellAtlasData() 
     BlueprintEncodeData()
     DatabaseImmuneCellExpressionData()
     NovershternHematopoieticData()
     MonacoImmuneData()
鼠 -- MouseRNAseqData() / ImmGenData()

常用数据库: cellMarker/PanglaoDB/CancerSEA
"
}



library(SingleR)
library(celldex)
library(Seurat)
library(ggplot2)
library(patchwork)
source("../scProcess_SCT.R")


counts = Read10X(data.dir = "pbmc3k/",gene.column = 2)   ## 1:ENSG00xxx   2:geneName
pbmc3k = CreateSeuratObject(counts,  min.cells=3, min.features = 200)
pbmc3k = scProcess_SCT(pbmc3k)


testMtx = GetAssayData(pbmc3k@assays$SCT, slot="data")  ## Normalized Counts
refSE = celldex::HumanPrimaryCellAtlasData() 
labels = refSE$label.main            ## label.main/fine/ont
clusters = pbmc3k@meta.data$seurat_clusters


SingleR_Anno = SingleR(
        test = testMtx,
        ref = refSE,
        labels = labels,
        clusters = clusters)


## Idents(pbmc3k)
pbmc3k@meta.data$Anno_1 = SingleR_Anno$labels[match(clusters,rownames(SingleR_Anno))] 



# > SingleR_Anno
# DataFrame with 5 rows and 4 columns
#                           scores      labels delta.next pruned.labels
#                         <matrix> <character>  <numeric>   <character>
# 0 0.215783:0.573543:0.515419:...     T_cells  0.0794487       T_cells
# 1 0.158448:0.515736:0.576116:...    Monocyte  0.4425706      Monocyte
# 2 0.203510:0.704584:0.563297:...      B_cell  0.0588963        B_cell
# 3 0.207598:0.534162:0.497759:...     NK_cell  0.3000568       NK_cell
# 4 0.146295:0.262188:0.343418:...   Platelets  0.1515570     Platelets


p1 = plotScoreHeatmap(SingleR_Anno)   ## plotDeltaDistribution(SingleR_Anno, ncol = 5)
ggsave('../img/03_1.png', p2, width= 8 , height= 5)
