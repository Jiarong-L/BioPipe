### SoupX的步骤和原理
if(FALSE) {
"
需要 raw_feature_bc_matrix/ 文件夹; 污染率需要用户基于组织坏死状况预估!!


一些简写
b:background--(空液滴<N_emp UMIs)
C:cell
G:gene
o(G,C) observed counts of gene_G in cell_C
o(C) = o(AllG,C) observed counts in cell_C


原理
1. 依据背景液滴计算 Frac_b(G) = b(G)/b(AllG) = background_counts_of_gene_G / background_counts_All
2. 聚类后, 计算每一种Cell的污染比例 Frac_contam(C)
    - 假设: o(G,C) = endo(G,C) + b(G,C)
    - 假设: b(G,C) = o(C) * Frac_contam(C) * Frac_b(G)
    - 已知: 一些细胞一定不会含有某些基因, 例如 HB 基因不应该在红细胞之外出现 (Auto/用户设定基因集)
        o(G,C) = b(G,C)  <==  endo(G,C) = 0
        o(G,C) = o(C) * Frac_contam(C) * Frac_b(G)
        Frac_contam(C) = o(G,C) / [ o(C) * Frac_b(G) ]         /// for gene_G only, below avg for all genes
        Frac_contam(C) = o(AllG,C) / [ o(C) * sum_of{Frac_b(G) for G in AllG} ] 
3. 依次从每一种细胞的表达中去除背景噪音 b(G,C), 留下 endo(G,C)


数据
wget https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_raw_gene_bc_matrices.tar.gz
tar -xzf pbmc3k_raw_gene_bc_matrices.tar.gz; ls raw_gene_bc_matrices/hg19
mkdir filtered_gene_bc_matrices; ln -s ../pbmc3k filtered_gene_bc_matrices/hg19


步骤参考
https://cran.r-project.org/web//packages/SoupX/SoupX.pdf
https://github.com/constantAmateur/SoupX
https://blog.csdn.net/qq_44933752/article/details/129392901
https://www.jianshu.com/p/b9799b401f8f
"
}


library(Seurat)
library(SoupX)
library(stringr)
library(ggplot2)
source("../scProcess_SCT.R")


scProcess_SoupX <- function(counts){
    SeuratObj = CreateSeuratObject(counts)
    SeuratObj = NormalizeData(SeuratObj)
    SeuratObj = FindVariableFeatures(SeuratObj, selection.method = "vst", nfeatures = 3000)
    SeuratObj = ScaleData(SeuratObj)#, features = rownames(SeuratObj))
    SeuratObj = scRoutine(SeuratObj)
    return(SeuratObj)
}





## load data , not suggest sc = load10X('.') 
raw_counts = Seurat::Read10X('./raw_gene_bc_matrices/hg19')
filt_counts = Seurat::Read10X('./filtered_gene_bc_matrices/hg19')
raw_counts <- raw_counts[rownames(filt_counts),]
rownames(raw_counts) = str_replace_all(rownames(raw_counts), "_", "-")
rownames(filt_counts) = str_replace_all(rownames(filt_counts), "_", "-")

## get cluster info
pbmc3k_filt = scProcess_SoupX(filt_counts)
cluster_labels = setNames(pbmc3k_filt@meta.data$seurat_clusters, rownames(pbmc3k_filt@meta.data))


######################### AutoPipe: Create SoupChannel object, it fits into dim(sc$toc)==>dim(filt_counts)
sc = SoupChannel(raw_counts,filt_counts)
sc = setClusters(sc,cluster_labels)
sc = autoEstCont(sc)
# // setContaminationFraction()
SoupX_adjusted = adjustCounts(sc)


######################### AdvancePipe: User set up the geneSet        names(sc) 
sc = SoupChannel(raw_counts,filt_counts)
sc = setClusters(sc,cluster_labels)

genes = c('CD79A','MS4A1','TCL1A','CST7','LINC00926','GZMA','HLA-DQA1','VPREB3','CD79B','PRF1','HLA-DQB1','NKG7','GZMH','FGFBP2','CFD','LGALS2','GZMB','BANK1','CD68','HLA-DQA2','IFI30','FCN1','S100A8','CTSW','CFP','SPI1','SERPINA1','HVCN1','CCL4','IFITM3','MS4A6A','CDA','GRN')  # get from codes of autoEstCont tgts, as -- str_subset(rownames(sc$toc), "^HB") too extreme for pbmc
gene_set_User = as.list(genes)
names(gene_set_User) = genes

pickup_cells = estimateNonExpressingCells(sc,gene_set_User) 
outStr = capture.output({
    calculateContaminationFraction(sc,gene_set_User,pickup_cells,forceAccept=TRUE)
}, type = "message")

Frac_contam_str = gsub("Estimated global contamination fraction of .*?(\\d+\\.\\d+).*", "\\1", paste(outStr, collapse = ", "))
Frac_contam = as.numeric(Frac_contam_str)/100
sc = setContaminationFraction(sc,Frac_contam)   

SoupX_adjusted = adjustCounts(sc)


