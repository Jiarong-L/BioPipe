
library(Seurat)


scFilter <- function(SeuratObj){        
    # 管家基因：受环境因素影响较小，不过此处不指定

    # 临近死亡的细胞会有很多线粒体
    SeuratObj[["pct.mt"]] <- PercentageFeatureSet(SeuratObj, pattern = "^MT-")
    # 识别红细胞，有时组织样本会有残留的血液
    SeuratObj[["pct.hb"]] <- PercentageFeatureSet(SeuratObj, pattern = "^HB[^(P)]")
    # 核糖体负责蛋白翻译，有些文章为了消除批次效应会直接删掉这些基因
    SeuratObj[["pct.rp"]] <- PercentageFeatureSet(SeuratObj, pattern = "^RP[SL]")
    # 事实上阈值应当根据图来估计 VlnPlot(SeuratObj, features = c("nFeature_RNA"), ncol = 1)
    SeuratObj <- subset(SeuratObj, subset = nFeature_RNA > 200 & pct.mt < 15 & pct.hb < 1)  ##  & pct.rp < 5& nFeature_RNA < 2500 
    return(SeuratObj)
}


scCellCycle <- function(SeuratObj){          
    s.genes <- CaseMatch(search=cc.genes$s.genes,match=rownames(SeuratObj))
    g2m.genes <- CaseMatch(search=cc.genes$g2m.genes,match=rownames(SeuratObj))
    SeuratObj <- CellCycleScoring(SeuratObj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
    SeuratObj$CC.Difference <- SeuratObj$S.Score - SeuratObj$G2M.Score
    ## RidgePlot(SeuratObj, features = g2m.genes[1:2], ncol = 2)
    return(SeuratObj)
}

scRoutine <- function(SeuratObj){     # 事实上PC数量应当根据图来估计，此处偷懒   
    SeuratObj <- RunPCA(SeuratObj, verbose = FALSE)    ##  , assay = "SCT"   
    SeuratObj <- RunUMAP(SeuratObj, reduction = "pca", dims = 1:30)             
    SeuratObj <- RunTSNE(SeuratObj, reduction = "pca", dims = 1:30)          
    SeuratObj <- FindNeighbors(SeuratObj, reduction = "pca", dims = 1:30, ) 
              ## FindNeighbors(SeuratObj, features = VariableFeatures(object = SeuratObj))    
              ## pbmc3k@neighbors/pbmc3k@graphs 只能返回一个，因为 return.neighbor =TRUE 与SNN互斥
    SeuratObj <- FindClusters(SeuratObj, verbose = FALSE, resolution = 0.1)  
    return(SeuratObj)
}



scProcess_SCT <- function(SeuratObj){
    ## 1. Filter off Cells
    SeuratObj <- scFilter(SeuratObj) 

    ## 2. SCT
    SeuratObj <- SCTransform(SeuratObj, assay = "RNA", vars.to.regress = c("pct.mt","pct.hb","pct.rp"), verbose = FALSE) 
    DefaultAssay(SeuratObj) <- "SCT"

    ## 3. CellCycle Effect: 仅消除非周期/周期细胞群内部的区别，c("S.Score", "G2M.Score")消除二者之间的区别
    SeuratObj <- scCellCycle(SeuratObj)
    SeuratObj <- ScaleData(SeuratObj, vars.to.regress = c("CC.Difference"), verbose = FALSE) ##, assay = "SCT"

    ## 4. 
    SeuratObj <- scRoutine(SeuratObj) ##, assay="SCT"
    return(SeuratObj)
}





