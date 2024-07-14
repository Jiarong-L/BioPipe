library(Seurat)
library(dplyr)
library(stringr)



## Load scRNA Data with cell_types
scRNAlist <- list()
for(i in 1:9){
    dir = paste("scRNA/GSM411953",i,sep="")
    counts <- Read10X(data.dir = dir,gene.column = 1)
    scRNAlist[[i]] <- CreateSeuratObject(counts,  min.cells=3, min.features = 200)
    scRNAlist[[i]]$Source <- paste("Sample",i,sep="")

    ## Replace QC here！！！
    scRNAlist[[i]][["mt_percent"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-")
    scRNAlist[[i]] <- subset(scRNAlist[[i]], subset = nFeature_RNA > 200 & nFeature_RNA  < 5000 & mt_percent < 10) 
}

## Merge: stack df, 9 layers present
## 10370 features across 1891 samples(cells) ==> 16291 features across 19318 samples(cells)
scRNA <- merge(scRNAlist[[1]], y = scRNAlist[2:9])


##### QC --- subset() may cause error if one of the samples returns NULL, replace QC above!!
#### scRNA[["mt_percent"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
#### scRNA <- subset(scRNA, subset = nFeature_RNA > 200 & nFeature_RNA  < 500 & mt_percent < 10) 


## SCTransform: Normalize & HVGs & Scale, then pca
scRNA <- NormalizeData(scRNA) %>% 
  FindVariableFeatures(selection.method = "vst",nfeatures = 3000) %>% 
  ScaleData() %>% 
  RunPCA(npcs = 30, verbose = T) 


## Load scRNA cell_types
myfunc <- function(x,sep,col){return(  str_split(x,sep)[[1]][col]     )}
cell_types <- read.table(gzfile("scRNA/GSE138794_scRNA_Seq_cell_types.txt.gz"))
cell_types_list <- cell_types$V2
names(cell_types_list) <- lapply(cell_types$V1, myfunc, sep='_', col=2)

anno <- cell_types_list[unlist(lapply(rownames(scRNA[[]]), myfunc, sep='_',col=1))] 
names(anno) <- rownames(scRNA[[]])
for(i in 1:length(anno)){
  if (is.na(anno[i])){
    anno[i]<- 'Unknown'      ## no, perhaps I should remove these NA cells? 
    }  
}
scRNA$cell_type <- anno

save(scRNA,file = 'sc.rdata')