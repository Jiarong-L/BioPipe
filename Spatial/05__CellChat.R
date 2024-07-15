library(ggplot2)
library(Seurat)
library(CellChat)
library(ggpubr)
load('GBM4.rdata')



GBM4@meta.data$cell_type <- paste("C",GBM4$seurat_clusters,sep = "")
Idents(GBM4) <- "cell_type"    ## Previous is GBM4@meta.data$Region

## Prepare CellChatObj
data.input = Seurat::GetAssayData(GBM4, slot = "data", assay = "SCT") 
meta = data.frame(meta_labels = Idents(GBM4), row.names = names(Idents(GBM4))) 
coordinates = Seurat::GetTissueCoordinates(GBM4, scale = NULL, cols = c("imagerow", "imagecol"))[1:2] 
scale.factors = jsonlite::fromJSON(txt = "GBM4_spaceranger_out/spatial/scalefactors_json.json")
scale.factors = list(spot.diameter = 65, spot = scale.factors$spot_diameter_fullres)
CellChatObj <- createCellChat(object = data.input,
                              meta = meta, 
                              group.by = "meta_labels", 
                              datatype = "spatial", 
                              coordinates = coordinates, 
                              scale.factors = scale.factors)



## Use Secreted Signaling subset of CellChatDB
CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation")
CellChatObj@DB <- CellChatDB.use


## Preprocess
CellChatObj <- subsetData(CellChatObj)     ## always need, regardless of subsetDB() or not
future::plan("multisession", workers = 1)  ## multi-threads
CellChatObj <- identifyOverExpressedGenes(CellChatObj)
CellChatObj <- identifyOverExpressedInteractions(CellChatObj)   ## Ligand-Receptors of OverExpressedGenes


## cell-cell network
### Prob by gene = mean_CellA_Ligand * mean_CellB_Receptor
CellChatObj <- computeCommunProb(CellChatObj, 
                                 type = "truncatedMean", trim = 0.1,   ## or 'triMean' fewer but stronger interactions
                                 distance.use = TRUE, 
                                 scale.distance = 0.01)
CellChatObj <- filterCommunication(CellChatObj, min.cells = 10)
### Prob by gene networks
CellChatObj <- computeCommunProbPathway(CellChatObj)   ## for each pathway
CellChatObj <- aggregateNet(CellChatObj)   ## Aggregate pathways results



## Plot pairs and their abd, pairLR[1] is TGFB1_TGFBR1_TGFBR2
pairLR = CellChatObj@LR$LRsig$interaction_name
p1 <- spatialFeaturePlot(CellChatObj, pairLR.use = pairLR[1], point.size = 1.5, 
                   do.binary = TRUE, cutoff = 0.05, enriched.only = F, 
                   color.heatmap = "Reds", direction = 1)
p2 <- spatialFeaturePlot(CellChatObj, features = c("TGFB1"), point.size = 0.8, color.heatmap = "Reds", direction = 1)
p3 <- spatialFeaturePlot(CellChatObj, features = c("TGFBR1"), point.size = 0.8, color.heatmap = "Reds", direction = 1)
p4 <- netVisual_heatmap(CellChatObj, signaling = pairLR[1], slot.name = "net",measure = "count", color.heatmap = "Blues")




## Plot Single Pathway: aggregating all L-R pairs in that pathway
signalings = CellChatObj@netP$pathways
p5 <- netVisual_aggregate(CellChatObj, signaling = signalings[2], layout = "spatial", ## “circle”/“hierarchy”/“chord”/“spatial”
                          edge.width.max = 2, vertex.size.max = 1, 
                          alpha.image = 0.2, vertex.label.cex = 3.5)
p6 <- netVisual_heatmap(CellChatObj, signaling = signalings[2], slot.name = "netP",measure = "count", color.heatmap = "Blues") 


## Plot Aggregated : signaling = c(.............)



## Save
png("../img/05_1.png", width=900, heigh=300)
ggarrange(p1,p2,p3,ncol=3)
dev.off() 

png("../img/05_2.png", width=900, heigh=300)
p5
dev.off() 

png("../img/05_3.png", width=900, heigh=300)
p4+p6
dev.off() 



