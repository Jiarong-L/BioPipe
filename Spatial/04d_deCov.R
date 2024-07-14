library(ggplot2)
library(Seurat)
library(dplyr)
library(tidyverse)
library(CARD)
library(MuSiC)
load('sc.rdata')
load('GBM4.rdata')


sc_count <- scRNA@assays$RNA$counts
spatial_count <- GBM4@assays$Spatial$counts
spatial_location <- GetTissueCoordinates(GBM4)[,1:2]    ## colnames() = c("x","y")
sc_meta <- scRNA@meta.data %>% 
  rownames_to_column("cellID") %>%
  dplyr::select(cellID,orig.ident,cell_type) %>% 
  mutate(CB = cellID) %>% 
  column_to_rownames("CB")

# > head(sc_meta)
#                                cellID    orig.ident       cell_type
# AAACCTGAGTCAAGGC_1 AAACCTGAGTCAAGGC_1 SeuratProject       Astrocyte
# AAACCTGTCAGGCAAG_1 AAACCTGTCAGGCAAG_1 SeuratProject         Myeloid
# AAACCTGTCCTGCCAT_1 AAACCTGTCCTGCCAT_1 SeuratProject Neoplastic_cell
# AAACCTGTCGGTCCGA_1 AAACCTGTCGGTCCGA_1 SeuratProject       Astrocyte
# AAACCTGTCTATGTGG_1 AAACCTGTCTATGTGG_1 SeuratProject         Unknown
# AAACGGGAGCTCAACT_1 AAACGGGAGCTCAACT_1 SeuratProject         Unknown



CARD_obj = createCARDObject( 
  sc_count = sc_count, 
  sc_meta = sc_meta, 
  spatial_count = spatial_count, 
  spatial_location = spatial_location, 
  ct.varname = "cell_type", 
  ct.select = unique(sc_meta$cell_type), ## All cell_type
  sample.varname = "orig.ident")   ## sample name: we have only one sample here, so didn't set up

CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)   ## deConvolution!!











save(scRNA,file = 'scRNA.rdata')

