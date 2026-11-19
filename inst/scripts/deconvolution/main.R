# library(VoltRon)
library(Seurat)
library(spacexr)
library(testthat)

###
# Visium niche clustering ####
###

####
## Import data ####
####

# import Visium data
MBrain_Sec <- importVisium("../../../../../data/10X_Visium_Mouse_Brain/Sagittal_Anterior/Section1/",
                         sample_name = "Anterior1")

####
## test old Seurat object ####
####

allen_reference <- readRDS("../../../../../data/10X_Visium_Mouse_Brain/allen_cortex_analyzed.rds")

# get assay data 
expect_error(
  tmp <- GetAssayData(allen_reference, assay = "RNA", slot = "counts") 
)
tmp <- GetAssayData(allen_reference, 
                    assay = "RNA", 
                    layer = "counts")

# visualize
Idents(allen_reference) <- "subclass"
gsubclass <- DimPlot(allen_reference, reduction = "umap", label = T) + NoLegend()
Idents(allen_reference) <- "class"
gclass <- DimPlot(allen_reference, reduction = "umap", label = T) + NoLegend()
gsubclass | gclass

# deconv
MBrain_Sec2 <- getDeconvolution(MBrain_Sec, sc.object = allen_reference, sc.cluster = "subclass", max_cores = 6)

####
## test new Seurat object ####
####

# get seurat 5 data
allen_reference_new <- CreateSeuratObject(counts = GetAssayData(allen_reference, layer = "counts"), 
                                          meta.data = allen_reference@meta.data)

# get assay data 
expect_error(
  tmp <- GetAssayData(allen_reference_new, assay = "RNA", slot = "counts") 
)
tmp <- GetAssayData(allen_reference_new, 
                    assay = "RNA", 
                    layer = "counts")

# deconv
MBrain_Sec2 <- getDeconvolution(MBrain_Sec, sc.object = allen_reference_new, sc.cluster = "subclass", max_cores = 6)