setwd("C:/Users/vitkauskm/OneDrive - A STAR/Hackathon_spatial")
library(readxl)
library(Seurat)
library(ggplot2)
library(readr)
library(stringr)

# load unfiltered feature- barcode matrix # and x,y coordinates
xenium.obj <- LoadXenium(
  "outs",
  fov = "fov"
)

# load cell-type annotations per barcode
cell_type <- read_xlsx(path = "Cell_Barcode_Type_Matrices.xlsx", sheet = 4)

index <- match(cell_type$Barcode, colnames(xenium.obj))
xenium.obj$cell_type <- "test"
xenium.obj$cell_type <- cell_type$Cluster[index]

# filter out cells with no spots
xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 0)
# VlnPlot(xenium.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)
# visualize the tissue map
# ImageDimPlot(xenium.obj, fov = "fov", molecules = c("LUM", "AQP1", "MKI67", "MYH11"), nmols = 20000)
# Idents(xenium.obj) <- xenium.obj$cell_type
# ImageDimPlot(xenium.obj, fov = "fov", split.by = "cell_type")
# figure out size of grid
g_n <- 250
print("grid of size (microns): ")
print(g_n)

# get metadata per cell
coord <- GetTissueCoordinates(xenium.obj)
coord$cell_type <- xenium.obj$cell_type
coord$gridGroup <- interaction(cut(coord$x, breaks=seq(0, max(coord$x) + g_n,by=g_n)),
                               cut(coord$y, breaks=seq(0, max(coord$y) + g_n, by=g_n)), sep=";")

# ggplot(data = coord, aes(x=x, y=y, color = gridGroup)) +
#   geom_hex() + guides(fill="none")
coord$Categorize <- paste0(coord$gridGroup, "_", coord$cell_type)
freq_table <- table(coord$cell_type, coord$gridGroup)
normalize_columns <- function(data) {
  # Calculate the sum of each column
  col_sums <- colSums(data)
  
  # Divide each element in the column by its respective column sum
  normalized_data <- t(t(data) / col_sums)
  
  return(normalized_data)
}
normalized_data <- normalize_columns(freq_table)