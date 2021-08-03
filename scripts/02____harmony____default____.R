library(harmony)
library(Seurat)
library(dplyr)

# all.pbmcs = readRDS(file = "data/01____2_2_merged_filteredBySubset____normalized__w_latent_vars_dimRed__res__0.5.Rds")
# DimPlot(all.pbmcs)
# all.pbmcs.meta = all.pbmcs@meta.data

all.pbmcs = RunHarmony(object = all.pbmcs, group.by.vars = "batch", plot_convergence = T)

Reductions(all.pbmcs)

saveRDS(all.pbmcs, "data/02____seurat__harmonized__default.Rds")


all.pbmcs <- RunUMAP(all.pbmcs, reduction = "harmony", dims = 1:30)

all.pbmcs <- FindNeighbors(all.pbmcs, reduction = "harmony", dims = 1:30) %>% FindClusters(res = 0.5)
DimPlot(all.pbmcs, group.by = c("batch", "ident"), ncol = 2, raster = F)

saveRDS(all.pbmcs, "data/seurat__harmonized__default__findNeighbors_dim_30__clust_res_0.5.Rds")
