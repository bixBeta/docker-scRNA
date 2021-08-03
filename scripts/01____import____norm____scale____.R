# import cca dataset

library(Seurat)
all.cca = readRDS("/Project____MECFS/Results/06_01____scRNAseq_118_Samples_PreProcessing/5_4____Integrated_18_batches____cluster_res0.5_.rds")
all.cca.meta = all.cca@meta.data


# import RDS
all.pbmcs <- readRDS("/Project____MECFS/Results/06_01____scRNAseq_118_Samples_PreProcessing/2_2____Merged_118_FilteredBySubset.rds")
VlnPlot(all.pbmcs, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0, group.by = "condition")

# extract metadata 
all.pbmcs.meta = all.pbmcs@meta.data


all.pbmcs <- NormalizeData(all.pbmcs, normalization.method = "LogNormalize", scale.factor = 10000)
all.pbmcs <- FindVariableFeatures(all.pbmcs, selection.method = "vst", nfeatures = 2000)

all.pbmcs <- ScaleData(all.pbmcs,vars.to.regress = c("nFeature_RNA", "percent.mt", "batch"))
all.pbmcs <- RunPCA(all.pbmcs, features = VariableFeatures(object = all.pbmcs))
#DimPlot(all.pbmcs, reduction = "pca")

#ElbowPlot(all.pbmcs)

all.pbmcs <- FindNeighbors(all.pbmcs, dims = 1:30)
all.pbmcs <- FindClusters(all.pbmcs, resolution = 0.5)

all.pbmcs <- RunUMAP(all.pbmcs, dims = 1:30)
png(filename = "figures/umap__non-integrated__group__by__batch.png", width = 1920, height = 1080)
DimPlot(all.pbmcs, group.by = "batch", raster = F)
dev.off()
saveRDS(object = all.pbmcs, file = "data/01____2_2_merged_filteredBySubset____normalized__w_latent_vars_dimRed__res__0.5.Rds")

all.cca = readRDS("/Project____MECFS/Results/06_01____scRNAseq_118_Samples_PreProcessing/5_4____Integrated_18_batches____cluster_res0.5_.rds")
all.cca.meta = all.cca@meta.data

