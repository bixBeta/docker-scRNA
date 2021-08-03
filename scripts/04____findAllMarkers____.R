all.pbmcs <- readRDS("data/seurat__harmonized__default__findNeighbors_dim_30__clust_res_0.5.Rds")
all.pbmcs.meta = all.pbmcs@meta.data

all.pbmcs.markers <- FindAllMarkers(all.pbmcs, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

write.csv(all.pbmcs.markers, "results/all.pbmcs.markers.csv")

FeaturePlot(all.pbmcs, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                               "CD8A"))
