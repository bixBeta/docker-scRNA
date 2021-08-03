# DE by condition/phenoday

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
caseD1.caseD2.list = list()

for (i in 1:length(levels(Idents(all.pbmcs)))) {
    caseD1.caseD2.list[[i]] <- FindMarkers(all.pbmcs, ident.1 = "case_D1", 
                                        ident.2 = "case_D2", group.by = "condition",
                                        subset.ident = levels(Idents(all.pbmcs))[i], 
                                        only.pos = F, logfc.threshold = 0.1)
  
    names(caseD1.caseD2.list)[[i]] <- paste0("CLUSTER__", levels(Idents(all.pbmcs))[i])
  
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
controlD1.controlD2.list = list()

for (i in 1:length(levels(Idents(all.pbmcs)))) {
  controlD1.controlD2.list[[i]] <- FindMarkers(all.pbmcs, ident.1 = "control_D1", 
                                         ident.2 = "control_D2", group.by = "condition",
                                         subset.ident = levels(Idents(all.pbmcs))[i], 
                                         only.pos = F, logfc.threshold = 0.1)
  
  names(controlD1.controlD2.list)[[i]] <- paste0("CLUSTER__", levels(Idents(all.pbmcs))[i])
  
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
caseD1.controlD1.list = list()

for (i in 1:length(levels(Idents(all.pbmcs)))) {
  caseD1.controlD1.list[[i]] <- FindMarkers(all.pbmcs, ident.1 = "case_D1", 
                                            ident.2 = "control_D1", group.by = "condition",
                                            subset.ident = levels(Idents(all.pbmcs))[i], 
                                            only.pos = F, logfc.threshold = 0.1)
  
  names(caseD1.controlD1.list)[[i]] <- paste0("CLUSTER__", levels(Idents(all.pbmcs))[i])
  
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
caseD2.controlD2.list = list()

for (i in 1:length(levels(Idents(all.pbmcs)))) {
  caseD2.controlD2.list[[i]] <- FindMarkers(all.pbmcs, ident.1 = "case_D2", 
                                            ident.2 = "control_D2", group.by = "condition",
                                            subset.ident = levels(Idents(all.pbmcs))[i], 
                                            only.pos = F, logfc.threshold = 0.1)
  
  names(caseD2.controlD2.list)[[i]] <- paste0("CLUSTER__", levels(Idents(all.pbmcs))[i])
  
}

save(caseD1.caseD2.list, controlD1.controlD2.list ,caseD1.controlD1.list, caseD2.controlD2.list , file = "results/FindMarkers_list_all.pbmcs.latenet.var.Rdata"  )


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
x = list()
for (i in 1:length(caseD1.caseD2.list)) {
  x[[i]] = length(which(pluck(caseD1.caseD2.list, i, "p_val_adj") < 0.05))
  names(x)[[i]] <- names(caseD1.caseD2.list)[[i]]
}

nDE.caseD1.caseD2 = do.call("rbind",x)
colnames(nDE.caseD1.caseD2 ) = "nDE.caseD1.caseD2"

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
y = list()
for (i in 1:length(controlD1.controlD2.list)) {
  y[[i]] = length(which(pluck(controlD1.controlD2.list, i, "p_val_adj") < 0.05))
  names(y)[[i]] <- names(controlD1.controlD2.list)[[i]]
}

nDE.controlD1.controlD2 = do.call("rbind",y)
colnames(nDE.controlD1.controlD2 ) = "nDE.controlD1.controlD2"

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
z = list()
for (i in 1:length(caseD1.controlD1.list)) {
  z[[i]] = length(which(pluck(caseD1.controlD1.list, i, "p_val_adj") < 0.05))
  names(z)[[i]] <- names(caseD1.controlD1.list)[[i]]
}

nDE.caseD1.controlD1 = do.call("rbind",z)
colnames(nDE.caseD1.controlD1 ) = "nDE.caseD1.controlD1"

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
w = list()
for (i in 1:length(caseD2.controlD2.list)) {
  w[[i]] = length(which(pluck(caseD2.controlD2.list, i, "p_val_adj") < 0.05))
  names(w)[[i]] <- names(caseD2.controlD2.list)[[i]]
}

nDE.caseD2.controlD2 = do.call("rbind",w)
colnames(nDE.caseD2.controlD2 ) = "nDE.caseD2.controlD2"

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
nDE__harmony.latent.vars =  cbind(nDE.caseD1.caseD2, nDE.caseD1.controlD1, nDE.caseD2.controlD2, nDE.controlD1.controlD2)
DimPlot(all.pbmcs, raster = F, label.size = 6, label = T)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
markers = c("CD3D", "CD3E", "CD3G", "CD4", "FOXP3", "IL2RA", "CCR7", 
            "IL7R", "S100A4", "CD14", "LYZ", "CD8A", "CD8B", "CD34", "THY1", "ENG", "KIT", "PROM1" ,
            "CD19", "MS4A1", "CD19", "MS4A1", "CD79A",
            "FCGR3A", "MS4A7", "HBB", "FCGR3A", "CD14", 
            "FCGR1A", "CD68", "S100A12", "GP9", "PPBP", 
            "TYMS", "MKI67", "NCAM1", "FCGR3A", "NKG7", 
            "GNLY", "IL3RA", "CD1C", "BATF3", "THBD", "CD209" , 
            "FCER1A", "CST3", "CD19", "IGHD", "IL4R")

markers = sort(unique(markers))
library(ggplot2)
DotPlot(all.pbmcs, features = markers, dot.scale = 8) + RotatedAxis()

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
library(reshape2)
nde.gg = melt(nDE__harmony.latent.vars)
colnames(nde.gg) = c("cluster","condition", "nDE")

ggplot(nde.gg, aes(x=cluster, y=nDE , size = nDE , color = condition)) + scale_colour_viridis_d("condition") +
   facet_wrap(~condition) + ggtitle("Number of DE genes - Harmony Latent Vars") +
  geom_point(alpha=1) +  theme(axis.text.x = element_text(angle = 90, hjust = 1))


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

save(nDE__harmony.latent.vars, 
     nDE.caseD1.controlD1, 
     nDE.caseD1.caseD2, 
     nDE.caseD2.controlD2, 
     nDE.controlD1.controlD2, 
     file = "data/nDE_harmony_latent_var_scaled_objects.Rdata" )

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

