library(viridis)
library(reshape2)
library(ggplot2)

DimPlot(all.pbmcs)

n_cells_harmony <- FetchData(all.pbmcs, 
                             vars = c("ident", "orig.ident")) %>%
  dplyr::count(ident, orig.ident) %>%
  tidyr::spread(ident, n)


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sum 100 to patient Harmony nclust 

n_cells.matrix = n_cells_harmony %>% select(-1)
n_cells.matrix[is.na(n_cells.matrix)] <- 0

n_cells.matrix = n_cells.matrix / rowSums(n_cells.matrix) * 100
n_cells.matrix$sample = n_cells_harmony$orig.ident

n_cells.matrix.joined = left_join(n_cells.matrix, all.pbmcs.meta, by = c("sample" = "orig.ident")) 

gg3.join = n_cells.matrix.joined %>% select(1:22,batch)
gg3.melt = reshape2::melt(gg3.join)


png("figures/latentVar__all.pbmcs.defaultHarmony.SUM-by-patient.png", width = 9000, height = 10000, res = 200)
ggplot(gg3.melt,
       aes(fill=batch, y=value, x=sample, label = batch)) + 
  geom_bar(position="dodge", stat="identity") + theme(legend.position = "none", 
                                                      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_viridis(discrete = T) +
  #facet_wrap(~variable + batch , ncol = 18, scales = "free_y") 
  facet_grid(variable ~ batch , scales = "free") + theme(axis.line=element_line())
#facet_rep_wrap(~ interaction(variable, batch), scales='free', repeat.tick.labels = c('top','left'))
dev.off()

