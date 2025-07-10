pacman::p_load(tidyverse, Seurat)



# scRNA-seq GSE159677 -----------------------------------------------------

sc_GSE159677 <- readRDS("data/sc_GSE159677.rds")

axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(5, "cm"))
DimPlot(sc_GSE159677, reduction = "umap_harmony", group.by = "cell_annotation", 
        cols = sc_GSE159677_col,
        pt.size = .1, alpha = .2, raster = F, label = T, label.size = 10) + 
  theme(aspect.ratio = 1,
        plot.title = element_blank(),
        axis.line = element_line(arrow = arrow(type = "closed", length = unit(0.40, "cm"))),
        axis.title = element_text(size = 20, hjust = 0.06)) +
  guides(color = FALSE, x = axis, y = axis) + 
  xlab("UMAP-1") + ylab("UMAP-2") +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL)
ggsave("results/sc/umap_all_cells_GSE159677.pdf", width = 8, height = 7.5, dpi = 600)


marker_features <- c("CD79A", "MS4A1", # B cell
                     "JCHAIN", "MZB1", # Plasma cell
                     "CD3D", "CD3E", "CD4", "CD8A", "CD8B", # T cell
                     "IL7R", # T cell
                     "GZMA", "GZMK", # T cell
                     "NKG7", "GNLY", "KLRF1", # NK
                     "CD14", "CD68", "CLEC10A",  # Monocyte
                     "NAMPT", "CSF3R", # Neutrophil
                     "PECAM1", "VWF", # EC
                     "ACTA2", "TAGLN", "FBLN1", "LUM", # "MYH11", # SMC
                     "TPSAB1", "CPA3", "TOP2A", "MKI67") # Mast & proliferating

p <- DotPlot(object = sc_GSE159677, features = marker_features, 
             group.by = "cell_annotation", dot.min = 0.01, cols = c("white", "red"), col.min = 0) +
  theme_bw() + theme(panel.border = element_rect(colour = "black", size = 1), axis.title = element_blank(), 
                     axis.text.x = element_text(size = 12, colour = "black", angle = 315, vjust = 1, hjust = 0),
                     panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
                     legend.direction='horizontal', legend.position = "bottom", legend.justification = "left",
                     plot.title = element_blank(), axis.text.y = element_text(size = 12, colour = "black")) +
  scale_y_discrete(limits = rev(names(sc_GSE159677_col)))
p
ggsave("results/sc/marker_dotplots_GSE159677.pdf", width = 8, height = 4, units = "in", dpi = 600)


sc_GSE159677 <- AddModuleScore(sc_GSE159677, features = list(genes$gene_symbol[moduleColors_UP == "greenyellow"]))

axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(5, "cm"))
FeaturePlot(sc_GSE159677, reduction = "umap_harmony", features = "Cluster1", pt.size = 0.1) + 
  ggtitle("T cell (UP-32)") +
  theme(aspect.ratio = 1,
        plot.title = element_blank(),
        axis.line = element_line(arrow = arrow(type = "closed", length = unit(0.40, "cm"))),
        axis.title = element_text(size = 20, hjust = 0.06)) + 
  guides(color = FALSE, x = axis, y = axis) + 
  xlab("UMAP-1") + ylab("UMAP-2") + 
  scale_x_continuous(breaks = NULL) + 
  scale_y_continuous(breaks = NULL)
ggsave("results/sc/t_cell_cluster_GSE159677.pdf", width = 8, height = 7.5, dpi = 600)


VlnPlot(sc_GSE159677, features = "Cluster1", group.by = "cell_annotation", pt.size = 0) + 
  ggtitle("T cell (UP-32)") + scale_fill_manual(values = sc_GSE159677_col) + 
  scale_x_discrete(limits = names(sc_GSE159677_col)) + ylab("Module score") + 
  theme(legend.position = "none", 
        plot.title = element_text(size = 15, face = "bold"),
        axis.line = element_line(),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 14, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_blank())
ggsave("results/sc/violin_all_cell_GSE159677.pdf", width = 5.5, height = 4)


axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(4, "cm"))
DimPlot(sc_GSE159677, reduction = "umap_harmony", pt.size = 0.1, alpha = .2, raster = F, group.by = "group", cols = c(PA = "#4DBBD5BB", AC = "#E64B35BB")) + 
  ggtitle("GSE159677") + 
  theme(aspect.ratio = 1, legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 30),
        axis.line = element_line(arrow = arrow(type = "closed", length = unit(0.4, "cm"))),
        axis.title = element_text(size = 20, hjust = 0.04)) + 
  guides(x = axis, y = axis) + 
  xlab("UMAP-1") + ylab("UMAP-2") + 
  scale_x_continuous(breaks = NULL) + 
  scale_y_continuous(breaks = NULL)
ggsave("results/sc/umap_all_group_GSE159677.pdf", width = 5.5, height = 6)


axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(4, "cm"))
p <- FeaturePlot(sc_GSE159677, reduction = "umap_harmony", pt.size = 0.1,
                 features = c("PRDM1", "RUNX3", "IRF7"), raster = F, ncol = 3)
for(i in 1:length(p)){
  p[[i]] <- p[[i]] + 
    theme(aspect.ratio = 1,
          plot.title = element_text(hjust = 0.5, size = 30, face = "bold.italic"),
          axis.line = element_line(arrow = arrow(type = "closed", length = unit(0.4, "cm"))),
          axis.title = element_text(size = 20, hjust = 0.04)) +
    guides(color = FALSE, x = axis, y = axis) + 
    xlab("UMAP-1") + ylab("UMAP-2") +
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL)}
p
ggsave("results/sc/t_cell_tfs_GSE159677.pdf", width = 16.5, height = 6)



# GSE159677 T cell subset -------------------------------------------------

sc_GSE159677_t <- readRDS("data/sc_GSE159677_t.rds")
sc_GSE159677_t$Cluster1 <- sc_GSE159677@meta.data[colnames(sc_GSE159677_t), ] %>% dplyr::select(Cluster1) %>% pull

axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(5, "cm"))
DimPlot(sc_GSE159677_t, reduction = "umap_harmony", group.by = "cell_annotation", 
        cols = sc_GSE159677_tsub_col,
        pt.size = .3, alpha = .4, raster = F, label = T, label.size = 10) + 
  theme(aspect.ratio = 1,
        plot.title = element_blank(),
        axis.line = element_line(arrow = arrow(type = "closed", length = unit(0.40, "cm"))),
        axis.title = element_text(size = 20, hjust = 0.06)) +
  guides(color = FALSE, x = axis, y = axis) + 
  xlab("UMAP-1") + ylab("UMAP-2") +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL)
ggsave("results/sc/umap_t_cells_GSE159677.pdf", width = 8, height = 7.5, dpi = 600)


VlnPlot(sc_GSE159677_t, features = "Cluster1", group.by = "cell_annotation", pt.size = 0) + 
  ggtitle("T cell (UP-32)") + scale_fill_manual(values = sc_GSE159677_tsub_col) + 
  scale_x_discrete(limits = names(sc_GSE159677_tsub_col)) + ylab("Module score") + 
  theme(legend.position = "none", 
        plot.title = element_text(size = 15, face = "bold"),
        axis.line = element_line(),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 14, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_blank())
ggsave("results/sc/violin_t_cell_GSE159677.pdf", width = 4.5, height = 4)


VlnPlot(sc_GSE159677_t, features = "PRDM1", group.by = "cell_annotation", pt.size = 0) + 
  ggtitle("PRDM1") + scale_fill_manual(values = sc_GSE159677_tsub_col) + 
  scale_x_discrete(limits = names(sc_GSE159677_tsub_col)) + 
  theme(legend.position = "none", 
        plot.title = element_text(size = 15, face = "bold"),
        axis.line = element_line(),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 14, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_blank())
ggsave("results/sc/violin_PRDM1_GSE159677.pdf", width = 6, height = 4)


VlnPlot(sc_GSE159677_t, features = "PRDM1", group.by = "group", pt.size = 0) +
  scale_y_continuous(limits = c(0, 5.2), expand = c(0, 0)) +
  scale_fill_manual(values = c(PA = "#4DBBD5BB", AC = "#E64B35BB")) +
  scale_x_discrete(limits = c("PA", "AC"), labels = c("PA", "AC")) +
  theme(legend.position = "none", 
        plot.title = element_text(size = 15, face = "bold"),
        axis.line = element_line(),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 14, angle = 0, colour = "black", hjust = 0.5),
        axis.title.y = element_text(size = 14, colour = "black"),
        axis.title.x = element_blank())
ggsave("results/sc/violin_PRDM1_GSE159677_t.pdf", width = 3, height = 4)


VlnPlot(sc_GSE159677_t, features = "RUNX3", group.by = "group", pt.size = 0) +
  scale_y_continuous(limits = c(0, 5.2), expand = c(0, 0)) +
  scale_fill_manual(values = c(PA = "#4DBBD5BB", AC = "#E64B35BB")) +
  scale_x_discrete(limits = c("PA", "AC"), labels = c("PA", "AC")) +
  theme(legend.position = "none", 
        plot.title = element_text(size = 15, face = "bold"),
        axis.line = element_line(),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 14, angle = 0, colour = "black", hjust = 0.5),
        axis.title.y = element_text(size = 14, colour = "black"),
        axis.title.x = element_blank())
ggsave("results/sc/violin_RUNX3_GSE159677_t.pdf", width = 3, height = 4)

FindMarkers(sc_GSE159677_t, group.by = "group", ident.1 = "AC", ident.2 = "PA") %>% 
  rownames_to_column("gene") %>% filter(gene %in% c("RUNX3", "PRDM1"))


VlnPlot(sc_GSE159677_t, 
        features = "Cluster1", group.by = "group", pt.size = 0) +
  scale_y_continuous(limits = c(-0.1, 0.52), expand = c(0, 0)) + 
  stat_summary(fun = "median", geom = "crossbar", width = 0.5, colour = "black") +
  scale_fill_manual(values = c(PA = "#4DBBD5BB", AC = "#E64B35BB")) +
  scale_x_discrete(limits = c("PA", "AC"), labels = c("PA", "AC")) + 
  ggtitle("T cell (UP-32)") + ylab("Module score") + 
  theme(legend.position = "none", 
        plot.title = element_text(size = 15, face = "bold"),
        axis.line = element_line(),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 14, angle = 0, colour = "black", hjust = 0.5),
        axis.title.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_blank())
ggsave("results/sc/violin_UP32_GSE159677_t.pdf", width = 3, height = 4)


mean(sc_GSE159677_t$Cluster1[sc_GSE159677_t$group == "PA"])
mean(sc_GSE159677_t$Cluster1[sc_GSE159677_t$group == "AC"])
wilcox.test(sc_GSE159677_t$Cluster1[sc_GSE159677_t$group == "PA"], sc_GSE159677_t$Cluster1[sc_GSE159677_t$group == "AC"])$p.value



# scRNA-seq GSE224273 -----------------------------------------------------

sc_GSE224273 <- readRDS("data/sc_GSE224273.rds")

axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(5, "cm"))
DimPlot(sc_GSE224273, reduction = "umap_harmony", group.by = "cell_annotation", 
        cols = sc_GSE224273_col,
        pt.size = .1, alpha = .4, raster = F, label = T, label.size = 10) + 
  theme(aspect.ratio = 1,
        plot.title = element_blank(),
        axis.line = element_line(arrow = arrow(type = "closed", length = unit(0.40, "cm"))),
        axis.title = element_text(size = 20, hjust = 0.06)) +
  guides(color = FALSE, x = axis, y = axis) + 
  xlab("UMAP-1") + ylab("UMAP-2") +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL)
ggsave("results/sc/umap_all_cells_GSE224273.pdf", width = 8, height = 7.5, dpi = 600)


marker_features <- c("CD79A", "MS4A1", # B cell
                     "JCHAIN", "MZB1", # Plasma cell
                     "CD3D", "CD3E", "CD4", "CD8A", "CD8B", # T cell
                     "IL7R", # T cell
                     "GZMA", "GZMK", # T cell
                     "NKG7", "GNLY", "KLRF1", # NK
                     "CD14", "CD68", "CLEC10A",  # Monocyte
                     "NAMPT", "CSF3R", # Neutrophil
                     "PECAM1", "VWF", # EC
                     "ACTA2", "TAGLN", "FBLN1", "LUM", # "MYH11", # SMC
                     "TPSAB1", "CPA3", "TOP2A", "MKI67") # Mast & proliferating

p <- DotPlot(object = sc_GSE224273, features = marker_features, 
             group.by = "cell_annotation", dot.min = 0.01, cols = c("white", "red"), col.min = 0) +
  theme_bw() + theme(panel.border = element_rect(colour = "black", size = 1), axis.title = element_blank(), 
                     axis.text.x = element_text(size = 12, colour = "black", angle = 315, vjust = 1, hjust = 0),
                     panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
                     legend.direction='horizontal', legend.position = "bottom", legend.justification = "left",
                     plot.title = element_blank(), axis.text.y = element_text(size = 12, colour = "black")) +
  scale_y_discrete(limits = rev(names(sc_GSE224273_col)))
p
ggsave("results/sc/marker_dotplots_GSE224273.pdf", width = 8, height = 4, units = "in", dpi = 600)


sc_GSE224273 <- AddModuleScore(sc_GSE224273, features = list(genes$gene_symbol[mergedColors_UP == "greenyellow"]))

axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(5, "cm"))
FeaturePlot(sc_GSE224273, reduction = "umap_harmony", features = "Cluster1", pt.size = 0.1) + 
  ggtitle("T cell (UP-32)") +
  theme(aspect.ratio = 1,
        plot.title = element_blank(),
        axis.line = element_line(arrow = arrow(type = "closed", length = unit(0.40, "cm"))),
        axis.title = element_text(size = 20, hjust = 0.06)) + 
  guides(color = FALSE, x = axis, y = axis) + 
  xlab("UMAP-1") + ylab("UMAP-2") + 
  scale_x_continuous(breaks = NULL) + 
  scale_y_continuous(breaks = NULL)
ggsave("results/sc/t_cell_cluster_GSE224273.pdf", width = 8, height = 7.5, dpi = 600)


VlnPlot(sc_GSE224273, features = "Cluster1", group.by = "cell_annotation", pt.size = 0) + 
  ggtitle("T cell (UP-32)") + scale_fill_manual(values = sc_GSE224273_col) +
  scale_x_discrete(limits = names(sc_GSE224273_col)) + ylab("Module score") +
  theme(legend.position = "none", 
        plot.title = element_text(size = 15, face = "bold"),
        axis.line = element_line(),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 14, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_blank())
ggsave("results/sc/violin_all_cell_GSE224273.pdf", width = 5.5, height = 4)


axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(4, "cm"))
DimPlot(sc_GSE224273, reduction = "umap_harmony", pt.size = 0.1, alpha = .4,
        raster = F, group.by = "group", cols = c(asymptomatic = "#4DBBD5BB", symptomatic = "#E64B35BB")) +
  ggtitle("GSE224273") +
  theme(aspect.ratio = 1, legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 30),
        axis.line = element_line(arrow = arrow(type = "closed", length = unit(0.4, "cm"))),
        axis.title = element_text(size = 20, hjust = 0.04)) +
  guides(x = axis, y = axis) +
  xlab("UMAP-1") + ylab("UMAP-2") +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL)
ggsave("results/sc/umap_all_group_GSE224273.pdf", width = 5.5, height = 6)


axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(4, "cm"))
p <- FeaturePlot(sc_GSE224273, reduction = "umap_harmony", pt.size = 0.1,
                 features = c("PRDM1", "RUNX3", "IRF7"), raster = F, ncol = 3)
for(i in 1:length(p)){
  p[[i]] <- p[[i]] + 
    theme(aspect.ratio = 1,
          plot.title = element_text(hjust = 0.5, size = 30, face = "bold.italic"),
          axis.line = element_line(arrow = arrow(type = "closed", length = unit(0.4, "cm"))),
          axis.title = element_text(size = 20, hjust = 0.04)) +
    guides(color = FALSE, x = axis, y = axis) + 
    xlab("UMAP-1") + ylab("UMAP-2") +
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL)}
p
ggsave("results/sc/t_cell_tfs_GSE224273.pdf", width = 16.5, height = 6)



# GSE224273 T cell subset -------------------------------------------------

sc_GSE224273_t <- readRDS("data/sc_GSE224273_t.rds")
sc_GSE224273_t$Cluster1 <- sc_GSE224273@meta.data[colnames(sc_GSE224273_t), ] %>% dplyr::select(Cluster1) %>% pull

VlnPlot(sc_GSE224273_t, features = "PRDM1", group.by = "group", pt.size = 0) +
  scale_y_continuous(limits = c(0, 5.2), expand = c(0, 0)) +
  scale_fill_manual(values = c(asymptomatic = "#4DBBD5BB", symptomatic = "#E64B35BB")) +
  scale_x_discrete(limits = c("asymptomatic", "symptomatic"), labels = c("Asympt", "Sympt")) +
  theme(legend.position = "none", 
        plot.title = element_text(size = 15, face = "bold"),
        axis.line = element_line(),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 14, angle = 0, colour = "black", hjust = 0.5),
        axis.title.y = element_text(size = 14, colour = "black"),
        axis.title.x = element_blank())
ggsave("results/sc/violin_PRDM1_GSE224273_t.pdf", width = 3, height = 4)


VlnPlot(sc_GSE224273_t, features = "RUNX3", group.by = "group", pt.size = 0) +
  scale_y_continuous(limits = c(0, 5.2), expand = c(0, 0)) +
  scale_fill_manual(values = c(asymptomatic = "#4DBBD5BB", symptomatic = "#E64B35BB")) +
  scale_x_discrete(limits = c("asymptomatic", "symptomatic"), labels = c("Asympt", "Sympt")) +
  theme(legend.position = "none", 
        plot.title = element_text(size = 15, face = "bold"),
        axis.line = element_line(),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 14, angle = 0, colour = "black", hjust = 0.5),
        axis.title.y = element_text(size = 14, colour = "black"),
        axis.title.x = element_blank())
ggsave("results/sc/violin_RUNX3_GSE224273_t.pdf", width = 3, height = 4)

FindMarkers(sc_GSE224273_t, group.by = "group", ident.1 = "symptomatic", ident.2 = "asymptomatic") %>% 
  rownames_to_column("gene") %>% filter(gene %in% c("RUNX3", "PRDM1"))


VlnPlot(sc_GSE224273_t, 
        features = "Cluster1", group.by = "group", pt.size = 0) +
  scale_y_continuous(limits = c(-0.1, 0.52), expand = c(0, 0)) + 
  stat_summary(fun = "median", geom = "crossbar", width = 0.5, colour = "black") +
  scale_fill_manual(values = c(asymptomatic = "#4DBBD5BB", symptomatic = "#E64B35BB")) +
  scale_x_discrete(limits = c("asymptomatic", "symptomatic"), labels = c("Asympt", "Sympt")) + 
  ggtitle("T cell (UP-32)") + ylab("Module score") + 
  theme(legend.position = "none", 
        plot.title = element_text(size = 15, face = "bold"),
        axis.line = element_line(),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 14, angle = 0, colour = "black", hjust = 0.5),
        axis.title.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_blank())
ggsave("results/sc/violin_UP32_GSE224273_t.pdf", width = 3, height = 4)


mean(sc_GSE224273_t$Cluster1[sc_GSE224273_t$group == "asymptomatic"])
mean(sc_GSE224273_t$Cluster1[sc_GSE224273_t$group == "symptomatic"])
wilcox.test(sc_GSE224273_t$Cluster1[sc_GSE224273_t$group == "asymptomatic"], sc_GSE224273_t$Cluster1[sc_GSE224273_t$group == "symptomatic"])$p.value

