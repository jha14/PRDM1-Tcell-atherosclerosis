pacman::p_load(GENIE3, minet, WGCNA, VennDiagram)


tfs_all <- read.csv(file = "data/human_transcription_factors.csv")
tfs_all <- tfs_all[, -c(8:11)]
colnames(tfs_all)[1] <- "ID"
tfs_all <- tfs_all[tfs_all$isTF == "Yes", ]
tfs_all <- tfs_all[tfs_all$Name %in% genes$gene_symbol, ]


set.seed(123) # For reproducibility of results
genie3_mtx_tfs <- GENIE3(t(rbind(exp_SP, exp_UP)), regulators = tfs_all$Name, nCores = 12)

aracne_mtx <- minet(rbind(exp_SP, exp_UP), method = "aracne")
aracne_mtx_tfs <- aracne_mtx[tfs_all$Name, ]


genie3_mtx_tfs <- readRDS("data/GENIE3.rds")
aracne_mtx_tfs <- readRDS("data/ARACNe.rds")


scaleFreeFitIndex(rowSums(genie3_mtx_tfs))
scaleFreePlot(rowSums(genie3_mtx_tfs))
scaleFreeFitIndex(rowSums(aracne_mtx_tfs))
scaleFreePlot(rowSums(aracne_mtx_tfs))


top_tfs <- 100
tfs_list_T_cell <- list()
tfs_list_T_cell$genie3_tfs = names(sort(rowSums(genie3_mtx_tfs[, moduleColors_UP == "greenyellow"]), decreasing = T))[1: top_tfs]
tfs_list_T_cell$aracne_tfs = names(sort(rowSums(aracne_mtx_tfs[, moduleColors_UP == "greenyellow"]), decreasing = T))[1: top_tfs]

venn.diagram(list(tfs_list_T_cell$genie3_tfs, tfs_list_T_cell$aracne_tfs), 
             category.names = c("GENIE3", "ARACNe"), 
             imagetype = "png", height = 1200 , width = 1200, resolution = 600, compression = "lzw",
             fill = c("blue", "red"),
             cex = 1, fontface = "bold", fontfamily = "sans",
             cat.cex = 1, cat.fontface = "bold", cat.default.pos = "outer",
             cat.pos = c(-20, 20), cat.dist = c(0.055, 0.055),
             cat.fontfamily = "sans",
             filename = "results/TFs_overlapping.png")


