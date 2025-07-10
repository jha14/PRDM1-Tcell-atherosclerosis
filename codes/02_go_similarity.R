pacman::p_load(tidyverse, clusterProfiler, org.Hs.eg.db, GOSemSim, enrichplot, viridis, ggpubr, VennDiagram, ggrepel)


# GSOA -----------------------------------------------------------

GO_result_SP <- list()
GO_result_UP <- list()
for(i in unique(mergedColors_SP))
  GO_result_SP[[i]] <- enrich_go_cp(i, mergedColors_SP)
  
for(i in unique(mergedColors_UP))
  GO_result_UP[[i]] <- enrich_go_cp(i, mergedColors_UP)


GO_result_SP <- 
  bind_rows(GO_result_SP) %>% 
  filter(module != "grey")

GO_result_UP <- 
  bind_rows(GO_result_UP) %>% 
  filter(module != "grey")

GO_result_SP$module <- module_name$SP[GO_result_SP$module, ]$cluster
GO_result_UP$module <- module_name$UP[GO_result_UP$module, ]$cluster

GO_result_SP <- GO_result_SP[order(as.numeric(str_split(GO_result_SP$module, pattern = "-", simplify = T)[, 2])), ]
GO_result_UP <- GO_result_UP[order(as.numeric(str_split(GO_result_UP$module, pattern = "-", simplify = T)[, 2])), ]

saveRDS(GO_result_SP, "data/GO_result_SP.rds")
saveRDS(GO_result_UP, "data/GO_result_UP.rds")

GO_result_SP <- readRDS("data/GO_result_SP.rds")
GO_result_UP <- readRDS("data/GO_result_UP.rds")


ggarrange(plotlist = list(plot_GO_clusters(GO_result_UP, "UP clusters"), 
                          plot_GO_clusters(GO_result_SP, "SP clusters")), ncol = 1, align = "hv")
ggsave(file = "results/CES_UP_SP.png", width = 4.5, height = 8)



# GSEA UP clusters -----------------------------------------------------------

protein_atlas <- 
  read_tsv("data/proteinatlas.tsv") %>% 
  dplyr::select(gene_symbol = Gene, immune_specificity = "RNA blood cell specificity", immune_distribution = "RNA blood cell distribution", 
                immune_specificity_score = "RNA blood cell specificity score", immune_specificity_nTPM = "RNA blood cell specific nTPM",
                sc_specificity = "RNA single cell type specificity", sc_distribution = "RNA single cell type distribution", 
                sc_specificity_score = "RNA single cell type specificity score", sc_specificity_nTPM = "RNA single cell type specific nTPM")

protein_atlas_neut <- 
  protein_atlas %>% 
  filter(grepl("neutrophil", immune_specificity_nTPM)) %>% 
  filter(grepl("Immune cell enriched", immune_specificity))

protein_atlas_mf <- 
  protein_atlas %>% 
  filter(grepl("Macrophage", sc_specificity_nTPM)) %>% 
  filter(grepl("Cell type enhanced", sc_specificity))

protein_atlas_t <- 
  protein_atlas %>% 
  filter(grepl("T-cells", sc_specificity_nTPM)) %>% 
  filter(!grepl("Group enriched", sc_specificity))


gsea_wgcna_terms_to_gene <- 
  data.frame(cluster = module_name$UP[unique(moduleColors_UP)[1], 2], 
             gene = genes$gene_symbol[moduleColors_UP == unique(moduleColors_UP)[1]])

for(i in unique(moduleColors_UP)[-1])
  gsea_wgcna_terms_to_gene <- rbind(gsea_wgcna_terms_to_gene, 
                                    data.frame(cluster = module_name$UP[i, 2], 
                                               gene = genes$gene_symbol[moduleColors_UP == i]))

gsea_wgcna_terms_to_gene <- 
  bind_rows(gsea_wgcna_terms_to_gene, 
            data.frame(cluster = "HPA Neutrophil signature", gene = protein_atlas_neut$gene_symbol),
            data.frame(cluster = "HPA Macrophage signature", gene = protein_atlas_mf$gene_symbol),
            data.frame(cluster = "HPA T cell signature", gene = protein_atlas_t$gene_symbol))


gsea_wgcna_degs <- genes$log2fc
names(gsea_wgcna_degs) <- genes$gene_symbol
gsea_wgcna_degs <- gsea_wgcna_degs[order(gsea_wgcna_degs, decreasing = T)]
gsea_wgcna <- GSEA(geneList = gsea_wgcna_degs, TERM2GENE = gsea_wgcna_terms_to_gene, eps = 0, pvalueCutoff = 1, maxGSSize = 500)
View(gsea_wgcna@result)

saveRDS(gsea_wgcna, file = "data/gsea_wgcna.rds")
gsea_wgcna <- readRDS("data/gsea_wgcna.rds")

gseaplot2(gsea_wgcna, geneSetID = 1, title = gsea_wgcna$Description[1])
ggsave(file = "results/GSEA_HPA_mf.pdf", width = 3.5, height = 4)
gseaplot2(gsea_wgcna, geneSetID = 2, title = gsea_wgcna$Description[2])
ggsave(file = "results/GSEA_HPA_t.pdf", width = 3.5, height = 4)
gseaplot2(gsea_wgcna, geneSetID = 12, title = gsea_wgcna$Description[12])
ggsave(file = "results/GSEA_HPA_neut.pdf", width = 3.5, height = 4)

gseaplot2(gsea_wgcna, geneSetID = 4, title = gsea_wgcna$Description[4])
ggsave(file = "results/GSEA_UP_22.pdf", width = 3.5, height = 4)
gseaplot2(gsea_wgcna, geneSetID = 5, title = gsea_wgcna$Description[5])
ggsave(file = "results/GSEA_UP_32.pdf", width = 3.5, height = 4)
gseaplot2(gsea_wgcna, geneSetID = 11, title = gsea_wgcna$Description[11])
ggsave(file = "results/GSEA_UP_11.pdf", width = 3.5, height = 4)

gsea_wgcna@result$core_enrichment[4] %>% str_split("/")
gsea_wgcna@result$core_enrichment[5] %>% str_split("/")
gsea_wgcna@result$core_enrichment[11] %>% str_split("/")



gsea_wgcna@result %>% 
  dplyr::select(ID, NES, p.adjust) %>% 
  filter(str_detect(ID, "UP")) %>% 
  arrange(-NES) %>% 
  mutate(color = ifelse(NES > 0, "up", "down"),
         logp = -log10(p.adjust),
         cluster = 1:31,
         ID_show = ifelse(p.adjust < 0.05, ID, "")) %>%
  mutate(color = factor(color, levels = c("up", "down")),
         cluster = as.character(cluster)) %>%
  ggplot(aes(x = cluster, y = NES)) + 
  geom_point(aes(color = color, size = logp, alpha = logp), shape = 16) +
  geom_text_repel(aes(label = ID_show)) +
  scale_x_discrete(limits = as.character(c(1:31))) +
  scale_color_manual(values = c("up" = "#B2182B", "down" = "#2166AC")) +
  theme_bw() + xlab("UP clusters") + ylab("Normalized enrichment score (NES)") +
  theme(panel.grid = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_text(color = "black"), 
        legend.position = "bottom")
ggsave("results/GSEA_summary.pdf", height = 4, width = 3.5)



# GO similarity -----------------------------------------------------------

hsGO <- godata('org.Hs.eg.db', ont = "BP")

go_sim <- matrix(0, nrow = length(unique(GO_result_SP$module)), ncol = length(unique(GO_result_UP$module)))
rownames(go_sim) <- unique(GO_result_SP$module)
colnames(go_sim) <- unique(GO_result_UP$module)

for(i in rownames(go_sim))
  for(j in colnames(go_sim))
    go_sim[i,j] <- mgoSim(GO_result_SP$ID[GO_result_SP$module == i], GO_result_UP$ID[GO_result_UP$module == j], hsGO)



ces_SP <- c()
for(i in unique(GO_result_SP$module))
  ces_SP <- c(ces_SP, -log10(GO_result_SP$p.adjust[GO_result_SP$module == i][1]))
names(ces_SP) <- unique(GO_result_SP$module)
ces_SP <- sort(ces_SP, decreasing = T)

ces_UP <- c()
for(i in unique(GO_result_UP$module))
  ces_UP <- c(ces_UP, -log10(GO_result_UP$p.adjust[GO_result_UP$module == i][1]))
names(ces_UP) <- unique(GO_result_UP$module)
ces_UP <- sort(ces_UP, decreasing = T)

go_sim_top10 <- go_sim[names(ces_SP)[1:10], names(ces_UP)[1:10]]


go_sim_top10 <- readRDS("data/GO_sim_top10.rds")
pheatmap::pheatmap(t(go_sim_top10), cluster_rows = F, cluster_cols = F, color = viridis(100), cellwidth = 15, cellheight = 15,
                   border_color = F, angle_col = 270, filename = "results/go_sim_top10.pdf", width = 3.5, height = 3)



# venn diagram ------------------------------------------------------------

# module_name$SP
# module_name$UP

p <- venn.diagram(list(UP = genes$gene_symbol[moduleColors_UP == "greenyellow"],
                       SP = genes$gene_symbol[moduleColors_SP == "darkorange"]),
                  category.names = c("UP-32", "SP-7"),  fill = c("greenyellow", "darkorange"),
                  filename = NULL,
                  cex = 1, fontface = "bold", fontfamily = "sans",
                  cat.cex = 1, cat.fontface = "bold", cat.default.pos = "outer",
                  cat.pos = c(-160, 160), cat.dist = c(0.055, 0.055),
                  cat.fontfamily = "sans")
pdf("results/venn1.pdf", height = 2, width = 2)
grid.draw(p)
dev.off()

p <- venn.diagram(list(UP = genes$gene_symbol[moduleColors_UP == "darkorange"],
                       SP = genes$gene_symbol[moduleColors_SP == "darkorange"]),
                  category.names = c("UP-5", "SP-7"), fill = c(scales::muted("darkorange"), "darkorange"),
                  filename = NULL,
                  cex = 1, fontface = "bold", fontfamily = "sans",
                  cat.cex = 1, cat.fontface = "bold", cat.default.pos = "outer",
                  cat.pos = c(-160, 160), cat.dist = c(0.055, 0.055),
                  cat.fontfamily = "sans")
pdf("results/venn2.pdf", height = 2, width = 2)
grid.draw(p)
dev.off()

p <- venn.diagram(list(UP = genes$gene_symbol[moduleColors_UP == "magenta"],
                       SP = genes$gene_symbol[moduleColors_SP == "antiquewhite2"]),
                  category.names = c("UP-34", "SP-14"), fill = c("magenta", "antiquewhite2"),
                  filename = NULL,
                  cex = 1, fontface = "bold", fontfamily = "sans",
                  cat.cex = 1, cat.fontface = "bold", cat.default.pos = "outer",
                  cat.pos = c(-160, 160), cat.dist = c(0.055, 0.055),
                  cat.fontfamily = "sans")
pdf("results/venn3.pdf", height = 2, width = 2)
grid.draw(p)
dev.off()

p <- venn.diagram(list(UP = genes$gene_symbol[moduleColors_UP == "grey60"],
                       SP = genes$gene_symbol[moduleColors_SP == "salmon"]),
                  category.names = c("UP-11", "SP-18"), fill = c("grey60", "salmon"),
                  filename = NULL,
                  cex = 1, fontface = "bold", fontfamily = "sans",
                  cat.cex = 1, cat.fontface = "bold", cat.default.pos = "outer",
                  cat.pos = c(-160, 160), cat.dist = c(0.055, 0.055),
                  cat.fontfamily = "sans")
pdf("results/venn4.pdf", height = 2, width = 2)
grid.draw(p)
dev.off()


hyper_test(genes$gene_symbol[moduleColors_UP == "greenyellow"], genes$gene_symbol[moduleColors_SP == "darkorange"])
hyper_test(genes$gene_symbol[moduleColors_UP == "darkorange"], genes$gene_symbol[moduleColors_SP == "darkorange"])
hyper_test(genes$gene_symbol[moduleColors_UP == "magenta"], genes$gene_symbol[moduleColors_SP == "antiquewhite2"])
hyper_test(genes$gene_symbol[moduleColors_UP == "grey60"], genes$gene_symbol[moduleColors_SP == "salmon"])




jaccard_index <- matrix(0, nrow = length(module_name$SP$cluster), ncol = length(module_name$UP$cluster))
hyper_p <- matrix(0, nrow = length(module_name$SP$cluster), ncol = length(module_name$UP$cluster))
rownames(jaccard_index) <- module_name$SP$cluster
colnames(jaccard_index) <- module_name$UP$cluster
rownames(hyper_p) <- module_name$SP$cluster
colnames(hyper_p) <- module_name$UP$cluster

for(i in 1:nrow(jaccard_index))
  for(j in 1:ncol(jaccard_index)){
    SP_genes <- genes$gene_symbol[moduleColors_SP == module_name$SP$module_color[i]]
    UP_genes <- genes$gene_symbol[moduleColors_UP == module_name$UP$module_color[j]]
    jaccard_index[i, j] <- length(unique(intersect(SP_genes, UP_genes)))/length(unique(union(SP_genes, UP_genes)))
    hyper_p[i, j] <- hyper_test(SP_genes, UP_genes)}

hyper_padj <- p.adjust(hyper_p, method = "BH") %>% matrix(nrow = nrow(hyper_p), ncol = ncol(hyper_p))
hyper_p[hyper_padj < 1e-10] <- "*"
hyper_p[hyper_padj < 1e-20] <- "**"
hyper_p[hyper_padj > 1e-10] <- ""

pheatmap::pheatmap(jaccard_index, cluster_rows = F, cluster_cols = F, color = viridis(100), 
                   cellwidth = 12, cellheight = 12,
                   display_numbers = hyper_p, fontsize = 10, fontsize_number = 15, 
                   border_color = F, angle_col = 270, number_color = "black",
                   labels_col = paste0(module_name$UP$cluster, " [", module_name$UP$size, "]"), 
                   labels_row = paste0(module_name$SP$cluster, " [", module_name$SP$size, "]"), 
                   filename = "results/jaccard_index.pdf", width = 8.3, height = 5.5)


writexl::write_xlsx(jaccard_index %>% as.data.frame, path = "jaccard_index.xlsx")
writexl::write_xlsx(hyper_padj %>% as.data.frame, path = "hyper_padj.xlsx")


