pacman::p_load(tidyverse, WGCNA, ape, pheatmap, ggpubr)


# 1. output network to cytoscape

exportNetworkToCytoscape(adjMat = adjacency_UP,
                         edgeFile = "results/network/edge_0.3_UP.tsv",
                         nodeFile = "results/network/node_0.3_UP.tsv",
                         weighted = T,
                         threshold = 0.3,
                         nodeNames = colnames(adjacency_UP),
                         altNodeNames = genes$gene_symbol)

# 2. create network in cytoscape using edge table, and process the network by removing genes that are not connected

network_nodes_UP <- read.csv("results/network/WGCNA_UP_network_nodes.csv")
network_nodes_UP <- t(exp_UP[, genes$gene_symbol %in% network_nodes_UP$name])
write.csv(network_nodes_UP, file = "results/network/network_UP_for_tSNE.csv")

# 3. run PRESTO

# 4. visualization

tsne_UP <- list()
tsne_UP$all <- readxl::read_excel("data/PRESTO_network_UP.xlsx")

#tsne all points
ggplot(tsne_UP$all, aes(x = PRESTO_X, y = PRESTO_Y)) +
  geom_point(aes(color = color), alpha = 0.6, shape = 16) + theme_bw() + 
  scale_color_manual(values = sort(names(table(tsne_UP$all$color)))) +
  ylim(-70, 75) + xlim(-85, 80) + xlab("tSNE-1") + ylab("tSNE-2") + 
  theme(axis.title.x = element_text(size = 20), axis.text.x = element_text(size = 15), 
        axis.title.y = element_text(size = 20),axis.text.y = element_text(size = 15),axis.ticks = element_blank(),
        aspect.ratio = 1,legend.position = "none", panel.background = element_rect(colour = "black",size = 1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(file = "results/network/tsne_all_UP.png", dpi = 600, width = 4, height = 4)

#tsne selected point
tsne_UP$highlight <- tsne_UP$all[tsne_UP$all$color %in% c("greenyellow", "grey60", "lightsteelblue1"),]
ggplot(tsne_UP$highlight, aes(x = PRESTO_X, y = PRESTO_Y)) +
  geom_point(aes(color = color), alpha = 0.6, shape = 16) + theme_bw() +#add colors
  scale_color_manual(values = sort(names(table(tsne_UP$highlight$color)))) +
  ylim(-70,75) + xlim(-85, 80) + xlab("tSNE-1") + ylab("tSNE-2") +
  theme(axis.title.x = element_text(size = 20), axis.text.x = element_text(size = 15), 
        axis.title.y = element_text(size = 20),axis.text.y = element_text(size = 15),axis.ticks = element_blank(),
        aspect.ratio = 1,legend.position = "none", panel.background = element_rect(colour = "black",size = 1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(file = "results/network/tsne_selected_UP.png", dpi = 600, width = 4, height = 4)




# heatmap of UP green yellow T cell cluster
anno_col <- data.frame(Plaque = factor(c(rep("SP", 16), rep("UP", 27)), levels = c("SP", "UP")))
rownames(anno_col) <- rownames(rbind(exp_SP, exp_UP))
anno_color <- list(Plaque = c(SP = "deepskyblue", UP = "coral"))
pdf("results/heatmap_UP32.pdf", width = 9, height = 18)
pheatmap(t(rbind(exp_SP, exp_UP)[, str_split(gsea_wgcna@result[5, 11], "/")[[1]]]), 
         scale = "row", cluster_cols = F, 
         clustering_distance_cols = "euclidean", 
         border_color = F,
         clustering_distance_rows = "euclidean", 
         annotation_col = anno_col, 
         cellwidth = 10, cellheight = 10,
         annotation_colors = anno_color,
         color = colorRampPalette(colors = c("blue","white","red"))(100),
         main = "UP-32 gene expression",
         show_rownames = T,
         show_colnames = F,
         gaps_col = 16,
         breaks = unique(c(seq(-2.5, 2.5, length = 100))))
dev.off()

pdf("results/heatmap_UP22.pdf", width = 9, height = 18)
pheatmap(t(rbind(exp_SP, exp_UP)[, str_split(gsea_wgcna@result[4, 11], "/")[[1]]]), 
         scale = "row", cluster_cols = F, 
         clustering_distance_cols = "euclidean", 
         border_color = F,
         clustering_distance_rows = "euclidean", 
         annotation_col = anno_col, 
         cellwidth = 10, cellheight = 10,
         annotation_colors = anno_color,
         color = colorRampPalette(colors = c("blue","white","red"))(100),
         main = "UP-22 gene expression",
         show_rownames = T,
         show_colnames = F,
         gaps_col = 16,
         breaks = unique(c(seq(-2.5, 2.5, length = 100))))
dev.off()

pdf("results/heatmap_UP11.pdf", width = 9, height = 18)
pheatmap(t(rbind(exp_SP, exp_UP)[, str_split(gsea_wgcna@result[10, 11], "/")[[1]]]), 
         scale = "row", cluster_cols = F, 
         clustering_distance_cols = "euclidean", 
         border_color = F,
         clustering_distance_rows = "euclidean", 
         annotation_col = anno_col, 
         cellwidth = 10, cellheight = 10,
         annotation_colors = anno_color,
         color = colorRampPalette(colors = c("blue","white","red"))(100),
         main = "UP-11 gene expression",
         show_rownames = T,
         show_colnames = F,
         gaps_col = 16,
         breaks = unique(c(seq(-2.5, 2.5, length = 100))))
dev.off()




# UP cluster distance
p <- list()

adjacency_UP <- readRDS("data/adjacency_UP.rds")

adj_similarity_UP32 <- list()
adj_similarity_UP32$UP22 <- adjacency_UP[moduleColors_UP == "greenyellow",moduleColors_UP=="lightsteelblue1"]
adj_similarity_UP32$UP11 <- adjacency_UP[moduleColors_UP == "greenyellow",moduleColors_UP=="grey60"]
adj_similarity_UP32$rest <- adjacency_UP[moduleColors_UP == "greenyellow",!(moduleColors_UP %in% c("greenyellow", "lightsteelblue1", "grey60"))]
adj_similarity_ggplot <- as.data.frame(cbind(c(as.vector(adj_similarity_UP32$UP22), 
                                               as.vector(adj_similarity_UP32$UP11), 
                                               as.vector(adj_similarity_UP32$rest)),
                                             c(rep("lightsteelblue1", length(adj_similarity_UP32$UP22)),
                                               rep("grey60", length(adj_similarity_UP32$UP11)),
                                               rep("rest", length(adj_similarity_UP32$rest)))))
colnames(adj_similarity_ggplot) <- c("Adjacency", "clusters")
adj_similarity_ggplot$Adjacency <- as.numeric(adj_similarity_ggplot$Adjacency)

adj_similarity_UP32$barplot <- data.frame(Adjacency = c(mean(adj_similarity_UP32$UP22), 
                                                         mean(adj_similarity_UP32$UP11), 
                                                         mean(adj_similarity_UP32$rest)),
                                           Cluster = factor(c("IFN (UP-22)","Angiogenesis (UP-11)","rest"), 
                                                            levels = c("IFN (UP-22)","Angiogenesis (UP-11)","rest")))
adj_similarity_UP32$barplot$Adjacency <- adj_similarity_UP32$barplot$Adjacency/adj_similarity_UP32$barplot$Adjacency[3]

p[[1]] <- 
  ggplot(adj_similarity_UP32$barplot, aes(x = Cluster, y = Adjacency, fill = Cluster)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.6, colour = "black") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 5.3)) + ylab("Relative cluster adjacency") +
  scale_fill_manual(values = c("lightsteelblue1", "grey60", "grey")) + ggtitle("UP network adjacency") +  
  scale_x_discrete(labels = c("UP-22 genes (IFN)", "UP-11 genes (Angiogenesis)", "rest")) +
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        plot.margin = margin(0.5, 1.8, 0.5, 0.5, "cm"), 
        axis.text.x = element_text(colour = "black", size = 12, angle = 315, vjust = 1, hjust = 0),
        axis.title.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour="black", size = 10),
        axis.title.y = element_text(size = 12),
        legend.position = "none", 
        plot.title = element_text(face = "bold", size = 12))
# ggsave(file = "results/barplot_cluster_distance_UP.pdf", dpi = 600, width = 3.5, height = 3.5)



# SP cluster distance
adjacency_SP <- readRDS("data/adjacency_SP.rds")

adj_similarity_UP32 <- list()
adj_similarity_UP32$UP22 <- adjacency_SP[moduleColors_UP == "greenyellow",moduleColors_UP=="lightsteelblue1"]
adj_similarity_UP32$UP11 <- adjacency_SP[moduleColors_UP == "greenyellow",moduleColors_UP=="grey60"]
adj_similarity_UP32$rest <- adjacency_SP[moduleColors_UP == "greenyellow",!(moduleColors_UP %in% c("greenyellow", "lightsteelblue1", "grey60"))]
adj_similarity_ggplot <- as.data.frame(cbind(c(as.vector(adj_similarity_UP32$UP22), 
                                               as.vector(adj_similarity_UP32$UP11), 
                                               as.vector(adj_similarity_UP32$rest)),
                                             c(rep("lightsteelblue1", length(adj_similarity_UP32$UP22)),
                                               rep("grey60", length(adj_similarity_UP32$UP11)),
                                               rep("rest", length(adj_similarity_UP32$rest)))))
colnames(adj_similarity_ggplot) <- c("Adjacency", "clusters")
adj_similarity_ggplot$Adjacency <- as.numeric(adj_similarity_ggplot$Adjacency)

adj_similarity_UP32$barplot <- data.frame(Adjacency = c(mean(adj_similarity_UP32$UP22), 
                                                         mean(adj_similarity_UP32$UP11), 
                                                         mean(adj_similarity_UP32$rest)),
                                           Cluster = factor(c("IFN (UP-22)","Angiogenesis (UP-11)","rest"), 
                                                            levels = c("IFN (UP-22)","Angiogenesis (UP-11)","rest")))
adj_similarity_UP32$barplot$Adjacency <- adj_similarity_UP32$barplot$Adjacency/adj_similarity_UP32$barplot$Adjacency[3]

p[[2]] <- 
  ggplot(adj_similarity_UP32$barplot, aes(x = Cluster, y = Adjacency, fill = Cluster)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.6, colour = "black") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 5.3)) + ylab("Relative cluster adjacency") +
  scale_fill_manual(values = c("lightsteelblue1", "grey60", "grey")) + ggtitle("SP network adjacency") +  
  scale_x_discrete(labels = c("UP-22 genes (IFN)", "UP-11 genes (Angiogenesis)", "rest")) +
  theme(axis.line = element_line(size = 1, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        plot.margin = margin(0.5, 1.8, 0.5, 0.5, "cm"),
        axis.text.x = element_text(colour = "black", size = 12, angle = 315, vjust = 1, hjust = 0),
        axis.title.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour="black", size = 10),
        axis.title.y = element_text(size = 12),
        legend.position = "none",
        plot.title = element_text(face = "bold", size = 12))

ggarrange(plotlist = p, ncol = 2)
ggsave(file = "results/barplot_cluster_distance.pdf", dpi = 600, width = 7, height = 4.5)



#hierarchical tree
MEs_UP_hc <- MEs_UP[, -39]
colnames(MEs_UP_hc) <- substr(colnames(MEs_UP_hc), 3, 100)
colnames(MEs_UP_hc) <- module_name$UP[colnames(MEs_UP_hc), 2]
pdf("results/hierarchical_tree.pdf", width = 9, height = 5)
plot(hclust(as.dist(1-cor(MEs_UP_hc)), method = "average"), main = "Hierarchical clustering of UP WGCNA eigengenes", xlab = "", sub = "")
dev.off()


MEs_SP_hc <- MEs_SP[, -26]
colnames(MEs_SP_hc) <- substr(colnames(MEs_SP_hc), 3, 100)
colnames(MEs_SP_hc) <- module_name$SP[colnames(MEs_SP_hc), 2]
pdf("results/hierarchical_tree_SP.pdf", width = 9, height = 5)
plot(hclust(as.dist(1-cor(MEs_SP_hc)), method = "average"), main = "Hierarchical clustering of SP WGCNA eigengenes", xlab = "", sub = "")
dev.off()







