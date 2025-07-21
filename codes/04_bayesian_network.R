pacman::p_load(tidyverse, bnlearn, igraph, RCy3)



# bnlearn UP eigengenes --------------------------------------------------

MEs_UP_bn <- MEs_UP
colnames(MEs_UP_bn) <- str_replace(colnames(MEs_UP_bn), pattern = "ME", replacement = "")
colnames(MEs_UP_bn) <- module_name$UP[colnames(MEs_UP_bn),]$cluster

rownames(UP_cluster_functions) <- UP_cluster_functions$clusters
key_cluster_UP <- table(GO_result_UP$module)
for(i in names(key_cluster_UP))
  key_cluster_UP[i] <- -log10(min(GO_result_UP$p.adjust[GO_result_UP$module == i]))
MEs_UP_bn <- MEs_UP_bn[, names(key_cluster_UP)[key_cluster_UP > 3]]
dim(MEs_UP_bn)



# bnlearn UP -------------------------------------------------------------

bn_UP_bootstrap <- boot.strength(MEs_UP_bn, R = 1000, algorithm = "hc")
saveRDS(bn_UP_bootstrap, "data/bn_UP_bootstrap.rds")

bn_UP_bootstrap <- readRDS("data/bn_UP_bootstrap.rds")
avg_diff <- averaged.network(bayesian.network.bootstrap)

# bn_UP_edges <- bn_UP_bootstrap[bn_UP_bootstrap$direction > 0.5 & bn_UP_bootstrap$strength > avg_diff$learning$args$threshold, ]
bn_UP_edges <- bn_UP_bootstrap[bn_UP_bootstrap$direction > 0.5 & bn_UP_bootstrap$strength > 0.504, ]
bn_UP_edges$color <- "grey60"
bn_UP_edges[bn_UP_edges$from == "UP-32" & bn_UP_edges$to == "UP-11", 5] = "red"
bn_UP_edges[bn_UP_edges$from == "UP-32" & bn_UP_edges$to == "UP-22", 5] = "red"

bn_UP_matrix <- matrix(0, length(table(c(bn_UP_edges[, 1], bn_UP_edges[, 2]))),
                       length(table(c(bn_UP_edges[, 1], bn_UP_edges[, 2]))))
colnames(bn_UP_matrix) <- names(table(c(bn_UP_edges[, 1], bn_UP_edges[, 2])))
rownames(bn_UP_matrix) <- names(table(c(bn_UP_edges[, 1], bn_UP_edges[, 2])))

for(i in 1:nrow(bn_UP_edges))
  bn_UP_matrix[bn_UP_edges[i, 1],bn_UP_edges[i, 2]] <- bn_UP_edges[i, 3]
bn_UP_igraph <- graph_from_adjacency_matrix(bn_UP_matrix, weighted = T)


plot(bn_UP_igraph)

cytoscapePing()
createNetworkFromIgraph(bn_UP_igraph, new.title = 'Bayesian network')



