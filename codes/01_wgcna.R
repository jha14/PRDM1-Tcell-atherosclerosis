pacman::p_load(tidyverse, WGCNA)


exp_SP <- readRDS("data/exp_SP.rds")
exp_UP <- readRDS("data/exp_UP.rds")
genes <- readRDS("data/genes.rds")

powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))



# SP WGCNA ---------------------------------------------------------------

sft_SP <- pickSoftThreshold(exp_SP, powerVector = powers, verbose = 5, blockSize = 15000)
adjacency_SP <- adjacency(exp_SP, power = 6)
TOM_SP <- TOMsimilarity(adjacency_SP)
dissTOM_SP <- 1-TOM_SP
saveRDS(adjacency_SP, "data/adjacency_SP.rds")
saveRDS(TOM_SP, "data/TOM_SP.rds")

geneTree_SP <- hclust(as.dist(dissTOM_SP), method = "average")
plot(geneTree_SP, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
minModuleSize <- 30
dynamicMods_SP <- cutreeDynamic(dendro = geneTree_SP, distM = dissTOM_SP, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize)
dynamicColors_SP <- labels2colors(dynamicMods_SP)
plotDendroAndColors(geneTree_SP, dynamicColors_SP, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")


MEList <- moduleEigengenes(exp_SP, colors = dynamicColors_SP)
MEs <- MEList$eigengenes
MEDiss <- 1-cor(MEs)
METree <- hclust(as.dist(MEDiss), method = "average")
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
MEDissThres <- 0.25
abline(h = MEDissThres, col = "red")
rm(MEList, MEs, MEDiss, METree)

merge_SP <- mergeCloseModules(exp_SP, dynamicColors_SP, cutHeight = MEDissThres, verbose = 3)
mergedColors_SP <- merge_SP$colors
mergedMEs_SP <- merge_SP$newMEs

plotDendroAndColors(geneTree_SP, cbind(dynamicColors_SP, mergedColors_SP), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)

moduleColors_SP <- mergedColors_SP
colorOrder <- c("grey", standardColors(50))
moduleLabels_SP <- match(moduleColors_SP, colorOrder)-1
MEs_SP <- mergedMEs_SP



# UP WGCNA ---------------------------------------------------------------

sft_UP <- pickSoftThreshold(exp_UP, powerVector = powers, verbose = 5, blockSize = 15000)
adjacency_UP <- adjacency(exp_UP, power = 6)
TOM_UP <- TOMsimilarity(adjacency_UP)
dissTOM_UP <- 1-TOM_UP
saveRDS(adjacency_UP, "data/adjacency_UP.rds")
saveRDS(TOM_UP, "data/TOM_UP.rds")

geneTree_UP <- hclust(as.dist(dissTOM_UP), method = "average")
plot(geneTree_UP, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
minModuleSize <- 30
dynamicMods_UP <- cutreeDynamic(dendro = geneTree_UP, distM = dissTOM_UP, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize)
dynamicColors_UP <- labels2colors(dynamicMods_UP)
plotDendroAndColors(geneTree_UP, dynamicColors_UP, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")


MEList <- moduleEigengenes(exp_UP, colors = dynamicColors_UP)
MEs <- MEList$eigengenes
MEDiss <- 1-cor(MEs)
METree <- hclust(as.dist(MEDiss), method = "average")
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
MEDissThres <- 0.25
abline(h = MEDissThres, col = "red")
rm(MEList, MEs, MEDiss, METree)

merge_UP <- mergeCloseModules(exp_UP, dynamicColors_UP, cutHeight = MEDissThres, verbose = 3)
mergedColors_UP <- merge_UP$colors
mergedMEs_UP <- merge_UP$newMEs

plotDendroAndColors(geneTree_UP, cbind(dynamicColors_UP, mergedColors_UP), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)

moduleColors_UP <- mergedColors_UP
colorOrder <- c("grey", standardColors(50))
moduleLabels_UP <- match(moduleColors_UP, colorOrder)-1
MEs_UP <- mergedMEs_UP



# module color-name conversion --------------------------------------------

module_name <- list()
module_name$SP <- data.frame(module_color = setdiff(unique(mergedColors_SP), "grey"), cluster = paste0("SP-", 1:25))
rownames(module_name$SP) <- setdiff(unique(mergedColors_SP), "grey")
module_name$UP <- data.frame(module_color = setdiff(unique(mergedColors_UP), "grey"), cluster = paste0("UP-", 1:38))
rownames(module_name$UP) <- setdiff(unique(mergedColors_UP), "grey")
module_name$SP$size <- table(moduleColors_SP)[module_name$SP$module_color]
module_name$UP$size <- table(moduleColors_UP)[module_name$UP$module_color]


