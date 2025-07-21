

# load LINCS L1000 data ---------------------------------------------------

gse70138_info <- list()
gse70138_info$inst_info <- read.delim("LINCS_data/GSE70138/GSE70138_Broad_LINCS_inst_info.txt", sep="\t", stringsAsFactors=F)
gse70138_info$pert_info <- read.delim("LINCS_data/GSE70138/GSE70138_Broad_LINCS_pert_info.txt", sep="\t", stringsAsFactors=F)
gse70138_info$sig_info <- read.delim("LINCS_data/GSE70138/GSE70138_Broad_LINCS_sig_info.txt", sep="\t", stringsAsFactors=F)
gse70138_info$sig_metrics <- read.delim("LINCS_data/GSE70138/GSE70138_Broad_LINCS_sig_metrics.txt", sep="\t", stringsAsFactors=F)
gse70138_info$cell_info <- read.delim("LINCS_data/GSE70138/GSE92742_Broad_LINCS_cell_info.txt", sep="\t", stringsAsFactors=F)
gse70138_info$gene_info <- read.delim("LINCS_data/GSE70138/GSE92742_Broad_LINCS_gene_info.txt", sep="\t", stringsAsFactors=F)

gse92742_info = list()
gse92742_info$inst_info <- read.delim("LINCS_data/GSE92742/GSE92742_Broad_LINCS_inst_info.txt", sep="\t", stringsAsFactors=F)
gse92742_info$pert_info <- read.delim("LINCS_data/GSE92742/GSE92742_Broad_LINCS_pert_info.txt", sep="\t", stringsAsFactors=F)
gse92742_info$pert_info_metrics <- read.delim("LINCS_data/GSE92742/GSE92742_Broad_LINCS_pert_metrics.txt", sep="\t", stringsAsFactors=F)
gse92742_info$sig_info <- read.delim("LINCS_data/GSE92742/GSE92742_Broad_LINCS_sig_info.txt", sep="\t", stringsAsFactors=F)
gse92742_info$sig_metrics <- read.delim("LINCS_data/GSE92742/GSE92742_Broad_LINCS_sig_metrics.txt", sep="\t", stringsAsFactors=F)
gse92742_info$cell_info <- read.delim("LINCS_data/GSE92742/GSE92742_Broad_LINCS_cell_info.txt", sep="\t", stringsAsFactors=F)
gse92742_info$gene_info <- read.delim("LINCS_data/GSE92742/GSE92742_Broad_LINCS_gene_info.txt", sep="\t", stringsAsFactors=F)
gse92742_info$gene_info_landmark <- read.delim("LINCS_data/GSE92742/GSE92742_Broad_LINCS_gene_info_delta_landmark.txt", sep="\t", stringsAsFactors=F)



dys_genes_UP32 <- genes[moduleColors_UP == "greenyellow", ]
dys_genes_UP32_up <- dys_genes_UP32[dys_genes_UP32$log2fc > 0, ]
dys_genes_UP32_down <- dys_genes_UP32[dys_genes_UP32$log2fc < 0, ]
dys_genes_UP32_up <- dys_genes_UP32_up[order(dys_genes_UP32_up$log2fc, decreasing = T), ]
dys_genes_UP32_down <- dys_genes_UP32_down[order(dys_genes_UP32_down$log2fc, decreasing = T), ]
dys_genes_UP32_up <- dys_genes_UP32_up[dys_genes_UP32_up$gene_symbol %in% gse92742_info$gene_info$pr_gene_symbol, ]
dys_genes_UP32_down <- dys_genes_UP32_down[dys_genes_UP32_down$gene_symbol %in% gse92742_info$gene_info$pr_gene_symbol, ]


drugs_JURKAT_KS <- getDrugMatrix("JURKAT", dys_genes_UP32_up, dys_genes_UP32_down)

customized_UP32 <- list()
customized_UP32$up <- dys_genes_UP32_up[, 3:4]
colnames(customized_UP32$up) <- c("up", "gene")
customized_UP32$up$up <- "up_UP32"
customized_UP32$down <- dys_genes_UP32_down[,3:4]
colnames(customized_UP32$down) <- c("down", "gene")
customized_UP32$down$down <- "down_UP32"


drugs_JURKAT_KS$up_ES <- getEnrichedScores(drugs_JURKAT_KS, customized_UP32$up)
drugs_JURKAT_KS$down_ES <- getEnrichedScores(drugs_JURKAT_KS, customized_UP32$down)

drugs_JURKAT_KS$results <- drugs_JURKAT_KS$drugs_info
drugs_JURKAT_KS$results$weighted_connectivity_scores <- rep(0, nrow(drugs_JURKAT_KS$results))

for(i in 1:nrow(drugs_JURKAT_KS$results)){
  if(drugs_JURKAT_KS$up_ES$enrichmentScore[i] * drugs_JURKAT_KS$down_ES$enrichmentScore[i] <= 0)
    drugs_JURKAT_KS$results$weighted_connectivity_scores[i] <- (drugs_JURKAT_KS$up_ES$enrichmentScore[i] - drugs_JURKAT_KS$down_ES$enrichmentScore[i])/2 else
      drugs_JURKAT_KS$results$weighted_connectivity_scores[i] <- 0}

drugs_JURKAT_KS$results <- drugs_JURKAT_KS$results[order(drugs_JURKAT_KS$results$weighted_connectivity_scores, decreasing = T), ]
mean_up <- mean(drugs_JURKAT_KS$results$weighted_connectivity_scores[drugs_JURKAT_KS$results$weighted_connectivity_scores > 0])
mean_down <- mean(drugs_JURKAT_KS$results$weighted_connectivity_scores[drugs_JURKAT_KS$results$weighted_connectivity_scores < 0])


drugs_JURKAT_KS$results$KS <- drugs_JURKAT_KS$results$weighted_connectivity_scores
drugs_JURKAT_KS$results$KS[drugs_JURKAT_KS$results$KS > 0] <- 
  drugs_JURKAT_KS$results$KS[drugs_JURKAT_KS$results$KS > 0]/max(drugs_JURKAT_KS$results$KS)
drugs_JURKAT_KS$results$KS[drugs_JURKAT_KS$results$KS < 0] <- 
  drugs_JURKAT_KS$results$KS[drugs_JURKAT_KS$results$KS < 0]/-min(drugs_JURKAT_KS$results$KS)


