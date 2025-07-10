pacman::p_load(tidyverse)


enrich_go_cp <- function(x, module_color){
  
  require(clusterProfiler)
  require(org.Hs.eg.db)

  entrezid <- genes$gene_entrez[module_color == x]
  
  ego <- enrichGO(gene = entrezid,
                  OrgDb = org.Hs.eg.db,
                  keyType = "ENTREZID",
                  ont = "ALL",
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.05,
                  readable = TRUE)
  
  ego <- ego@result %>% 
    mutate(module = all_of(x)) %>% 
    relocate(module, .before = ONTOLOGY)
  
  return(ego)}


plot_GO_clusters <- function(x, name){
  
  require(ggplot2)
  require(ggthemes)
  
  x$log10_p_value <- -log10(x$p.adjust)
  data_table <- names(table(x$module))
  data_table <- cbind(data_table, rep(0, length(data_table)))
  for(i in 1:nrow(data_table))
    data_table[i, 2] <- x[x$module == data_table[i, 1], 12][1]
  
  colnames(data_table) <- c("cluster", "adj.p.value")
  data_table <- as.data.frame(data_table)
  data_table$adj.p.value <- as.numeric(data_table$adj.p.value)
  data_table$cluster <- factor(data_table$cluster, levels = data_table$cluster[order(data_table$adj.p.value)])
  data_table <- data_table[order(data_table$adj.p.value, decreasing = T), ][1:10, ]
  
  p <- ggplot(data_table, aes(x = cluster, y = adj.p.value)) +
    geom_bar(stat = "identity", fill = "grey40") + theme_minimal() + #add colors
    ylab(expression('-Log'[10]*'(adjusted p-value)')) +
    theme(plot.title = element_text(hjust = 0.5, size = 30, face = "bold")) +
    scale_y_continuous(limits = c(0, 33), expand = c(0, 0)) +
    theme(axis.text.x = element_text(size = 15, colour = "black"), 
          axis.line = element_line(colour = "black", linewidth = 0.5),
          axis.title.x = element_text(size = 15, colour = "black"), 
          axis.title.y = element_blank(), 
          axis.text.y = element_text(size = 15, colour = "black"),
          plot.title = element_text(size = 20),
          axis.ticks.y = element_blank(),
          legend.position = "none",
          aspect.ratio = 0.9,
          plot.background=element_blank()) +
    ggtitle(name) +
    coord_flip()
  
  return(p)}


hyper_test <- function(set1, set2, all = 13641)
  return(phyper(sum(set1 %in% set2)-1, length(set1), all-length(set1), length(set2), lower.tail = F))



loadGCTXData = function(GSE, cell_line){
  
  library(cmapR)
  if(GSE == "GSE70138"){
    GSE_info <- gse70138_info
    ds_path <- "LINCS_data/GSE70138/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328.gctx"}
  
  if(GSE == "GSE92742"){
    GSE_info <- gse92742_info
    ds_path <- "LINCS_data/GSE92742/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"}
  gene_info <- GSE_info$gene_info
  
  idx <- which(GSE_info$sig_info$cell_id == cell_line)
  sig_ids <- GSE_info$sig_info$sig_id[idx]
  drug_ds <- parse.gctx(ds_path, cid = sig_ids, rid = as.character(GSE_info$gene_info$pr_gene_id))
  
  drug_matrix <- list()
  drug_matrix$gene_info_drug <- gene_info
  drug_matrix$drug_matrix <- drug_ds
  return(drug_matrix)}


getEnrichedScores <- function(drugs, target_list_up){
  
  results <- drugs$drugs_info
  results$enrichmentScore <- rep(0, nrow(results))
  results$NES <- rep(0, nrow(results))
  results$pvalue <- rep(0, nrow(results))
  results$p.adjust <- rep(0, nrow(results))
  
  for(i in results$sig_id){
    gene_list <- drugs$drug_FC_matrix[,i]
    names(gene_list) <- drugs$gene_info$pr_gene_symbol
    gene_list <- sort(gene_list, decreasing = T)
    gsea_results <- GSEA(gene_list, TERM2GENE = target_list_up, pvalueCutoff = 1, verbose=FALSE)
    results[results$sig_id == i,9:12] <- gsea_results@result[4:7]}
  
  return(results)}


getDrugMatrix <- function(cell_line, target_list_up, target_list_down){
  
  results <- list()
  if(cell_line %in% gse92742_info$sig_info$cell_id)
    gse92742_drug_data <- loadGCTXData("GSE92742", cell_line)
  
  if(cell_line %in% gse70138_info$sig_info$cell_id)
    gse70138_drug_data <- loadGCTXData("GSE70138", cell_line)
  
  if(exists("gse92742_drug_data") & exists("gse70138_drug_data")){
    results$drug_FC_matrix <- cbind(gse92742_drug_data$drug_matrix@mat, gse70138_drug_data$drug_matrix@mat)
    results$drugs_info <- rbind(gse92742_info$sig_info[,c(1:5, 8, 11, 12)], gse70138_info$sig_info)
    results$gene_info <- gse92742_drug_data$gene_info_drug}
  
  if(!exists("gse92742_drug_data") & exists("gse70138_drug_data")){
    print("only gse70138")
    results$drug_FC_matrix <- gse70138_drug_data$drug_matrix@mat
    results$drugs_info <- gse70138_info$sig_info
    results$gene_info <- gse70138_drug_data$gene_info_drug}
  
  if(exists("gse92742_drug_data") & !exists("gse70138_drug_data")){
    print("only gse92742")
    results$drug_FC_matrix <- gse92742_drug_data$drug_matrix@mat
    results$drugs_info <- gse92742_info$sig_info[,c(1:5, 8, 11, 12)]
    results$gene_info <- gse92742_drug_data$gene_info_drug}
  
  results$drugs_info <- results$drugs_info[results$drugs_info$cell_id == cell_line,]
  results$drugs_info$pert_idose <- gsub(results$drugs_info$pert_idose, pattern = "ÂµM", replacement = "µM")
  results$target_list_up <- target_list_up
  results$target_list_down <- target_list_down
  return(results)}


calc_ab <- function(target_list, drug_index, gene_rank){
  
  num_genes <- nrow(target_list)
  target_list_rank <- gene_rank[gene_rank %in% target_list$SYMBOL]
  aa <- rep(0, num_genes)
  bb <- rep(0, num_genes)
  
  for(j in 1:num_genes){
    aa[j] <- j/num_genes - which(gene_rank == target_list_rank[j])/length(gene_rank)
    bb[j] <- which(gene_rank == target_list_rank[j])/length(gene_rank) - (j-1)/num_genes}
  aa <- max(aa)
  bb <- max(bb)
  
  if(aa > bb)
    result <- aa
  if(bb > aa)
    result <- -bb
  
  return(result)}



set.seed(123)
sc_GSE159677_col <- scales::hue_pal()(11) %>% sample
names(sc_GSE159677_col) <-
  c("B cell", "Plasma cell", "T cell", "NK cell", "Mono/Macro/DC", "Neutrophil", "EC", "Fibroblast/SMC", "Fibroblast", "Mast cell", "Proliferating cell")

set.seed(123)
sc_GSE159677_tsub_col <- scales::hue_pal()(7) %>% sample
names(sc_GSE159677_tsub_col) <-
  c("CD4 Tnaive", "CD4 Tmem", "CD4 Th17", "CD4 Treg", "CD8 Tem", "CD8 Trm", "CD8 CTL")


set.seed(123)
sc_GSE224273_col <- scales::hue_pal()(10) %>% sample
names(sc_GSE224273_col) <-
  c("B cell", "Plasma cell", "CD4 T cell", "CD8 T cell", "NK cell", "Myeloid", "EC", "Fibroblast/SMC", "Mast cell", "Proliferating cell")

set.seed(123)
sc_GSE224273_tsub_col <- scales::hue_pal()(7) %>% sample
names(sc_GSE224273_tsub_col) <-
  c("CD4 Tnaive", "CD4 Tmem", "CD4 Treg", "CD8 Tnaive", "CD8 Tem", "CD8 MAIT", "CD8 CTL")

