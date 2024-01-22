# remotes::install_github("YuLab-SMU/clusterProfiler")
getOption("clusterProfiler.download.method")
R.utils::setOption("clusterProfiler.download.method","auto")

library(clusterProfiler)
library(msigdbr)
library(data.table)
library("org.Mmu.eg.db")

set.seed(10)

full_data <- fread("stabilized_lasso_SecondaryGenes_all_age_converted.csv")
# ID_conversion <- fread("DE_sleuth/mm10_ref2gene_gene.txt")
# ID_conversion_2 <- fread("DE_sleuth/mm10_ref2gene_ensembl.txt")
# all_gene_sets = msigdbr(species = "Macaca mulatta")
# msigdbr_df <- msigdbr(species = "Macaca mulatta")
# msigdbr_t2g = msigdbr_df %>% dplyr::distinct(gs_name, ensembl_gene) %>% as.data.frame()

core_genes = unique(full_data$'Core Gene')

final_results = NULL

for(currCoreGene in core_genes){
  
  curr_gene_list <- subset(full_data, full_data$'Core Gene' == currCoreGene)
  curr_gene_list = curr_gene_list$`Secondary Gene`
  
  # gene_list = bitr(curr_gene_list, fromType='ENSEMBL', toType='ENTREZID', OrgDb = 'org.Mmu.eg.db')
  gene_list = bitr(curr_gene_list, fromType='SYMBOL', toType='ENTREZID', OrgDb = 'org.Mmu.eg.db')
  gene_list = gene_list$ENTREZID
  
  enrichment_KEGG = NULL
  
  enrichment_KEGG = as.data.frame(enrichKEGG(gene_list, organism='mcc', pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.05))
  KEGG = enrichKEGG(gene_list, organism='mcc', pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.05)
  
  if(is.null(enrichment_KEGG) | nrow(enrichment_KEGG) == 0){
    print("No KEGG pathways enriched. Skipping...")
  }
  else{
    enrichment_KEGG$'Core Gene' <- currCoreGene
    enrichment_KEGG$'Method' <- 'KEGG'
    
    if(exists("final_results")){
      final_results <- rbindlist(list(final_results, enrichment_KEGG))
    }
    
    else{
      final_results <- enrichment_KEGG
    }
    
    file_name = paste(currCoreGene, 'KEGG_all_age_converted.pdf', sep='_')
    pdf(file_name)
    print(dotplot(KEGG, showCategory=30))
    dev.off()
  }
  
  
  
  enrichment_GO = as.data.frame(enrichGO(gene = gene_list, OrgDb = org.Mmu.eg.db, ont = 'CC', pAdjustMethod='BH', pvalueCutoff = 0.01, qvalueCutoff=0.05))
  if(is.null(enrichment_GO) | nrow(enrichment_GO) == 0){
    print("No GO pathways enriched. Skipping...")
    next
  } 
  GO = enrichGO(gene = gene_list, OrgDb = org.Mmu.eg.db, ont = 'CC', pAdjustMethod='BH', pvalueCutoff = 0.01, qvalueCutoff=0.05)
  enrichment_GO$'Core Gene' <- currCoreGene
  enrichment_GO$'Method' <- 'GO'
  
  if(exists("final_results")){
    final_results <- rbindlist(list(final_results, enrichment_GO))
  }else{
    final_results <- enrichment_GO
  }
  
  file_name = paste(currCoreGene, 'GO_dot_all_age_converted.pdf', sep='_')
  pdf(file_name)
  print(dotplot(GO, showCategory=30))
  dev.off()
  
  file_name = paste(currCoreGene, 'GO_graph_all_age_converted.pdf', sep='_')
  pdf(file_name)
  print(goplot(GO))
  dev.off()
}

fwrite(final_results, 'scope_all_age_converted.csv')





