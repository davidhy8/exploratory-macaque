getOption("clusterProfiler.download.method")
install.packages("R.utils")


library(clusterProfiler)
library(msigdbr)
library(data.table)
library("org.Mmu.eg.db")
# library(tidyr)

R.utils::setOption("clusterProfiler.download.method","auto")

full_data <- fread("kallisto_wald_test_lrt_gene_passed_all_age.csv") 

# ID_conversion <- fread("genesymbol_entrezid_key.txt")

curr_gene_list = full_data$target_id

# gene_list = match(curr_gene_list, ID_conversion$gene_symbol)

# gene_list <- gene_list[!is.na(gene_list)]

# gene_list <- unlist(strsplit(as.character(gene_list), " "))

# df <- data.frame(gene_list)

# gene_list = df$gene_list

# gene_list = unique(gene_list)

gene_list = bitr(curr_gene_list, fromType='SYMBOL', toType='ENTREZID', OrgDb = 'org.Mmu.eg.db')

gene_list = gene_list$ENTREZID

# all_gene_sets = msigdbr(species = "Macaca mulatta")

# msigdbr_df <- msigdbr(species = "Macaca mulatta")

# msigdbr_t2g = msigdbr_df %>% dplyr::distinct(gs_name, ensembl_gene) %>% as.data.frame()


# KEGG analysis
# enriched = enricher(gene = curr_gene_list, TERM2GENE = msigdbr_t2g, qvalueCutoff = 0.05)

KEGG = enrichKEGG(gene_list, organism='mcc', pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.05)
summary(KEGG)

file_name = 'KEGG_dot_all_age_sleuth.pdf'
pdf(file_name)
print(dotplot(KEGG, showCategory=30))
dev.off()

# Can choose to view the various KEGG pathways enriched
# browseKEGG(KEGG, 'mcc04015')
# browseKEGG(KEGG, 'mcc04724')
# browseKEGG(KEGG, 'mcc04713')
# pdf(file="keggDot_CF.pdf")
# dotplot(enriched, showCategory=30)
# dev.off()


# GO analysis
GO = enrichGO(gene = gene_list, OrgDb = org.Mmu.eg.db, ont = 'CC', pAdjustMethod='BH', pvalueCutoff = 0.05, qvalueCutoff=0.05)
head(summary(GO))

file_name = 'GO_dot_all_age_sleuth.pdf'
pdf(file_name)
print(dotplot(GO, showCategory=30))
dev.off()

file_name = 'GO_graph_all_age_sleuth.pdf'
pdf(file_name)
print(goplot(GO))
dev.off()

go_df = as.data.frame(GO)
go_df$"Method" = "GO"
kegg_df = as.data.frame(KEGG)
kegg_df$'Method' = "KEGG"
final_results = rbind(go_df, kegg_df)
fwrite(final_results, 'all_age_PE.csv')

