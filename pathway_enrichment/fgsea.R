library(fgsea)
library(msigdbr)
library(ggplot2)
library(data.table)

#data(examplePathways)
#data(exampleRanks)

msigdbr_df <- msigdbr(species = "Macaca mulatta", category = "H")
msigdbr_df_all <- msigdbr(species = "Macaca mulatta")

pathwaysH = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)

pathways = split(x = msigdbr_df_all$gene_symbol, f = msigdbr_df_all$gs_name)


### FOR LRT RESULTS
data = read.csv('all_age_lrt.csv')
#remove rows with NA
data = data[complete.cases(data), ]

rankData <- data$log2FoldChange
names(rankData) <- data$gene

head(rankData)

fgseaRes <- fgsea(pathways=pathwaysH, stats=rankData, minSize=15, maxSize = 500, nperm = 100)
fgseaRes_all <- fgsea(pathways=pathways, stats=rankData, minSize=15, maxSize = 500, nperm = 100)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()

pdf("fgsea_ggplot_all_age_lrt_hallmark.pdf")
print(ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
        geom_col(aes(fill=padj<0.05)) +
        coord_flip() +
        labs(x="Pathway", y="Normalized Enrichment Score",
             title="Pathways NES from GSEA") + 
        theme_minimal())
dev.off()

topPathwaysUp <- fgseaRes_all[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes_all[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

pdf("gsea_table_all_age_lrt.pdf")
print(plotGseaTable(pathways[topPathways], rankData, fgseaRes_all, 
              gseaParam=1))
dev.off()

fwrite(fgseaRes_all, file="fgseaRes_all_age_lrt.txt", sep="\t", sep2=c("", " ", ""))


### FOR WALD TEST RESULTS
data = read.csv('all_age_wald.csv')
#remove rows with NA
data = data[complete.cases(data), ]

rankData <- data$log2FoldChange
names(rankData) <- data$gene

head(rankData)

fgseaRes <- fgsea(pathways=pathwaysH, stats=rankData, minSize=15, maxSize = 500, nperm = 100)
fgseaRes_all <- fgsea(pathways=pathways, stats=rankData, minSize=15, maxSize = 500, nperm = 100)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()

pdf("fgsea_ggplot_all_age_wald_hallmark.pdf")
print(ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
        geom_col(aes(fill=padj<0.05)) +
        coord_flip() +
        labs(x="Pathway", y="Normalized Enrichment Score",
             title="Pathways NES from GSEA") + 
        theme_minimal())
dev.off()

topPathwaysUp <- fgseaRes_all[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes_all[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

pdf("gsea_table_all_age_wald.pdf")
print(plotGseaTable(pathways[topPathways], rankData, fgseaRes_all, 
                    gseaParam=1))
dev.off()

fwrite(fgseaRes_all, file="fgseaRes_all_age_wald.txt", sep="\t", sep2=c("", " ", ""))





