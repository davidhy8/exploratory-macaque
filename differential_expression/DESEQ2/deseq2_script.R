library(DESeq2)
library(dplyr)

rawCounts <- read.csv("expression_table_raw.csv", header = TRUE, sep=',', row.names = NULL)
rawCounts <- rawCounts %>% mutate_if(is.numeric, round)
head(rawCounts)

sampleData <- read.table("sample_file.txt", header = TRUE, sep='\t', row.names = NULL)
head(sampleData)

dds <- DESeqDataSetFromMatrix(countData=rawCounts,colData=sampleData,design=~rnfl_class + sex + location + tissue + age, tidy = TRUE)
dds_lrt <- DESeq(dds, test='LRT', reduced =~rnfl_class + sex + location + tissue)
dds_wald <- DESeq(dds)

res_lrt <- results(dds_lrt)
res_wald <- results(dds_wald)

#head(results(dds, tidy=TRUE))
lrt = mcols(res_lrt, use.names = TRUE)
wald = mcols(res_wald, use.names = TRUE)

write.csv(res_lrt, file='all_age_lrt.csv')
write.csv(res_wald, file='all_age_wald.csv')

pdf('volcano_plot_age_lrt.pdf')
print(par(mfrow=c(1,1)))
# Make a basic volcano plot
print(with(res_lrt, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3))))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
print(with(subset(res_lrt, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue")))
print(with(subset(res_lrt, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red")))
dev.off()


pdf('volcano_plot_age_wald.pdf')
print(par(mfrow=c(1,1)))
# Make a basic volcano plot
print(with(res_wald, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3))))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
print(with(subset(res_wald, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue")))
print(with(subset(res_wald, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red")))
dev.off()
