lrt = read.csv("all_age_lrt.csv")
wald = read.csv("all_age_wald.csv")

lrt = lrt[lrt$padj < 0.05, ]
wald = wald[wald$padj < 0.05, ]

wald = wald[rowSums(is.na(wald)) == 0, ]
lrt = lrt[rowSums(is.na(lrt)) == 0, ]

# final = rbind(wald,lrt)
lrt.sig.ids <- lrt$gene
wd.sig.ids <- wald$gene
shared_ids <- intersect(lrt.sig.ids, wd.sig.ids)
final = lrt[lrt$gene %in% shared_ids,]

write.csv(final, "all_age_DE.csv", row.names = FALSE)
