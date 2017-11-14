library(clusterProfiler)
library(DOSE)

#ensure organism available
search_kegg_organism('vda', by='kegg_code')

gene_lists <- read.csv("sig-genes-list.txt", sep = '\t',header = F)
gene_id <- as.character(gene_lists[,1])

#ids <- bitr(x, fromType="SYMBOL", toType=c("UNIPROT", "ENSEMBL"), OrgDb="org.Hs.eg.db")
#head(ids)

kk <- enrichKEGG(gene = gene_id, organism = "vda", keyType = "kegg", pvalueCutoff = 0.05)

kk_df <- as.data.frame(kk) %>%
  dplyr::select(-ID)

write.table(kk_df, file="KEGG_out.txt")

##plot
pdf(file="KEGG_dotplot.pdf")
dotplot(kk)
dev.off()
