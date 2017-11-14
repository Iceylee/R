#Zebrafish
library(org.Dr.eg.db)
#mouse
# source("https://bioconductor.org/biocLite.R")
# biocLite("org.Mm.eg.db")
# library(org.Mm.eg.db)  
#human
#library(org.Hs.eg.db)


library(clusterProfiler)
library(DOSE)

gene_lists <- read.csv("sig-genes-deseq2.txt", sep = '\t',header = T)
gene_id <- as.character(gene_lists[,1])

info_go_BP <- enrichGO(gene = gene_id, 
                       OrgDb = org.Dr.eg.db, 
                       keytype = "ENSEMBL", 
                       ont = "BP", 
                       pAdjustMethod = "BH", 
                       pvalueCutoff = 0.01, 
                       qvalueCutoff = 0.05)
info_go_CC <- enrichGO(gene = gene_list, 
                       OrgDb = org.Dr.eg.db, 
                       keytype = "ENSEMBL", 
                       ont = "CC", 
                       pAdjustMethod = "BH", 
                       pvalueCutoff = 0.01, 
                       qvalueCutoff = 0.05)
info_go_MF <- enrichGO(gene = gene_list, 
                       OrgDb = org.Dr.eg.db, 
                       keytype = "ENSEMBL", 
                       ont = "MF", 
                       pAdjustMethod = "BH", 
                       pvalueCutoff = 0.01, 
                       qvalueCutoff = 0.05)

#output table
write.table(as.data.frame(info_go_BP@result), file="GO_BP_out.txt")
write.table(as.data.frame(info_go_CC@result), file="GO_CC_out.txt")
write.table(as.data.frame(info_go_MF@result), file="GO_MF_out.txt")

#get gene sum number
tep_num1 = as.data.frame(info_go_CC)[1,3]
tep_num2 = strsplit(tep_num1 , "/")
gene_num = as.numeric(tep_num2[[1]][2])

#p value top 15 select
ego_CC_df <- as.data.frame(info_go_CC)[1:15, c(2:3, 9)] %>%
  transform(onco = "Cellular component", GeneRatio = Count / gene_num * 100)

ego_BP_df <- as.data.frame(info_go_BP)[1:15, c(2:3, 9)] %>%
  transform(onco = "Biological process", GeneRatio = Count / gene_num * 100)

ego_MF_df <- as.data.frame(info_go_MF)[1:8, c(2:3, 9)] %>%
  transform(onco = "Molecular function", GeneRatio = Count / gene_num * 100)


##bind three onco
ego_three <- rbind(ego_CC_df, ego_BP_df, ego_MF_df)

##plot bar
#########NEED CHANGE limits
library(ggplot2)
p <- ggplot(ego_three, aes(y = GeneRatio, x = Description)) +
  geom_bar(stat = "identity", aes(fill = onco), alpha = 0.7) +
  facet_grid(onco ~ ., scales = "free") +
  coord_flip()  +
  scale_y_continuous(limits = c(0, 5))+
  scale_fill_discrete(name = "Ontology", labels = c("Biological process", "Molecular function", "Cellular component")) +
  theme_light() +
  theme(axis.text = element_text(size = 6, face = "bold"), legend.text = element_text(size = 6, face = "bold")) +
  labs(y = "Percentage of Genes", x = "Term")

pdf(file="GO_barplot.pdf")
p
dev.off()

#KEGG
########NEED CHANGE organism#########
#info_ids <- bitr(info[,1], fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Dr.eg.db)
kk <- enrichKEGG(gene = gene_lists, organism = "dre", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05)

kk_df <- as.data.frame(kk) %>%
  dplyr::select(-ID)

write.table(kk_df, file="KEGG_out.txt")

##plot
pdf(file="KEGG_dotplot.pdf")
dotplot(kk)
dev.off()