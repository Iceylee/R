library(clusterProfiler)
library(DOSE)

library(dplyr)

require(AnnotationHub)
hub <- AnnotationHub()
query(hub, "pseudomonas putida")
maize <- hub[['AH56235']]

length(keys(maize))
columns(maize)
keys(maize)[1:100]


gene_list <- read.csv("swissprot-trezID.txt", sep = '\t',header = F)
gene_id <- gene_list$V1

##enrichGO
info_go_BP <- enrichGO(gene = gene_id, 
                       OrgDb = maize, 
                       keytype = "ENTREZID", 
                       ont = "BP", 
                       pAdjustMethod = "BH", 
                       pvalueCutoff = 0.05, 
                       qvalueCutoff = 0.05)
info_go_CC <- enrichGO(gene = gene_id, 
                       OrgDb = maize, 
                       keytype = "ENTREZID", 
                       ont = "CC", 
                       pAdjustMethod = "BH", 
                       pvalueCutoff = 0.05, 
                       qvalueCutoff = 0.05)
info_go_MF <- enrichGO(gene = gene_id, 
                       OrgDb = maize, 
                       keytype = "ENTREZID", 
                       ont = "MF", 
                       pAdjustMethod = "BH", 
                       pvalueCutoff = 0.05, 
                       qvalueCutoff = 0.05)

#output table
write.table(as.data.frame(info_go_BP@result), file="GO_BP_out.txt")
write.table(as.data.frame(info_go_CC@result), file="GO_CC_out.txt")
write.table(as.data.frame(info_go_MF@result), file="GO_MF_out.txt")

#get gene sum number
tep_num1 = as.data.frame(ego_CC)[1,3]
tep_num2 = strsplit(tep_num1 , "/")
gene_num_CC = as.numeric(tep_num2[[1]][2])

tep_num1 = as.data.frame(ego_MF)[1,3]
tep_num2 = strsplit(tep_num1 , "/")
gene_num_MF = as.numeric(tep_num2[[1]][2])

tep_num1 = as.data.frame(ego_BP)[1,3]
tep_num2 = strsplit(tep_num1 , "/")
gene_num_BP = as.numeric(tep_num2[[1]][2])


ego_CC <- as.data.frame(info_go_CC)
ego_BP <- as.data.frame(info_go_BP)
ego_MF <- as.data.frame(info_go_MF)


#p value top 15 select
ego_CC_df <- as.data.frame(info_go_CC)[1:3, c(2:3, 9)] %>%
  transform(onco = "CC", GeneRatio = Count / gene_num_CC * 100)

ego_BP_df <- as.data.frame(info_go_BP)[1:10, c(2:3, 9)] %>%
  transform(onco = "BP", GeneRatio = Count / gene_num_BP * 100)

ego_MF_df <- as.data.frame(info_go_MF)[1:10, c(2:3, 9)] %>%
  transform(onco = "MF", GeneRatio = Count / gene_num_MF * 100)


##bind three onco
ego_three <- rbind(ego_CC_df, ego_BP_df, ego_MF_df)
View(ego_three)
##plot bar
#########NEED CHANGE limits
library(ggplot2)
p <- ggplot(ego_three, aes(y = GeneRatio, x = Description)) +
  geom_bar(stat = "identity", aes(fill = onco), alpha = 0.7) +
  facet_grid(onco ~ ., scales = "free", space = "free") +
  coord_flip()  +
  scale_y_continuous(limits = c(0, 70))+
  scale_fill_discrete(name = "Ontology", labels = c("Biological process", "Molecular function", "Cellular component")) +
  theme_light() +
  theme(axis.text = element_text(size = 6, face = "bold"), legend.text = element_text(size = 6, face = "bold")) +
  labs(y = "Percentage of Genes", x = "Term")

pdf(file="GO_barplot.pdf")
p
dev.off()