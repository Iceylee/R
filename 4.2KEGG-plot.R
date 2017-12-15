setwd("~/Downloads/work/10月/11.1 /")
library(dplyr)


kegg_list <- read.csv("sig-genes-deseq2_ID_identify.out", sep = '\t',header = F,stringsAsFactors = F)
colnames(kegg_list) <- c("Term","Database","ID", "Input_number",    "Background_number", "PValue"," CorrectedPValue", "Input","Hyperlink")



kegg_list$Input_number <- as.numeric(kegg_list$Input_number)
kegg_list$PValue <- as.numeric(kegg_list$PValue)

y = kegg_list %>% arrange(PValue) 

kegg_top <- y[1:10,]



total <- sum(kegg_list$Input_number)
kegg_top <- kegg_top %>% mutate(GeneRatio = Input_number/total) %>%
  arrange(GeneRatio)
kegg_top$Term <- as.factor(kegg_top$Term)


library(ggplot2)
p <- ggplot(kegg_top, aes(x = GeneRatio,y=reorder(Term,GeneRatio))) +
  geom_point(aes(size=Input_number,color = PValue)) +
  scale_size("Count") +
  scale_color_continuous(low='red', high='grey') +
  theme_light() +
  labs(x="Gene Ratio", y = "Term")

pdf(file="KEGG_dotplot.pdf")
p
dev.off()

###
##########用clusterprofiler的kegg画图

kegg_list <- read.table(file = "KEGG_out.txt", sep = '\t',header =T,stringsAsFactors = F)

#get gene sum number
#可能的问题：BP结果为空，无法计算
tep_num1 = as.data.frame(kegg_list)[1,2]
tep_num2 = strsplit(tep_num1 , "/")
gene_num = as.numeric(tep_num2[[1]][2])

kegg_list <- kegg_list %>% mutate(Gene_Ratio = Count/gene_num) %>%
  arrange(GeneRatio)



library(ggplot2)
p <- ggplot(kegg_list, aes(x = Gene_Ratio,y=reorder(Description,Gene_Ratio))) +
  geom_point(aes(size=Count,color = pvalue)) +
  scale_size("Count") +
  scale_color_continuous(low='red', high='grey') +
  theme_light() +
  labs(x="Gene Ratio", y = "Term",color="PValue")

#pdf(file="KEGG_dotplot2.pdf")
p = p + guides(
  color = guide_colorbar(order = 1),
  fill = guide_legend(order = 0)
)
# dev.off()
ggsave("KEGG_dotplot2.pdf",p, width=8, height=8, units="in")

