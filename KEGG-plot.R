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


