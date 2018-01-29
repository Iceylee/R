
library(dplyr)

df <- read.csv("kegg.count", sep = '\t',header = F)
colnames(df) <- c("Class","Description","Count")

df <- df %>%
  top_n(n = 20, wt = Count) %>%
  arrange(Class,Count)

df$Description<- factor(df$Description, order=TRUE, levels=df$Description)
df$onco<- factor(df$onco, order=T)
###
library(ggplot2)
p <- ggplot(df, aes(y = Count, x = Description)) +
  geom_bar(stat = "identity", aes(fill = Class), alpha = 1) +
  scale_fill_discrete(name = "Ontology") +
  theme_light() +
  theme(axis.text = element_text(size = 7), legend.text = element_text(size = 8)) +
  labs(y = "Number of Genes", x = "Term")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 
#palette = "Spectral"
p <- p +scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1, vjust = 1))

ggsave("KEGG_barplot.pdf",p, width=10, height=7, units="in")




