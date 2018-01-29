
library(dplyr)


df <- read.csv("go_count.out", sep = '\t',header = F)
colnames(df) <- c("GO","Count","onco","Description")

df_CC <- filter(df,onco == "cellular_component") %>%
  arrange(desc(Count))
df_MF <- filter(df,onco == "molecular_function") %>%
  arrange(desc(Count))
df_BP <- filter(df,onco == "biological_process") %>%
  arrange(desc(Count))

CC_sum = sum(df_CC$Count)
BP_sum = sum(df_BP$Count)
MF_sum = sum(df_MF$Count)


#p value top 15 select
ego_CC_df <- df_CC[1:10,] %>%
  transform(GeneRatio = Count / CC_sum * 100)
ego_MF_df <- df_MF[1:10,] %>%
  transform(GeneRatio = Count / MF_sum * 100)
ego_BP_df <- df_BP[1:10,] %>%
  transform(GeneRatio = Count / BP_sum * 100)


##bind three onco
ego_three <- rbind(ego_BP_df,ego_CC_df,ego_MF_df)

ego_three<- ego_three %>%
  group_by(onco) %>%
  top_n(n = 15, wt = Count) %>%
  arrange(onco,Count)

ego_three$Description<- factor(ego_three$Description, order=TRUE, levels=ego_three$Description)
#ego_three$onco<- factor(ego_three$onco, order=TRUE, levels=c("Molecular function","Cellular component","Biological process"))
ego_three$onco<- factor(ego_three$onco, order=TRUE)
levels(ego_three$onco) <- c("BP","CC","MF")

lable_name <- ego_three$onco[!duplicated(ego_three$onco)]
###
p <- ggplot(ego_three, aes(y = Count, x = Description)) +
  geom_bar(stat = "identity", aes(fill = onco), alpha = 1) +
  facet_grid(onco ~ ., scales = "free", space = "free",margins = F) +
  coord_flip()  +
  #scale_y_continuous(limits = c(0, 70))+
  scale_fill_discrete(name = "Ontology", labels = lable_name) +
  theme_light() +
  theme(axis.text = element_text(size = 9), legend.text = element_text(size = 8)) +
  labs(y = "Number of Genes", x = "Term")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 

pdf(file="GO_barplot.pdf")
p
dev.off()



