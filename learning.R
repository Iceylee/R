#得到class中每组的count的前10

go_id_top <- go_id %>%
  group_by(class) %>%
  top_n(n = 10, wt = count) %>%
  arrange(class)


#change F,C,P to MF,CC,BP
levels(go_id$class)[levels(go_id$class)=="C"]<- "CC"
levels(go_id$class)[levels(go_id$class)=="F"]<- "MF"
levels(go_id$class)[levels(go_id$class)=="P"]<- "BP"

ids_du <- ids[!duplicated(ids$UNIGENE),]



conda install r
conda install -y bioconductor-deseq bioconductor-deseq2 bioconductor-edger r-gplots


source("http://bioconductor.org/biocLite.R")
biocLite("DESeq")
biocLite("DESeq2")
biocLite("edgeR")