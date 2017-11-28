setwd("~/work/11月/1113blast/")
library(dplyr)

dir.create("out/")
genome.list = read.table(file="sorted.genome.list", header=F,stringsAsFactors = F)
result <- genome.list
#View(genome.list)
colnames(result) <- "Gene_ID"

anno_list = c("cog","go","kegg","nr","pfam","swissprot","trembl")


for (i in anno_list){
  filename = paste(i,".anno.out", sep="")
  outfilename = paste("out/",i,".annotation.csv",sep="")
  df = read.table(file = filename, sep = "\t", header=F,stringsAsFactors = F,quote="")
  
  if (i == "go" || i == "kegg"){
    colnames(df) <- c("Gene_ID", "Function")
    gene_func <- df
  }
  else{
    colnames(df) <- c("Gene_ID",paste(i,"_ID",sep=""),"Identity","Align_Length","Mismatch","Gap","Q_start","Q_end", "T_start", "T_end", "E_value", "Score", "Function")
    gene_func<- select(df,Gene_ID,Function)
  }
  write.csv(df,file = outfilename,row.names = F,quote = F)
  result <- left_join(result,gene_func,by = "Gene_ID") 
}

colnames(result) <- c("transcript_ID","cog","go","kegg","nr","pfam","swissprot","trembl")
write.table(result,file="out/result.txt",na="",row.names = F,quote = F,sep = "\t")
