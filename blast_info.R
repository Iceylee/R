library(dplyr)


genome.list = read.table(file="sorted.genome.list", header=F,stringsAsFactors = F)
#View(genome.list)
colnames(genome.list) <- "Gene_ID"


ardb = read.table(file="ardb.anno.out", sep = "\t", header=F,stringsAsFactors = F)
#View (ardb)
colnames(ardb) <- c("Gene_ID","Swissprot_ID","Identity","Align_Length","Mismatch","Gap","Q_start","Q_end", "T_start", "T_end", "E_value", "Score", "ardb_Function")
ardb_gene_func<- select(ardb,Gene_ID,ardb_Function)
write.csv(ardb,file="ardb.annotation.csv",row.names = F)
# View(ardb_gene_func)

cog = read.table(file="cog.anno.out", sep = "\t", header=F,stringsAsFactors = F,quote="")
colnames(cog) <- c("Gene_ID","Swissprot_ID","Identity","Align_Length","Mismatch","Gap","Q_start","Q_end", "T_start", "T_end", "E_value", "Score", "cog_Function")
cog_gene_func<- select(cog,Gene_ID,cog_Function)
write.csv(cog,file="cog.annotation.csv",row.names = F)

phi = read.table(file="phi.anno.out", sep = "\t", header=F,stringsAsFactors = F,quote="")
colnames(phi) <- c("Gene_ID","Swissprot_ID","Identity","Align_Length","Mismatch","Gap","Q_start","Q_end", "T_start", "T_end", "E_value", "Score", "phi_Function")
phi_gene_func<- select(phi,Gene_ID,phi_Function)
write.csv(phi,file="phi.annotation.csv",row.names = F)

trembl = read.table(file="trembl.anno.out", sep = "\t", header=F,stringsAsFactors = F,quote="")
colnames(trembl) <- c("Gene_ID","Swissprot_ID","Identity","Align_Length","Mismatch","Gap","Q_start","Q_end", "T_start", "T_end", "E_value", "Score", "trembl_Function")
trembl_gene_func<- select(trembl,Gene_ID,trembl_Function)
write.csv(trembl,file="trembl.annotation.csv",row.names = F)

vfdb = read.csv(file="vfdb.anno.out", sep = "\t", header=F,stringsAsFactors = F)
colnames(vfdb) <- c("Gene_ID","Swissprot_ID","Identity","Align_Length","Mismatch","Gap","Q_start","Q_end", "T_start", "T_end", "E_value", "Score", "vfdb_Function")
vfdb_gene_func<- select(vfdb,Gene_ID,vfdb_Function)
write.csv(vfdb,file="vfdb.annotation.csv",row.names = F,quote = F)


nr = read.table(file="nr.anno.out", sep = "\t", header=F,stringsAsFactors = F,quote="")
colnames(nr) <- c("Gene_ID","Swissprot_ID","Identity","Align_Length","Mismatch","Gap","Q_start","Q_end", "T_start", "T_end", "E_value", "Score", "nr_Function")
nr_gene_func<- select(nr,Gene_ID,nr_Function)

sp = read.table(file="swissprot.anno.out", sep = "\t", header=F,stringsAsFactors = F,quote="")
colnames(sp) <- c("Gene_ID","Swissprot_ID","Identity","Align_Length","Mismatch","Gap","Q_start","Q_end", "T_start", "T_end", "E_value", "Score", "swissprot_Function")
sp_gene_func<- select(sp,Gene_ID,swissprot_Function)

#######
result <- left_join(genome.list,ardb_gene_func,by = "Gene_ID")  %>%
   left_join(cog_gene_func,by = "Gene_ID")  %>%
   left_join(phi_gene_func,by = "Gene_ID") %>%
   left_join(trembl_gene_func,by = "Gene_ID")  %>%
   left_join(vfdb_gene_func,by = "Gene_ID") %>%
  left_join(nr_gene_func, by = "Gene_ID") %>%
  left_join(sp_gene_func, by = "Gene_ID")

# count import items
# sum(!is.na(result$phi_Function))

colnames(result) <- c("transcript_ID","ardb","cog","phi","trembl","vfdb","nr","swissprot")

write.csv(result,file="result.csv",na="",row.names = F)

