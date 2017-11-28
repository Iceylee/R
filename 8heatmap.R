library(dplyr)
library(pheatmap)
setwd("~/work/11月/1121cuffdiff/")

fpkm_data <- read.table(file = "Cuffdiff/genes.fpkm_tracking", header=TRUE, sep="\t", stringsAsFactors = F)

fpkm_data2 <- select(fpkm_data,Pal_12dpi_FPKM,Pal.103_12dpi_FPKM,gene_short_name)
colnames(fpkm_data2)[3] <- "gene" 


diff_data <- read.table(file = "gene_exp.diff", header=TRUE, sep="\t", stringsAsFactors = F)
sig_diff <- filter(diff_data,abs(log2.fold_change.) > 1,q_value < 0.05)
write.table(sig_diff, file="sig-genes-cuffdiff.txt", sep="\t",  row.name=F, col.names=TRUE,quote=FALSE)

result <- dplyr::left_join(sig_diff,fpkm_data2,by = "gene")

#get q value top 50
result_q <- arrange(result,q_value) 
result_q_top <- result_q[1:30,15:16] 
rownames(result_q_top) <- result_q[1:30,3]

#log变换
result_q_top_log <- log10(result_q_top + 1)


#添加了legend的标题 log10(FPKM+1)
pdf(file="heatmap_top30.pdf")
pheatmap(result_q_top_log,scale="column", legend_breaks = c(0, 1, 2,max(result_q_top_log)),legend_labels = c("0", "1", "2","log10(FPKM+1)\n"),
         legend = TRUE)
dev.off()


##
# test <- matrix(rexp(200, rate=.1), ncol=20)
# colnames(test) = paste("Room", 1:20, sep = "")
# rownames(test) = paste("Building", 1:10, sep = "")


# pheatmap(test, legend_breaks = c(10, 20, 30, 40, max(test)), 
#          main = "", legend_labels = c("10", "20", "30", "40", "title\n"),
#          legend = TRUE, cluster_rows = FALSE, cluster_cols = FALSE)
