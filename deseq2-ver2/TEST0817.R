Rscript /data1/script/deseq2+GO+KEGG/Rpipe/3GO-KEGG.R AH57973 hsa SYMBOL ENTREZID /data2/ClientData/2018_08/Lishuang/


dbname <- 'AH57973'
kegg_org <- "hsa"
GO_KEY <- "SYMBOL" 
KEGG_NEED_KEY <- "ENTREZID" 
output_path <- "/data2/ClientData/2018_08/XuLeLe/"

#跑前删除pathway文件夹。不然有可能报错
#模式生物：sig files 中有log2FoldChange列
#非模式生物：第四列为log2FoldChange的数值。手动加列名

Rscript /data1/script/deseq2+GO+KEGG/Rpipe/pathview.R ENSEMBL 0.05 FALSE ./nonModel/4.GO_KEGG_Enrichment/ ./nonModel/3.DiffExprGene/

Rscript /data1/script/deseq2+GO+KEGG/Rpipe/pathview.R ENSEMBL 0.05 TRUE ./Model/4.GO_KEGG_Enrichment/ ./Model/3.DiffExprGene/


Rscript pathview.R ENSEMBL 0.05 FALSE ./nonModel/4.GO_KEGG_Enrichment/ ./nonModel/3.DiffExprGene/

Rscript pathview.R ENSEMBL 0.05 TRUE ./result/4.GO_KEGG_Enrichment/ ./result/3.DiffExprGene/