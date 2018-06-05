#!/usr/bin/env Rscript

args<-commandArgs(T)

# test if there is at least one argument: if not, return an error
if (length(args)!=4) {
	print ("4 arguments must be supplied:")
	print ("matrixFile groupFile baseGroup outputPath")
	print ("CountMatrix4DESeq.csv colData.csv NC ./")
  stop(call.=FALSE)
} 


count_table <- args[1] #"CountMatrix4DESeq.csv"
coldata_file <- args[2] #"colData.csv"
base_group <- args[3] #"NC"
output_path <- args[4] 

#scriptPath="/data/script/deseq2+GO+KEGG/Rpipe/"
scriptPath="/Users/Icey/Documents/GitHub/R/deseq2-ver2/"
script1 = paste(scriptPath,"1deseq2-cor-heatmap.R",sep="")
script2 = paste(scriptPath,"2sigGene-volcano.R",sep="")


#1
path1 = paste(output_path,"3.DiffExprGene",sep="")
dir.create(path1,showWarnings = FALSE)
path3 = paste(output_path,"7.Correlation_analysis",sep="")
dir.create(path3,showWarnings = FALSE)
source(script1)


#2
path2 = path1
dir.create(path2,showWarnings = FALSE)
colData = read.csv(coldata_file, header=T)

for (exp_group in colData$condition[duplicated(colData$condition)])
  if (exp_group != base_group)
    source(script2)


# count_table <- "CountMatrix4DESeq.csv"
# coldata_file <- "colData.csv"
# base_group <- "BZ"
# dbname <- 'AH5955'
# kegg_org <- "vda"
# GO_KEY <- "SYMBOL"
# KEGG_NEED_KEY <- "SYMBOL"


