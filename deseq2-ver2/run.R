#!/usr/bin/env Rscript

args<-commandArgs(T)

# test if there is at least one argument: if not, return an error
if (length(args)!=5) {
	print ("5 arguments must be supplied:")
	print ("matrixFile groupFile baseGroup outputPath have_replicates?")
	print ("CountMatrix4DESeq.csv colData.csv NC ./ TRUE")
  stop(call.=FALSE)
} 


count_table <- args[1] #"CountMatrix4DESeq.csv"
coldata_file <- args[2] #"colData.csv"
base_group <- args[3] #"NC"
output_path <- args[4]
have_replicates <- args[5]

#scriptPath="/data/script/deseq2+GO+KEGG/Rpipe/"
scriptPath="/Users/Icey/Documents/GitHub/R/deseq2-ver2/"
script1 = paste(scriptPath,"1deseq2-cor-heatmap.R",sep="")
script2_1 = paste(scriptPath,"2.1deseq2-volcano.R",sep="")
script2_2 = paste(scriptPath,"2.2edgeR-volcano.R",sep="")

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



if (have_replicates){
for (exp_group in colData$condition[duplicated(colData$condition)]){
  if (exp_group != base_group) {source(script2_1)}}
} else {
for (exp_group in colData$condition){
  if (exp_group != base_group) {source(script2_2)}
}}



