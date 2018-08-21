#!/usr/bin/env Rscript

args<-commandArgs(T)

# test if there is at least one argument: if not, return an error
if (length(args)!=5) {
  print ("5 arguments must be supplied:")
  print ("GeneKeyType pSet ModelOrNot KEGGDir SigFileDir")
  print ("ENSEMBL 0.05 FALSE ./nonModel/4.GO_KEGG_Enrichment/ ./nonModel/3.DiffExprGene/")
  stop(call.=FALSE)
} 


###Arguments
GENE_KEY= args[1] #"ENSEMBL" 
PSET = args[2] #0.05
MODEL = args[3] #FALSE
KEGG_DIR = args[4] #"./nonModel/4.GO_KEGG_Enrichment/"
SIG_PATH = args[5] #"./nonModel/3.DiffExprGene/"

suppressMessages(library(pathview))
#############################################
######For Model Animal(eg.Human)#############
#############################################
PathViewModel <- function(sigfile,pSet){  
  #folddir
  groups = sapply(strsplit(sigfile, "_sig"), "[", 1)
  pathway_path = paste(KEGG_DIR,"/KEGG_Pathway","/",groups,"/",sep="")
  
  #check if exists，exists then delete
  if dir.exists(pathway_path) {unlink(pathway_path,recursive=TRUE}
  dir.create(pathway_path,showWarnings = FALSE,recursive=T)
  
  #read in data
  Gene <- read.table(paste(SIG_PATH,sigfile,sep="/"), sep="\t", header = T, row.names = 1)        
  path_info = paste(KEGG_DIR,pSet,"/",groups,"_KEGG_Enrichment.txt",sep="")
  Path <- read.table(path_info, header = T, row.names = 1, sep="\t")
  
  #Gene expr data
  Gene$geneID <- rownames(Gene)
  
  #去掉括号及括号内字符
  Gene$geneID = gsub("\\s*\\([^\\)]+\\)","",Gene$geneID)
  
  #去除重复
  Gene <- Gene[!duplicated(Gene$geneID),]
  row.names(Gene) <- Gene$geneID
  Gene4Path <- dplyr::select(Gene,log2FoldChange)
  
  #KEGG path id
  KEGG_path_id <- (rownames(Path))
  setwd(pathway_path)

  cc <- function(KEGG_path_id){
    pv.out <- pathview(gene.data = Gene4Path[,1,drop=FALSE], pathway.id = KEGG_path_id, 
                       species = "hsa", kegg.native = TRUE,
                       gene.idtype = GENE_KEY)
  }
  
  sapply(KEGG_path_id, cc) 
  setwd("../../../")
  
}

#############################################
######For Non-Model Animal##################
#############################################
PathViewNoModel <- function(sigfile,pSet){  
  #folddir
  groups = sapply(strsplit(sigfile, "_sig"), "[", 1)
  pathway_path = paste(KEGG_DIR,"/KEGG_Pathway","/",groups,"/",sep="")
  dir.create(pathway_path,showWarnings = FALSE,recursive=T)
  
  #read in data
  Gene <- read.table(paste(SIG_PATH,sigfile,sep="/"), sep="\t", header = T)        
  path_info = paste(KEGG_DIR,pSet,"/",groups,"_KEGG_Enrichment.txt",sep="")
  Path <- read.table(path_info, header = T, row.names = 1, sep="\t")
  
  #colnames
  colnames(Gene)[1:4] = c("koID","geneID","baseMean","log2FoldChange")
  
  #去除重复
  Gene <- Gene[!duplicated(Gene$koID),]
  row.names(Gene) <- Gene$koID
  Gene4Path <- dplyr::select(Gene,log2FoldChange)
  
  #KEGG path id
  KEGG_path_id <- (rownames(Path))
  setwd(pathway_path)
  cc <- function(KEGG_path_id){
    pv.out <- pathview(gene.data = Gene4Path[,1,drop=FALSE], pathway.id = KEGG_path_id, 
                       species = "ko", kegg.native = TRUE,
                       gene.idtype = GENE_KEY)
  }
  
  sapply(KEGG_path_id, cc) 
  setwd("../../../")
  
}


#############################################
######Run Run Run############################
#############################################


if (MODEL == TRUE) {
  sigfiles = list.files(SIG_PATH,pattern="sig_genes_exprData.txt")
  sapply(sigfiles, PathViewModel, pSet=PSET)
}

if (MODEL == FALSE) {
  sigfiles = list.files(SIG_PATH,pattern="sig_pathview.txt")
  sapply(sigfiles, PathViewNoModel, pSet=PSET)
}
