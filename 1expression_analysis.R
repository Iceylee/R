

# load count table and coldata
count_table = "CountMatrix_F_genename.csv"
coldata_file = "colData.csv"


# Load the library.
library(DESeq2)
library(dplyr)


##########sep NEED CHANGE#################
countData = read.table(count_table, header=TRUE, sep=" ", row.names=1)
View(countData)

colData = read.csv(coldata_file, header=T, sep=",", row.names=1 )
View(colData)

colData$condition = as.factor(colData$condition)
colnames(countData) <- rownames(colData)
#check
#all(rownames(colData) %in% colnames(countData))

# Create DESEq2 dataset.
dds = DESeqDataSetFromMatrix(countData=countData, colData=colData, design = ~condition)

#pre filter
dds = dds[ rowSums(counts(dds)) > 1 ,] 

#Set the reference to be compared
####NEED CHANGE#######
dds$condition = relevel(dds$condition,"BZ")


##
# Run deseq
dds = DESeq(dds)

# Format the results.
res = results(dds)

# Sort the results data frame by the padj and foldChange columns.
sorted = res[with(res, order(padj, -log2FoldChange)), ]

# Turn it into a dataframe to have proper column names.
sorted.df = data.frame("id"=rownames(sorted),sorted)

# Write the table out.
write.table(sorted.df, file="deseq2.txt", row.names = FALSE,sep="\t", quote=FALSE)

# Get normalized counts and write this to a file
nc = counts(dds,normalized=TRUE)

# Turn it into a dataframe to have proper column names.
dt = data.frame("id"=rownames(nc),nc)

# Save the normalize data matrix.
write.table(dt, file="norm-matrix-deseq2.txt", sep="\t",  row.name=FALSE, col.names=TRUE,quote=FALSE)

########significantly different genes
#padj（<0.05)和log2 fold (>1)


regSig <- subset(res, padj < 0.05)
regSig2 <- subset(regSig, abs(log2FoldChange) > 2)

sorted_regSig2 = regSig2[with(regSig2, order(-log2FoldChange)), ]

sig = data.frame("id"=rownames(regSig2),regSig2)
write.table(sig, file="sig-genes-deseq2.txt", sep="\t",  row.name=FALSE, col.names=TRUE,quote=FALSE)


save(dds, res, sorted.df, dt, sig, file="1Data_Input.RData")


