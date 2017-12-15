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

#修改levels顺序，按照表格的顺序（逆序），而不是字母表。作图的轴线出现的顺序会受影响	
ego_three$Description<- factor(ego_three$Description, order=TRUE, levels=ego_three$Description)

p <- p + theme(panel.spacing = unit(0, "lines"))

ggplot
p+theme(axis.text = element_text(size = 8，face="bold"), legend.text = element_text(size = 8)) 





###可以按照表格显示的term顺序作图 ---序列是从小到大的
##如何让表格显示的顺序是BP CC MF呢（level需要修改）
ego_three$term<- factor(ego_three$term, order=TRUE, levels=rev(ego_three$term))
# ego_three$class<- factor(ego_three$class, order=TRUE, levels=c("MF","CC","BP"))



#分开string并计数
p = c('a|b|c|d','b|c|c','x|y')
pp = strsplit(as.character(p),'|',fixed=TRUE)
count = unlist(lapply(pp,length))


##分栏
#1.
p = c('a|b|','b|c','x|y')
pp = strsplit(as.character(p),'|',fixed=TRUE)
ppp = do.call('rbind',pp ) #function name; a list -output a list 函数只能接受两个参数，list中上一个和当前。
foo <- data.frame(ppp)
#2.
within(df, FOO<-data.frame(do.call('rbind', strsplit(as.character(FOO), '|', fixed=TRUE))))
#3.
require(reshape)
df <- data.frame(ID=11:13, FOO=c('a|b','b|c','x|y'))
df = transform(df, FOO = colsplit(FOO, split = "\\|", names = c('a', 'b')))

#
#引号中为一个object，该object是否存在
exists("ego_CC_df")


###研究
#诉求：三个df，检测是否存在。rbind合并存在的几个df
#ego_list <- list(BP=ego_BP_df,CC=ego_BP_CC,MF=ego_BP_MF)
ego_list <- c("ego_BP_df", "ego_MF_df", "ego_CC_df")
TorF <- sapply(ego_list,exists)
real_list <- sapply(ego_list[TorF],as.name)
ppp = do.call('rbind',real_list )


##remove colon:
gene_id <- sapply(gene_id,function (x) {unlist(strsplit(x,split=":"))[2]})


##变量1转列名 变量2转行名 统计变量3
library(reshape2)
rpkm_df <- read.table(file = "AllSamplesRPKMValue.txt",sep = '\t',header =F,stringsAsFactors=F)
dff <- rpkm_df[,1:3] #只能有三列
rpkm_tab <- dcast(dff, V2 ~ V1) #V2为左侧列 V1为最上的行
#转回去是melt和cast


