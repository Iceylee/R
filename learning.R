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
foo <- data.frame()
#2.
within(df, FOO<-data.frame(do.call('rbind', strsplit(as.character(FOO), '|', fixed=TRUE))))
#3.
require(reshape)
df <- data.frame(ID=11:13, FOO=c('a|b','b|c','x|y'))
df = transform(df, FOO = colsplit(FOO, split = "\\|", names = c('a', 'b')))

