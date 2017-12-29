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

rep(seq(1, 7, by = 2), times = 7)

#正则
# The emails vector has already been defined for you
emails <- c("john.doe@ivyleague.edu", "education@world.gov", "dalai.lama@peace.org",
            "invalid.edu", "quant@bigdatacollege.edu", "cookie.monster@sesame.tv")
# Use grepl() to match for "edu"   [1]  TRUE  TRUE FALSE  TRUE  TRUE FALSE
grepl("edu",emails)
# Use grep() to match for "edu", save result to hits  [1] 1 2 4 5
grep("edu",emails)   
# Subset emails using hits
emails[grepl("edu",emails)]
emails[grep("edu",emails)]
#匹配@
#匹配任何字符多次 .*
#反义圆点 \\
grepl("@.*\\.edu$",emails)
# Use sub() to convert the email domains to datacamp.edu
sub("@.*\\.edu$","@datacamp.edu",emails)

#报错的情况赋值
qobj <- tryCatch(qvalue(p=Over$pvalue, lambda=0.05, pi0.method="bootstrap"), error=function(e) NULL)

# Print the R version details using version
version
# Assign the variable `major` to the major component
major <-  version$major
# Assign the variable `minor` to the minor component
minor <- version$minor

# How long does it take to read movies from CSV?
system.time(read.csv("movies.csv"))
# How long does it take to read movies from RDS?
system.time(readRDS("movies.rds"))
#等号和箭头的区别 在函数内部使用箭头，变量也被赋值。 
median(x <- 1:10); x
#比较函数运行的时间
# Load the microbenchmark package
library(microbenchmark)
# Compare the two functions，10即各运行10次
compare <- microbenchmark(read.csv("movies.csv"), 
                          readRDS("movies.rds"), 
                          times = 10)
# Print compare
print(compare)

##function
#for循环不好 如果df是空 将会报错
df <- data.frame()
1:ncol(df)

for (i in 1:ncol(df)) {
  print(median(df[[i]]))
}
##用seq_along
for (i in seq_along(df)) {
  print(median(df[[i]]))
}
##将df的结果存在output里面
# Create new double vector: output
output = vector("double", length = ncol(df))
#
for (i in seq_along(df)) {
  output[[i]] = median(df[[i]])
}
# Print output
print (output)

##
# Define example vectors x and y
x <- c( 1, 2, NA, 3, NA)
y <- c(NA, 3, NA, 3,  4)
# Count how many elements are missing in both x and y
# x和y都是NA的，计1.sum求这样的个数
sum(is.na(x) & is.na(y))
