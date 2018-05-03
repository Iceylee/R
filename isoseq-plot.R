setwd("/Users/Icey/work/pacbio三代/iso-seq全长转录组/有参全长转录组/original data/")
library(ggplot2)
curve = read.table(file="saturation_curve.txt", header=T,stringsAsFactors = F)
options(scipen=200) #avoid 1e10
p <- ggplot(data=curve, aes(x=FL_reads_number, y=Isoform_number, group=1)) +
  geom_line()+
  scale_x_continuous(breaks=seq(0,200000,20000))+
  scale_y_continuous(breaks=seq(0,120000,10000)) 
pdf(file="reads_isofrom.pdf")
p
dev.off()

p <- ggplot(data=curve, aes(x=FL_reads_number, y=Gene_number, group=1)) +
  geom_line()+
  scale_x_continuous(breaks=seq(0,200000,20000))+
  scale_y_continuous(breaks=seq(0,10000,1000)) 
pdf(file="reads_gene.pdf")
p
dev.off()
  


#splice
as = read.table(file="AS_count.out", header=F,stringsAsFactors = F)
as = as[order(as$number, decreasing = TRUE),]



colnames(as) <- c("AS_event","number")
#myLabel = as.vector(as$AS_event)   
#myLabel = paste(myLabel, "(", round(as$number / sum(as$number) * 100, 2), "%)", sep = "")  
myLabel = paste(round(as$number / sum(as$number) * 100, 2),"%",sep = "")

p <- ggplot(as, aes(x = "", y = number, fill = AS_event)) + 
  geom_bar(stat = "identity",width =1) + 
  theme_bw()+
  theme(axis.ticks = element_blank())+theme(panel.grid=element_blank()) +    ## 去掉白色圆框和中间的坐标线
  theme(panel.border=element_blank()) +
  theme(axis.text.x = element_blank()) +
  geom_text(aes(y = number/3 + c(0, cumsum(number)[-length(number)]),x=1,label = myLabel), size = 4) +
  coord_polar(theta = "y") +
  labs(x = "", y = "", title = "") 


ggsave("splice.as.pdf",p, width=7, height=5, units="in")

#geom_text(aes(y = A/2 + c(0, cumsum(A)[-length(A)]), x = sum(A)/20, label = myLabel), size = 5)   ## 在图中加上百分比：x 调节标签到圆心的距离, y 调节标签的左右位置 