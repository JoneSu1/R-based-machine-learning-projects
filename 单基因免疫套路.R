setwd("F:/单基因套路/BJX-437-1")
########下载数据##############

rm(list=ls())  #清空环境内变量
options(stringsAsFactors = F)  #避免自动将字符串转换为R语言因子
library(GEOquery)
library(stringr)
library(ggplot2)
library(reshape2)
library(limma)
#使用GEOquery()函数来下载GEO数据然后使用exprs()函数获取表达矩阵.
gset = getGEO('GSE46517', destdir=".", AnnotGPL = F, getGPL = F)  #设置后两个为T比较耗时且作用不大
exp<-exprs(gset[[1]])  #exp即为表达矩阵
dim(exp)
head(exp)
#可以通过pData()函数来获得getGEO对象的临床信息。
#使用pData()函数提取临床信息
pdata<-pData(gset[[1]])
dim(pdata)
head(pdata)#
#write.csv(x = pdata,file = "pdata.csv")
#筛选出只含
#stringr包用于字符串的处理，str_detect是该包里的函数，
#用来确定一个字符向量能否匹配一种模式。它返回一个与输入向量具有同样长度的逻辑向量：
group<- pdata[,"tissue type:ch1"]
group<- cbind(group,pdata[,"tissue type:ch1"])
rownames(group)=group[,1]
group=group[,-1]
write.csv(x = group,file = "cl.csv")

#在excel中将转移癌和原发癌都定义为Tumor组，正常皮肤和痣定义为control，将正常组的去掉。
group1<- read.csv(file = "cl.csv",header = T,row.names = 1)
#统一样本顺序
exp=exp[,-105:-121]
dim(exp)
colnames(exp)
#进行重注释

gpl = getGEO("GPL96",destdir = ".")#这种取出的对象是s4
id_pre = gpl@dataTable@table#使用@调用s4对象中某个的值
colnames(id_pre)
head(id_pre)
ids = id_pre[,c("ID","Gene Symbol")]   #只需要"ID","Symbol"这两列
colnames(ids) = c("probe_id","symbol") #给这两列重新取名
tail(sort(table(ids$symbol)))#看每个基因用了些什么探针
head(ids)
#准备ID置换
#使得探针和symbol对应
table(rownames(exp) %in% ids$probe_id)#选取了exp中的行名和ids（s3）文件中的probe_id对应. table（）%in%（）的意思
#是进行逻辑判断a是不是在表b中。有26879个符合
exp = exp[rownames(exp) %in% ids$probe_id,]#把没有对应名字的过滤
table(rownames(exp) %in% ids$probe_id)#对比探针，探针数只剩下54675的了
dim(exp)
ids=ids[match(rownames(exp),ids$probe_id),]#使用match函数把ids里的探针顺序和表达矩阵的统一.
#match()函数返回的是一个位置向量，该向量记录着第一个参数中每个元素在第二个参数中的位置。所以，此时ids里的探针顺序与表达矩阵exp的探针顺序一一对应。
head(ids)
head(exp)#确定顺序对了
#完全对应上，我们就可以通过probe_id将表达矩阵exp进行分组，将同一个symbol所对应的多个探针分成不同的组，
#并对每组探针进行统计：计算每组中每行探针表达量的平均值（也就是每个探针在6个样本中表达量的均值rowMeans(x)），
#再取平均值最大的那个探针作为该symbol所对应的唯一探针，该组中的其它探针过滤掉，这样每个symbol就对应一个探针了：
tmp = by(exp,
         ids$symbol,
         function(x) rownames(x)[which.max(rowMeans(x))])
probes = as.character(tmp)
dim(exp)
exp = exp[rownames(exp) %in% probes,] # 过滤有多个探针的基因
dim(exp)
head(tmp)

#by()函数在这里发挥的功能就是将表达矩阵exp中的探针分组，同一个symbol所对应的多个探针分到一组，
#并对每组探针进行统计得到symbol所对应的唯一探针。所以tmp里放着by()函数的统计结果即每个symbol所对应的唯一探针IDprobe_id，
#用probes = as.character(tmp)将结果变身为纯字符型向量：
head(probes)
dim(exp)
dim(ids$symbol)
#by()函数就可以返回每个分组里的统计结果，即每个symbol所对应的唯一探针IDprobe_id。
#这时，探针ID和基因symbol就一一对应了，将表达矩阵探针ID即exp表达矩阵的行名（rownames(exp)）换为基因symbol。
rownames(exp)=ids[match(rownames(exp),ids$probe_id),2]
head(exp)
dim(exp)
exp=log2(exp+1)
write.csv(x = exp,file = "黑色素瘤exp.csv")
#计算每一列的平均表达值，再取这两列的表达值的中位数进行分类.5.974117
OVOL1=exp1[,"OVOL1"]
med.exp<-median(exp["OVOL1",])
OVOL1


more.med.exp.index<-which(OVOL1>=med.exp)
less.med.exp.index<-which(OVOL1<med.exp)
OVOL1[more.med.exp.index]<-paste0('High expression (',
                                        length(more.med.exp.index),')')
OVOL1[less.med.exp.index]<-paste0('Low expression (',
                                        length(less.med.exp.index),')')
table(OVOL1)
OVOL1=data.frame(OVOL1)
library(dplyr)
high= OVOL1%>%filter(OVOL1=="High expression (52)")
low=OVOL1%>%filter(OVOL1=="Low expression (52)")
Agroup= rbind(low,high)
table(Agroup$OVOL1)
#write.csv(x = Agroup,file = "OVOL1分组.csv")


#改顺序
class(exp)
colnames(exp)
rownames(Agroup)
dim(exp)
exp=exp[,rownames(Agroup)]
head(exp)
write.csv(x = exp,file = "高低分组的表达文件.csv")

#顺序统一了,Low52,High 52
group_list<- c(rep("Low",52),rep("High",52))
group_list=factor(group_list,
                  levels = c("Low","High"))

table(group_list)
#########进行差异分析##############


library(limma)
library(dplyr)
design=model.matrix(~factor(group_list))#把group设置成一个model matrix

colnames(design) <- c("Low", "High")
df.fit <- lmFit(exp, design)  ## 数据与list进行匹配

df.fit <- eBayes(df.fit)
tempOutput = topTable(df.fit,coef=2,n=Inf,adjust="BH")
nrDEG = na.omit(tempOutput)
mRNA.diff <- nrDEG

foldChange=0.5
padj=0.05

logFC <-mRNA.diff$logFC
deg.pad <- mRNA.diff$P.Value
p.j<- mRNA.diff$adj.P.Val
data <- data.frame(logFC=logFC,p=deg.pad,p.j=p.j)
data$mRNA.diffsig[(data$p > 0.05|data$p=="NA")|(data$logFC < foldChange)& data$logFC -foldChange] <- "no"
data$mRNA.diffsig[data$p <= 0.05 & data$logFC >= foldChange] <- "up"
data$mRNA.diffsig[data$p <= 0.05 & data$logFC <= -foldChange] <- "down"
x_lim <- max(logFC,-logFC)
rownames(data)<- rownames(mRNA.diff)
data["OVOL1",]
table(data$mRNA.diffsig)#down520, up562.186
deg=data
## 导出所有的差异结果
deg = na.omit(deg) ## 去掉数据中有NA的行或列

#为deg数据框添加几列
#1.加probe_id列，把行名变成一列
library(dplyr)
deg <- mutate(deg,probe_id=rownames(deg))
head(deg)

# write.csv(deg, "differ.csv")
#########绘制差异火山图#######
library(dplyr)
library(ggplot2)
dat  = deg


pdf(file = "low vs high_volcano.pdf",width=10,height=10)
dat  = deg
table(dat$mRNA.diffsig)

p <- ggplot(data = dat,
            aes(x = logFC,
                y = -log10(p))) +
  geom_point(alpha=0.4, size=3.5,
             aes(color=mRNA.diffsig)) +
  ylab("-log10(p_value)")+
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-0.5,0.5),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  theme_bw()
p <- p +theme(axis.text=element_text(size=20),axis.title=element_text(size=20))#是的labal变大

p
dev.off()
###热图

deg = na.omit(deg) ## 去掉数据中有NA的行或列


head(deg)
library(ggplot2)
x=deg$logFC
names(x)=deg$probe_id
cg = deg$probe_id[deg$mRNA.diffsig !="no"]
dim(exp)
table(rownames(exp)%in%cg)
n=exp[cg,]
dim(n)
library(pheatmap)
annotation_col=data.frame(group=group_list)
rownames(annotation_col)=colnames(n)
library(ggplotify)

pdf(file = "Low_vs_High_heatmap.pdf",width=10,height=10)
heatmap_plot<- as.ggplot(pheatmap(n,show_colnames =F,
                                  show_rownames = F,
                                  scale = "row",
                                  cluster_cols = F,
                                  annotation_col=annotation_col))


dev.off()

diffg= dat%>%filter(mRNA.diffsig!="no")
dim(diffg)
dim(dat)
##################差异基因富集分析##################
#利用clusterProfiler的内部函数enrichGo（）进行有参GO富集.
library(org.Hs.eg.db)#人类的H19
##clusterProfiler的 GO 富集
library(clusterProfiler)
library(stringr)
library(BiocGenerics)
require(DOSE)
library(ggplot2)
columns(org.Hs.eg.db)#查看数据库中有的类型
#准备文件

#一般来说一个基因可以对应多个GO terms，一个GO terms也可以对应多个基因。
#先将SYMBOL进行ID转换#=bitr()可以进行ID转换
x<- diffg$probe_id
eg= bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")#把id转换成ENTREZID
head(eg)
dim(eg)#967个基因有对应的ID
ego<-  enrichGO(
  gene  = eg$ENTREZID,
  keyType = "ENTREZID",
  OrgDb   = org.Hs.eg.db,
  ont     = "all",
  pAdjustMethod = "none",
  pvalueCutoff  = 0.05,
  qvalueCutoff = 1,
  readable      = TRUE)
dim(ego)#1263
head(ego)
ego_BP<-  enrichGO(
  gene  = eg$ENTREZID,
  keyType = "ENTREZID",
  OrgDb   = org.Hs.eg.db,
  ont     = "BP",
  pAdjustMethod = "none",
  pvalueCutoff  = 0.05,
  qvalueCutoff = 1,
  readable      = TRUE)
dim(ego_BP)#1042
head(ego_BP)
ego_CC<-  enrichGO(
  gene  = eg$ENTREZID,
  keyType = "ENTREZID",
  OrgDb   = org.Hs.eg.db,
  ont     = "CC",
  pAdjustMethod = "none",
  pvalueCutoff  = 0.05,
  qvalueCutoff = 1,
  readable      = TRUE)
dim(ego_CC)#101
head(ego_CC)

ego_MF<-  enrichGO(
  gene  = eg$ENTREZID,
  keyType = "ENTREZID",
  OrgDb   = org.Hs.eg.db,
  ont     = "MF",
  pAdjustMethod = "none",
  pvalueCutoff  = 0.05,
  qvalueCutoff = 1,
  readable      = TRUE)
dim(ego_MF)#120
head(ego_MF)
#write.csv(x = ego,file = "GO.csv") #储存GO信息表



#绘图
pdf(file = "GO.dot.pdf",width=8,height=7)
dotplot(ego, orderBy = "GeneRatio", color = "p.adjust",
        showCategory = 5, size = NULL,     font.size = 10, title = "GO Analysis",
        split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scales = "free")
dev.off()

#KEGG分析
#KEGG分析

ekegg <- enrichKEGG(
  gene = eg$ENTREZID,
  keyType = "kegg",
  organism  = 'hsa',#hsa是表示人
  pvalueCutoff  = 0.05,
  pAdjustMethod  = "none",
  qvalueCutoff = 1
)
head(ekegg)
dim(ekegg)#48
#绘图
pdf(file = "KEGG.dot.pdf",width=8,height=6)
dotplot(ekegg, orderBy='pvalue',decreasing=FALSE, showCategory = 10)
dev.off()
#write.csv(x = ekegg,file = "KEGG.csv")

#准备进行免疫浸润分析（sssGSVA）#############
setwd("F:/单基因套路/BJX-437-1/免疫")
#install.packages("pacman")
library(pacman)
p_load(tidyverse)
df <- read.csv("./data_set/store-reliab-data.csv", header = T) %>% as_tibble
geneset <- read.csv("immune_cell28.csv",header = F)%>% as_tibble #免疫细胞Marker基因文件从TISIDB数据库下载
#将免疫细胞的基因和名字分开，为了构建后续的列表
gene=read.table("genes.txt",header = F,sep = "\t",quote = "")
immune_cell=read.table("immune_cells.txt",header = F,sep = "\t",quote = "")
gene.set1=as.matrix(gene)
gene.set2=t(gene.set1)
dim(gene.set2)#77 28
gmt=list()#把基因构个网络
for (i in 1:28) {
  y=as.character(gene.set2[,i])
  b=which(y=="")#去掉列表中的空值
  gmt[[i]]=y[-b]
}
#给列表重新命名
names(gmt)=immune_cell$V1
View(gmt)
exp
geneset2=gmt
#准备完全了
#免疫文件准备好了，准备利用gsva来计算免疫细胞在不同样本的富集分数#####
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GSVA")
BiocManager::install("HDF5Array")
library(GSVA)
#ssGSEA分析/计算免疫细胞在不同样本的富集分数：
res <- gsva(
  as.matrix(exp), #表达量矩阵，可以是log-CPMs、log-RPKMs或log-TPMs
  geneset2, #Marker基因list
  method = "ssgsea", #单样本基因富集分数估算方法,可选"gsva", "ssgsea", "zscore"和"plage"
  kcdf = "Gaussian",
  mx.diff = F,
  verbose   = F #是否提供每个步骤详细信息
)

#Min-Max标准化：
#将每个样本中不同的免疫细胞富集分数标准化到0-1间；
resm <- res
for (i in colnames(res)) {
  resm[,i] <- (res[,i] -min(res[,i]))/(max(res[,i] )-min(res[,i] ))
}

#查看结果表(组织为矩阵)：
View(res)
write.csv(x = res,file = "immune_value.csv")
#绘制表达量直方图
##########直方图
library(RColorBrewer)
library(tidyr)
#添加样本注释(tumor和normal分组）：
group_list 
annotation <- data.frame(group_list)
rownames(annotation) <- colnames(res)
table(annotation)
head(annotation)
group2=c(rep("Low",52),rep("High",52))
resm=rbind(resm,group2)
group3=resm[28,]
resm=resm[-28,]
#宽数据转换为长数据(ggplot2绘图格式)：
dt <- resm %>% t() %>% as.data.frame() %>%
  rownames_to_column("sample") %>%
  gather(key = cell_type,
         value = value,
         -sample)
head(dt)

#将value根据样本转换为百分比形式(新增一列)：
dtt <- dt %>%
  group_by(sample) %>%
  mutate(proportion = round(value/sum(value),3))
head(dtt)

#新增分组列：
table(group_list)
dim(dtt)
group4=data.frame(group3[1:52])

Low=rownames(group4)
dtt$group <- ifelse(dtt$sample==Low,"Low","High")#利用ifelse函数赋值分组属性
##直方图尝试
library(RColorBrewer)
library(tidyr)
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))

ggplot(dtt,aes(sample,proportion,fill = cell_type)) + 
  geom_bar(stat = "identity") +
  labs(fill = "Cell Type",x = "",y = "Estiamted Proportion") + 
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = mypalette(27))

#分组箱线图展示：
p2 <- ggplot(dtt,
             aes(x = cell_type,y = proportion,fill = group)) +
  geom_boxplot(color = "black",alpha = 0.6,outlier.shape = 21,outlier.size = 1.2) +
  scale_fill_manual(values = c("#4979b6","#d9352a")) +
  labs(x = "cell type", y = "proportion") +
  theme_bw() + mytheme + theme(axis.text.x = element_text(angle = 45))
p2<- p2+stat_compare_means(method="t.test",hide.ns=F,comparisons=)

#重新指定箱线图排序（按相对丰度中位数从大到小）：
dtt_arrange <- dtt %>%
  group_by(cell_type) %>%
  summarise(de = median(proportion)) %>%
  arrange(desc(de)) %>%
  pull(cell_type)

dtt$cell_type <- factor(dtt$cell_type,levels = unique(dtt_arrange))


#ggpubr绘图：
library(ggpubr)
my_comparisons<-list(c("Low","High"))
p3 <- ggboxplot(dtt, x = "cell_type", y = "proportion",
                fill = "group", color = "black",
                width=0.7, alpha = 0.6,
                outlier.shape = 21, outlier.size = 1.2) +
  scale_fill_manual(values = c("#d9352a","#4979b6")) +
  labs(x = "cell type", y = "proportion") +
  theme_bw() + mytheme + theme(axis.text.x = element_text(angle = 75)) 
p3=p3+stat_compare_means(aes(group=group),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif")

pdf(file = "immune_boxplot.pdf",width = 10,height = 8)
print(p3)
dev.off()
####免疫细胞含量和目标基因的相关性
#准备免疫细胞打分文件和目标基因表达文件
#免疫细胞打分文件
resm

OVOL1=data.frame(exp["OVOL1",])
OVOL1$group=c(rep("Low",52),rep("High",52))
colnames(OVOL1)=c("OVOL1","group")
immu_data=t(resm)
immu_data1=immu_data[,c("Central memory CD4 T cell","Activated dendritic cell","Myeloid derived suppressor cell","Type 2 T helper cell","Memory B cell","Effector memeory CD4 T cell","Eosinophil","Neutrophil")]
colnames(immu_data1)
#8个免疫细胞有组间差异："Central memory CD4 T cell"            "Activated dendritic cell"           "Myeloid derived suppressor cell"               
#"Type 2 T helper cell"  "Memory B cell" "Effector memeory CD4 T cell"   "Eosinophil"   "Neutrophil" 


####准备免疫相关性分析和绘制棒棒糖图############
library(limma)
library(ggplot2)
library(ggpubr)
install.packages("ggExtra")
library(ggExtra)
library(dplyr)



# 打开两个想进行相关性分析的清洁数据
# 数据1为行名为样品名，列名为基因名
# 数据2为行名为样品名，列名为想做的和基因有相关性关系的，如免疫细胞

# 样品1:expr_data
# 样品2:immu_data

#进行相关性计算
y <- as.numeric(OVOL1$OVOL1)
#使用了R自带的cor.test包进行了计算
cor_data <- do.call(rbind,lapply(colnames(immu_data1),function(x){
  dd <- cor.test(as.numeric(immu_data1[,x]),y,method ="spearman",exact=FALSE)
  data.frame(cell=x,cor=dd$estimate,p.value=dd$p.value)
}))


# 画单基因与样本2中所有列名的相关性的棒棒糖图

#定义圆圈颜色的函数
p.col = c('gold','pink','orange','LimeGreen','darkgreen')
fcolor = function(x,p.col){
  color = ifelse(x>0.8,p.col[1],ifelse(x>0.6,p.col[2],ifelse(x>0.4,p.col[3],
                                                             ifelse(x>0.2,p.col[4], p.col[5])
  )))
  return(color)
}

#定义设置圆圈大小的函数
p.cex = seq(2.5, 5.5, length=5)
fcex = function(x){
  x=abs(x)
  cex = ifelse(x<0.1,p.cex[1],ifelse(x<0.2,p.cex[2],ifelse(x<0.3,p.cex[3],
                                                           ifelse(x<0.4,p.cex[4],p.cex[5]))))
  return(cex)
}

dat <- cor_data
gene <- 'OVOL1'
#根据pvalue定义圆圈的颜色
points.color = fcolor(x=dat$p.value,p.col=p.col)
dat$points.color = points.color

#根据相关系数定义圆圈的大小
points.cex = fcex(x=dat$cor)-0.8
dat$points.cex = points.cex
dat=dat[order(dat$cor),]

########绘制图形########
xlim = ceiling(max(abs(dat$cor))*10)/10         #x轴范围
pdf(file="Bangbangtang.pdf", width=9, height=7)      #输出图形
layout(mat=matrix(c(1,1,1,1,1,0,2,0,3,0),nc=2),width=c(8,2.2),heights=c(1,2,1,2,1))
par(bg="white",las=1,mar=c(5,18,2,4),cex.axis=1.5,cex.lab=2)
plot(1,type="n",xlim=c(-xlim,xlim),ylim=c(0.5,nrow(dat)+0.5),xlab="Correlation Coefficient",ylab="",yaxt="n",yaxs="i",axes=F)
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col="#F5F5F5",border="#F5F5F5")
grid(ny=nrow(dat),col="white",lty=1,lwd=2)
#绘制图形的线段
segments(x0=dat$cor,y0=1:nrow(dat),x1=0,y1=1:nrow(dat),lwd=4)
#绘制图形的圆圈
points(x=dat$cor,y = 1:nrow(dat),col = dat$points.color,pch=16,cex=dat$points.cex)+
  scale_size_continuous(range =c(2,4))
#展示免疫细胞的名称
text(par('usr')[1],1:nrow(dat),dat$cell,adj=1,xpd=T,cex=1.5)
#展示pvalue
pvalue.text=ifelse(dat$p.value<0.001,'<0.001',sprintf("%.03f",dat$p.value))
redcutoff_cor=0
redcutoff_pvalue=0.05
text(par('usr')[2],1:nrow(dat),pvalue.text,adj=0,xpd=T,col=ifelse(abs(dat$cor)>redcutoff_cor & dat$p.value<redcutoff_pvalue,"red","black"),cex=1.5)
axis(1,tick=F)

#绘制圆圈大小的图例
par(mar=c(0,4,3,4))
plot(1,type="n",axes=F,xlab="",ylab="")
legend("left",legend=c(0.1,0.2,0.3,0.4,0.5),col="black",pt.cex=p.cex,pch=16,bty="n",cex=2,title="abs(cor)")

#绘制圆圈颜色的图例
par(mar=c(0,6,4,6),cex.axis=1.5,cex.main=2)
barplot(rep(1,5),horiz=T,space=0,border=NA,col=p.col,xaxt="n",yaxt="n",xlab="",ylab="",main="pvalue")
axis(4,at=0:5,c(1,0.8,0.6,0.4,0.2,0),tick=F)
dev.off()
########进行免疫治疗分析#########
#先使用write.txt去输出以Table作为分隔符号的文件
# 这里为了得到比较好的结果，采用two direction median centered
#使用原始探针数据试试
setwd("F:/单基因套路/BJX-437-1/免疫")
exp00=data.frame(exp)
write.csv(x = exp00,file = "exp.csv")
TIDE=read.csv(file = "exp.csv",header = T)
#计算行平均值，按降序排列
index=order(rowMeans(TIDE[,-1]),decreasing = T)
#调整表达谱的基因顺序
expr_ordered=TIDE[index,]
#对于有重复的基因，保留第一次出现的那个，即行平均值大的那个
keep=!duplicated(expr_ordered$X)
#得到最后处理之后的表达谱矩阵
expr_max=expr_ordered[keep,]
expr_max
write.csv(x = expr_max,file = "exp2.csv")
TIDE=read.csv(file = "exp2.csv",header = T)
rownames(TIDE)=TIDE$X
TIDE=TIDE[,-1:-2]
TIDE <- sweep(TIDE, 2, apply(TIDE, 2, median, na.rm = T))
TIDE <- sweep(TIDE, 1, apply(TIDE, 1, median, na.rm = T))
write.table(TIDE,"TIDE_input.self_subtract",sep = "\t",row.names = T,col.names = NA,quote = F)
# 3 统计分析#######################
# 参照文件夹中TIDE使用教程得到输出文件TIDE_output.csv
TIDE.res <- read.csv("immune_treat.csv",header = T,row.names = 1,check.names = F,stringsAsFactors = F)
ann=Agroup
ann <- merge(ann,TIDE.res,by="row.names")
print(table(ann$Responder,ann$OVOL1))

# 检验免疫治疗响应性和亚型是否相关，p<0.05表示相关
print(fisher.test(table(ann$Responder,ann$OVOL1))) 


#画图
## 3.1堆叠柱状图#########################
library(ggplot2)
# pdf('histogram_TIDE.pdf',height = 6,width = 6)
R=colnames(ann)
colnames(ann)=c("Row.names","group",
                "No benefits" ,"Responder","TIDE","IFNG","MSI Expr Sig"
                ,"Merck18","CD274","CD8","CTL.flag","Dysfunction","Exclusion","MDSC","CAF","TAM M2") 
pdf(file = "堆叠图免疫治疗.pdf",width = 10,height = 8)
ggplot(data=ann,mapping=aes( x=group,fill=Responder))+
  geom_bar()+
  theme_bw()+
  scale_y_continuous(expand = c(0,0),limits=c(0,55))+ #(起始值，终止值，间隔),调整y轴属性，使柱子与X轴坐标接触,需要根据我们的样本量进行修改
  theme(
    #标题字体设置
    text=element_text(size=12),
    #设置标题居中hjust = 0.5
    plot.title = element_text(hjust = 0.5,vjust = 0.5), 
    #Y轴字体设置
    axis.text.y=element_text(size=12,color = "black"),
    #X轴字体设置
    #angle：调整横轴标签倾斜角度
    #hjust：上下移动横轴标签
    axis.text.x=element_text(size=12,  color = "black", hjust = 0.5,
                             vjust = 0.5),
    #图例的标题字体设置，可以修改colour、size
    #legend.title=element_text(size=12), 
    #图例字体设置
    legend.text=element_text(size=12)
    #legend.position = ' none' #删除图例
    #legend.position="bottom" ,#图例位置放中间，可选参数为“left”,“top”, “right”, “bottom”.
    #legend.background = element_rect(fill="lightblue",size=0.5, linetype="solid", colour ="darkblue")
  )
dev.off()


## 3.2 TIDE在高低风险组的差异###################

library(ggpubr)
sz=15

my_comparisons<-list(c("high","low"))

g<-ggboxplot(ann, x="group", y="TIDE", color = "group",
             ylab="TIDE",  xlab="group ", palette = 'nejm',add = "jitter",size = 1,axis.line =2)+
  theme(axis.text.y=element_text(size=sz),
        axis.title=element_text(size=sz),
        axis.text.x=element_text(size=sz)) +
  theme(plot.margin=unit(rep(3,4),'lines'),
        legend.text = element_text(size = sz),
        legend.title = element_text(size = sz)) +
  
  stat_compare_means(label.y =12,method="wilcox.test",size =4)+
  stat_compare_means(comparisons=my_comparisons,label="p.signif")
pdf(file = "TIDE在高低表达组中的差异.pdf",width = 10,height = 8)
g
dev.off()

# rotate_x_text(60)
ggsave(g,filename = "ImmuneTherapy_RiskGroup_boxplot.pdf",width = 10,height = 8)


##进行免疫治疗的submap分析，submap预测亚型的免疫治疗响应性
#输入表达矩阵和样品信息，整理成submap所需的格式，在线运行后输出SubMap_SubMapResult.txt用于画图。
# 自定义函数用来产生submap需要的数据格式
setwd("F:/单基因套路/BJX-437-1/免疫")
sbumap=read.csv(file = "exp.csv",header = T,row.names = 1)
in_gct=sbumap
generateInputFileForSubMap <- function(in_gct, gct_file, cls_file, sam_info, type_name = "type"){
  in_gct <- data.frame(GeneID=rownames(in_gct),
                       description="na",
                       in_gct, 
                       stringsAsFactors = F,
                       check.names = F)
  cat("#1.2\n", file = gct_file)
  cat(nrow(in_gct),"\t",ncol(in_gct)-2,"\n", file = gct_file, append = T)
  cat(paste(colnames(in_gct), collapse = "\t"),"\n", file = gct_file, append = T)
  for(i in 1:nrow(in_gct)) cat(paste(in_gct[i,], collapse = "\t"),"\n", file = gct_file, append = T)
  
  cat(nrow(sam_info),length(levels(factor(sam_info$rank))),1, "\n", file = cls_file )
  cat("#", paste0(levels(factor(sam_info[, type_name])), collapse = " " ), "\n", file = cls_file, sep = "", append = T)
  cat(as.numeric(factor(sam_info[, type_name])), file = cls_file, append = T)
}

# 创建submap需要的数据格式 (SKCM)
skcm.immunotherapy.logNC <- read.table("skcm.immunotherapy.47samples.log2CountsNorm.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F) #原文提供的log2转化的标准化count值
rownames(skcm.immunotherapy.logNC) <- toupper(rownames(skcm.immunotherapy.logNC)) # 基因大写，因为我使用的数据是把基因名都大写的
skcm.immunotherapy.info <- read.table("skcm.immunotherapy.47sampleInfo.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)

skcm.immunotherapy.info <- skcm.immunotherapy.info[order(skcm.immunotherapy.info$label),]
skcm.immunotherapy.info$rank <- rep(c(1,2,3,4),times=as.character(table(skcm.immunotherapy.info$label))) #1: CTLA4_noR 2: CTLA4_R 3:PD1_noR 4:PD1_R


# 创建submap需要的数据格式 (TCGA)
exp
dim(exp)
exp=exp[,rownames(Agroup)]
write.csv(x = exp,file = "exp3.csv")
tmp <- read.csv(file = "exp3.csv",header = T)
#去除重复基因
#计算行平均值，按降序排列
index=order(rowMeans(tmp[,-1]),decreasing = T)
#调整表达谱的基因顺序
expr_ordered=tmp[index,]
#对于有重复的基因，保留第一次出现的那个，即行平均值大的那个
keep=!duplicated(expr_ordered$id)
#得到最后处理之后的表达谱矩阵
expr_max=expr_ordered[keep,]
rownames(expr_max)=expr_max$id
expr_max=expr_max[,-1]
tmp=expr_max
#去除平均表达量<1 90%的原始数据
tmp <- tmp[apply(tmp, 1, function(x) sum(x > 1) > 104*0.9), ]
dim(tmp)
head(data.frame(exp["UBC",]))
head(data.frame(tmp["UBC",]))
identical(exp["ACTB",],tmp["ACTB",])


GENELIST <- intersect(rownames(tmp),rownames(skcm.immunotherapy.logNC)) # 取交集后的基因列表

sam_info <- skcm.immunotherapy.info
in_gct <- skcm.immunotherapy.logNC[GENELIST,rownames(skcm.immunotherapy.info)]

# 产生输出数据的文件名
gct_file <- "skcm.immunotherapy.for.SubMap.gct"
cls_file <- "skcm.immunotherapy.for.SubMap.cls"
generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")

# 提出亚型的样本，顺序排列
samples.C1 <- rownames(ann[which(ann$ImmClust == "C1"),])
samples.C2 <- rownames(ann[which(ann$ImmClust == "C2"),])

sam_info <- data.frame("ImmClust"=c(samples.C1,samples.C2),row.names = c(samples.C1,samples.C2))
sam_info$rank <- rep(c(1,2),times=c(length(samples.C1),length(samples.C2))) #1: C1,即HPV16-IMM 2: C2,即HPV16-KRT

# 产生输出数据的文件名
gct_file <- "Immune2.for.SubMap.gct"
cls_file <- "Immune2.for.SubMap.cls"

in_gct <- log2(tmp[GENELIST,rownames(sam_info)] + 1) # 产生和示例数据类似的形式，log2转化的标准化count值
generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")
