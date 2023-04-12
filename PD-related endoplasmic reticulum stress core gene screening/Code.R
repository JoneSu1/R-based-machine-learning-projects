##############进行差异分析（differ）##############

rm(list=ls()) #清空环境内变量 (Clear all variables in the environment)
options(stringsAsFactors = F) #避免自动将字符串转换为R语言因子 (Prevent automatic conversion of strings to R factors)
suppressMessages(library(GEOquery))
library(stringr)
library(ggplot2)
library(reshape2)
library(limma)
library("GEOquery")
##########下载数据集并转换ID，准备clinical文件(Download the dataset and convert the ID, prepare the clinical file)###########
setwd("E:/427-1WGCNA-Diff-ROC诊断/test")
#这个样本作为测试集，50个pd（帕金森），33个非PD，23个完全健康 (This sample is used as a test set, 50 PD (Parkinson's), 33 non-PD, 23 completely healthy)
#使用GEOquery()函数来下载GEO数据然后使用exprs()函数获取表达矩阵. (Use the GEOquery() function to download GEO data and then use the exprs() function to obtain the expression matrix.)
gset = getGEO('GSE6613', destdir=".", AnnotGPL = F, getGPL = F) #设置后两个为T比较耗时且作用不大 (Setting the latter two to T takes more time and has little effect)
exp<-exprs(gset[[1]]) #exp即为表达矩阵 (exp is the expression matrix)
#进行重注释 (Perform re-annotation)

gpl = getGEO('GPL96',destdir = ".")#这种取出的对象是s4 (The extracted object is s4)
id_pre = gpl@dataTable@table#使用@调用s4对象中某个的值 (Use @ to call a value in an s4 object)
colnames(id_pre)
head(id_pre)
#先都转换成GB_ACC(Convert all to GB_ACC first)
ids2 = id_pre[,c("ID","Gene Symbol")] #只需要"ID","GB_ACC"这两列 (Only need "ID" and "GB_ACC" columns)
colnames(ids2) = c("probe_id","symbol") #给这两列重新取名 (Rename these two columns)
tail(sort(table(ids2$symbol)))#看每个基因用了些什么探针 (See what probes each gene used)
ids=as.data.frame(ids2)
table(rownames(exp) %in% ids$probe_id)#选取了exp中的行名和ids（s3）文件中的probe_id对应. table（）%in%（）的意思 (Select the row names in exp and the probe_id in the ids (s3) file. The meaning of table() %in% ())
exp = exp[rownames(exp) %in% ids$probe_id,]#把没有对应名字的过滤 (Filter out those without corresponding names)
table(rownames(exp) %in% ids$probe_id)#对比探针，探针数只剩下26183的了 (Compare probes, only 26183 probes left)
ids=ids[match(rownames(exp),ids$probe_id),]#使用match函数把ids里的探针顺序和表达矩阵的统一. (Use the match function to unify the probe order in ids with the expression matrix.)
head(ids)
head(exp)#确定顺序对 (Confirm the order is correct)
tmp = by(exp,
ids$symbol,
function(x) rownames(x)[which.max(rowMeans(x))])
probes = as.character(tmp)
dim(exp)
exp = exp[rownames(exp) %in% probes,] # 过滤有多个探针的基因 (Filter genes with multiple probes)
dim(exp)
head(tmp)
rownames(exp)=ids[match(rownames(exp),ids$probe_id),2]
head(exp)

#准备临床数据 (Prepare clinical data)
#使用pData()函数提取临床信息 (Use the pData() function to extract clinical information)
pdata<-pData(gset[[1]])
dim(pdata)
head(pdata)#
cl=data.frame(pdata[,1])
rownames(cl)= rownames(pdata)

#write.csv(x = cl,file = "cl.csv")
#write.csv(x = pdata,file = "pdata.csv")
group=read.csv(file = "cl.csv",header = T,row.names = 1)
#统一样本顺序 (Unify the sample order)
table(colnames(exp)%in%rownames(group))
table(group$Type)
exp=exp[,rownames(group)]
dim(exp)#22 Control 50,Disease

####准备差异分析############## (Prepare differential analysis)
group_list<-c(rep("Control",22),rep("Disease",50))%>% factor(., levels = c("Control", "Disease"), ordered = F)
table(group_list)
exp=log2(exp+1)
class(exp)

#write.csv(x = exp,file = "exp.csv")
library(limma)
library(dplyr)

design<-model.matrix(~0+group_list)
colnames(design)<-c("Control","Disease")
fit<-lmFit(exp,design)
contrast.matrix <- makeContrasts(Disease-Control,levels=design)
fit2<-contrasts.fit(fit,contrast.matrix)
fit2<-eBayes(fit2)
result<-topTable(fit2,coef=1,adjust="BH",number=dim(exp) [1])
mRNA.diff=result
mRNA.diff["UBC",]
exp["UBC",]

logFC <-mRNA.diff$logFC
deg.padj <- mRNA.diff$P.Value
data <- data.frame(logFC=logFC,padj=deg.padj)
data$mRNA.diffsig[(data$padj > 0.05|data$padj=="NA")|(data$logFC < foldChange)& data$logFC -foldChange] <- "no"
data$mRNA.diffsig[data$padj <= 0.05 & data$logFC >= foldChange] <- "up"
data$mRNA.diffsig[data$padj <= 0.05 & data$logFC <= -foldChange] <- "down"
x_lim <- max(logFC,-logFC)
rownames(data)<- rownames(mRNA.diff)
head(data)
data["MOB1A",]#和原表达数据验证 UP (Verify with the original expression data UP)
exp["MOB1A",]
###做过基因箱图看差异 (Look at the differences with gene boxplots)
#先做箱图（Do the box diagram first）




ANKRD11<- exp[c("ANKRD11","UBC"),]
ANKRD11<-data.frame(t(ANKRD11))
ANKRD11$ID=rownames(ANKRD11)
ANKRD11=ANKRD11[,-2]
ANKRD11=ANKRD11[,c("ID","ANKRD11")]
ANKRD11$Type=group$Type
#绘图（plotting）

library(ggpubr)
p <- ggboxplot(ANKRD11, x = "Type", y = "ANKRD11",
               fill = "Type",legend=F,palette =c("#00AFBB", "#E7B800"),bxp.errorbar=T)+
  theme(legend.position='none')+
  ylab(label = 'MOB1A expression')

my_comparisons <- list( c("Disease", "Control") )
p + stat_compare_means(comparisons = my_comparisons)

#验证完成 (Validation complete)
table(data$mRNA.diffsig)#down117, up158 no13241
colnames(data)<- c("logFC","P_value","diffsig")

deg <- mutate(data,probe_id=rownames(data))
deg=na.omit(deg)
#绘制火山图 (Draw volcano plot)
diff=deg%>%select("diffsig")%>%filter(diffsig!="no")
dim(diff)
write.csv(x = diff,file = "diff.csv")

#########绘制差异火山图####### (Draw differential volcano plot)
library(dplyr)
library(ggplot2)
dat = deg

pdf(file = "volcano.pdf",width=8,height=6)
dat = deg

p <- ggplot(data = dat,
aes(x = logFC,
y = -log10(P_value))) +
geom_point(alpha=0.4, size=3.5,
aes(color=diffsig)) +
ylab("-log10(Pvalue)")+
scale_color_manual(values=c("blue", "grey","red"))+
geom_vline(xintercept=c(-0.5,0.5),lty=4,col="black",lwd=0.8) +
geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
theme_bw()
p<- p +theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20))
p
dev.off()
###热图 (Heatmap)


x=deg$logFC
names(x)=deg$probe_id
cg = deg$probe_id[deg$diffsig !="no"]
n=exp[cg,]
dim(n)
library(pheatmap)
annotation_col=data.frame(group=group_list)
rownames(annotation_col)=colnames(n)
library(ggplotify)
#列不聚类 (No column clustering)
pdf(file = "heatmap.pdf",width=8,height=5)
heatmap_plot <- as.ggplot(pheatmap(n,show_colnames =F,
show_rownames = F,
scale = "row",
cluster_cols = F,
annotation_col=annotation_col,
main = "differ_heatmap"))
options(repr.plot.width = 3, repr.plot.height =2)#调整图大小 (Adjust the image size)
heatmap_plot
dev.off()
#########WCGNA分析，TCGA的数据############# (WGCNA analysis, TCGA data)
#准备进行无标度拓扑分析的文件 (Prepare files for scale-free topology analysis)
library(WGCNA)
gene<- read.csv(file = "TCGA_Standardization_exp.csv",header = T,row.names = 1)
dim(gene)
#过滤低表达的，只保留平均表达值在 1 以上的,为过滤掉,15369 (Filter low-expression, keep only average expression value above 1, filtered out, 15369)
gene <- subset(gene, rowSums(gene)/ncol(gene) >= 1)
dim(gene)
#转置矩阵，行是基因，列是样本 (Transpose matrix, rows are genes, columns are samples)
gene<- t(gene)
dim(gene)
datExpr_filted <- gene
#检查离群样品#找Plot函数里面的字体，调小一点.cex= 用来调节字体大小。 (Adjust font size with cex)
pdf(file = "TCGAcluster.pdf",width=8,height=4)
gsg = goodSamplesGenes(datExpr_filted, verbose = 3)
gsg$allOK
sampleTree = hclust(dist(datExpr_filted), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",cex = 0.05)
abline(h = 320, col = 'red')#给聚类图加红线#样本无离群 (Add red line to the clustering diagram, no outliers in the samples)
dev.off()
#确定软阈值(soft threholding prowers)作为幂指数β (Determine the soft threshold as the power exponent β)
#pickSoftThreshold()可通过将给定的基因表达矩阵按表达相似度计算加权网络，并分析指定powers值（构建一组向量，
#传递给powerVector）下网络的无标度拓扑拟合指数以及连通度等信息，目的是帮助为网络构建选择合适的powers值。 (pickSoftThreshold() can calculate the weighted network according to the expression similarity of the given gene expression matrix and analyze the scale-free topology fitting index and connectivity of the network under the specified powers value (construct a set of vectors, pass to powerVector), the purpose is to help select the appropriate powers value for network construction.)
##power 值选择 (Select power value)
powers <- 1:20
sft <- pickSoftThreshold(gene, powerVector = powers, verbose = 5)
#拟合指数与 power 值散点图 (Fitting index and power value scatter plot)

pdf(file = "TCGApower_point.pdf",width=8,height=6)
par(mfrow = c(1, 2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], type = 'n',
     xlab = 'Soft Threshold (power)', ylab = 'Scale Free Topology Model Fit,signed R^2',
     main = paste('Scale independence'))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels = powers, col = 'red');
abline(h = 0.85, col = 'red')

#平均连通性与 power 值散点图 (Mean connectivity and power value scatter plot)
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab = 'Soft Threshold (power)', ylab = 'Mean Connectivity', type = 'n',
main = paste('Mean connectivity'))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, col = 'red')
dev.off()
sft#通过pickSoftThreshold函数估计的软阈值为8. (Soft threshold estimated by pickSoftThreshold function is 8)
powers <- sft$powerEstimate+1
#获得 TOM 矩阵 (Obtain TOM matrix)
#使用了TOMsimilarity()函数进行计算 (TOMsimilarity() function is used for calculation)
adjacency <- adjacency(gene, power = powers)
tom_sim <- TOMsimilarity(adjacency)
rownames(tom_sim) <- rownames(adjacency)
colnames(tom_sim) <- colnames(adjacency)
tom_sim[1:6,1:6]
#输出 TOM 矩阵 (Output TOM matrix)
#write.table(tom_sim, 'TOMsimilarity.txt', sep = '\t', col.names = NA, quote = FALSE)
#共表达模块识别 (Co-expression module identification)
#上述获得了基因间的共表达相似度矩阵，记录了基因间表达的相似性。
#可以进一步据此计算距离并绘制聚类树，查看相似度的整体分布特征 (The above obtained the co-expression similarity matrix between genes, which records the similarity of gene expression. Based on this, the distance can be further calculated and the clustering tree can be plotted to view the overall distribution characteristics of similarity)
#TOM 相异度 = 1 – TOM 相似度 (TOM dissimilarity = 1 - TOM similarity)
tom_dis <- 1 - tom_sim
#层次聚类树，使用中值的非权重成对组法的平均聚合聚类 (Hierarchical clustering tree, using the average agglomerative clustering of the non-weighted pair-group method with median)

#层次聚类树，使用中值的非权重成对组法的平均聚合聚类 (Hierarchical clustering tree, using the average agglomerative clustering of the non-weighted pair-group method with median)

geneTree <- hclust(as.dist(tom_dis), method = 'average')
plot(geneTree, xlab = '', sub = '', main = 'Gene clustering on TOM-based dissimilarity',
labels = FALSE, hang = 0.04)

#使用动态剪切树挖掘模块 (Use dynamic tree cutting to mine modules)
minModuleSize <- 100 #模块基因数目,基因样本量大可以把模块设大到100或者200. (Module gene number, if the gene sample size is large, the module can be set to 100 or 200)
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = tom_dis,
deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize)

table(dynamicMods)
#模块颜色指代 12 (Module color代表 12 个模块 (Represents 12 modules)
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)

plotDendroAndColors(geneTree, dynamicColors, 'Dynamic Tree Cut',
dendroLabels = FALSE, addGuide = TRUE, hang = 0.03, guideHang = 0.05,
main = 'Gene dendrogram and module colors')
dev.off()
#模块特征基因计算 (Module eigengene calculation)
#使用moduleEigengenes()可对基因表达值矩阵中的基因按模块划分后，执行特征分解，获取模块特征基因。 (Using moduleEigengenes(), genes in the gene expression matrix can be divided into modules, perform eigenvalue decomposition, and obtain module eigengenes)
#计算基因表达矩阵中模块的特征基因（第一主成分） (Calculate the eigengenes of the modules in the gene expression matrix (first principal component))
MEList <- moduleEigengenes(gene, colors = dynamicColors)
MEs <- MEList$eigengenes
head(MEs)[1:6]
#输出模块特征基因矩阵 (Output module eigengene matrix)
#write.table(MEs, 'moduleEigengenes.txt', sep = '\t', col.names = NA, quote = FALSE)
#共表达模块的进一步聚类 (Further clustering of co-expression modules)
#通过模块特征基因计算模块间相关性，表征模块间相似度 (Calculate the correlation between modules through the module eigengenes, characterizing the similarity between modules)
ME_cor <- cor(MEs)
ME_cor[1:6,1:6]
#绘制聚类树观察 (Plot the clustering tree to observe)
METree <- hclust(as.dist(1-ME_cor), method = 'average')
plot(METree, main = 'Clustering of module eigengenes', xlab = '', sub = '')
#探索性分析，观察模块间的相似性 (Exploratory analysis, observe the similarity between modules)
#height 值可代表模块间的相异度，并确定一个合适的阈值作为剪切高度 (The height value can represent the dissimilarity between modules and determine an appropriate threshold as the cutting height)
#以便为低相异度（高相似度）的模块合并提供依据 (To provide a basis for merging modules with low dissimilarity (high similarity))
abline(h = 0.2, col = 'blue')
abline(h = 0.3, col = 'red')
#模块特征基因聚类树热图 (Module eigengene clustering tree heatmap)
plotEigengeneNetworks(MEs, '', cex.lab = 0.8, xLabelsAngle= 90,
marDendro = c(0, 4, 1, 2), marHeatmap = c(3, 4, 1, 2))
#相似模块合并，以 0.25 作为合并阈值（剪切高度），在此高度下的模块将合并 (Merge similar modules with 0.25 as the merge threshold (cutting height), modules at this height will be merged)
#近似理解为相关程度高于 0.75 的模块将合并到一起 (It can be approximately understood that modules with a correlation higher than 0.75 will be merged together)
merge_module <- mergeCloseModules(gene, dynamicColors, cutHeight = 0.3, verbose = 3)
mergedColors <- merge_module$colors
table(mergedColors)
#基因表达和模块聚类树,重剪切后保留了11个模块，原来有12个模块 (Gene expression and module clustering tree, after re-cutting, 11 modules are retained, originally there were 12 modules)

pdf(file = "TCGA_genes_clustertree.pdf",width=10,height=8)

plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c('Dynamic Tree Cut', 'Merged dynamic'),
                    dendroLabels = FALSE, addGuide = TRUE, hang = 0.03, guideHang = 0.05)
dev.off()

#共表达模块与与临床表型的关联分析 (Co-expression module and clinical phenotype association analysis)
#患者的临床表型数据 (Patient's clinical phenotype data)
trait <- read.csv(file = "TCGA_GROUP.list.csv",header = T,row.names = 1)
dim(trait)

#使用上一步新组合的共表达模块的结果 (Using the co-expression module results from the previous step)
module <- merge_module$newMEs
dim(module)
dim(trait)
rownames(module)
rownames(trait)
table(rownames(module)%in%rownames(trait))

#转换trait里面的符号 (Convert symbols in trait)
rownames(trait) <- substring(rownames(trait),1,15) %>% gsub("-",".",.)#把影响输入的符号改成点 (Change the symbols affecting the input to dots)
table(rownames(module)%in%rownames(trait))
#使得临床数据和表达数据相和，不和的活cor函数会有NA。 (Make clinical data and expression data match, and cor function will have NA for mismatched.)
#使得探针和symbol对应 (Make probe and symbol correspond)
ids2<- module

ids=as.data.frame(ids2)
table(rownames(trait)%in% rownames(ids))#选取了exp中的行名和ids（s3）文件中的probe_id对应. table（）%in%（）的意思 (Select the row names in exp and the probe_id in the ids (s3) file. The meaning of table() %in% ())
datTraits<- trait
colnames(datTraits)<- c("Normal","Tumor")
head(datTraits)
#患者基因共表达模块和临床表型的相关性分析 (Patient gene co-expression module and clinical phenotype correlation analysis)

step 5 (最重要的) 模块和性状的关系 (The most important step 5, module and trait relationship)
这一步主要是针对于连续变量，如果是分类变量，需要转换成连续变量方可使用 (This step is mainly for continuous variables. If it is a categorical variable, it needs to be converted into a continuous variable for use)
#明确样本数和基因数 (Specify the number of samples and genes)
nGenes = ncol(gene)
nSamples = nrow(gene)

identical(rownames(trait),row.names(gene))
dim(module)
module=module[rownames(trait),]
dim(module)
moduleTraitCor <- cor(module, trait, use = 'p')
moduleTraitCor #相关矩阵 (Correlation matrix)
#相关系数的 p 值矩阵 (P value matrix of correlation coefficient)
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(module))
#相关图绘制 black,yellow,black (Correlation plot drawing black, yellow, black)
textMatrix <- paste(signif(moduleTraitCor, 2), '\n(', signif(moduleTraitPvalue, 1), ')', sep = '')
dim(textMatrix) <- dim(moduleTraitCor)
sizeGrWindow(6,6)#设置图片size的 (Set picture size)
par(mar = c(5, 10, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,main = paste('Module-trait relationships'),
xLabels = names(trait), yLabels = names(module), ySymbols = names(module),
colorLabels = FALSE, colors = greenWhiteRed(50), cex.text = 0.6, zlim = c(-1,1),
textMatrix = textMatrix, setStdMargins = FALSE)
###############输出################### (Output)
#site<-which(abs(as.numeric(moduleTraitCor[,2]))>=0.40)
site<-which(abs(as.numeric(moduleTraitPvalue[,1]))<1)#输出全部 (Output all)
site3<-gsub("ME","",rownames(moduleTraitCor))
site3<-site3[site]
moduleTraitCor<-moduleTraitCor[site,]
moduleTraitPvalue<-moduleTraitPvalue[site,]
setwd("/data/project-Sujunhua/YQS-419-12/YQS-419-12/TCGA_OUTPUT")

for (mod in site3)
{
modules = mod
probes = colnames(gene)
inModule = (mergedColors == modules)#这里的mergedColors是未经过切割的原颜色文件 (Here, mergedColors is the original color file that has not been cut)
modGenes = probes[inModule]
write.table(modGenes, file =paste0("./",modules,".txt"),sep=",",row.names=F,col.names=F,quote=F)
}
table(moduleTraitCor)

##############数据集结果不好，用新的GEO数据集，用DESEq2包处理#################### (The dataset results are not good, use a new GEO dataset, and process with DESEq2 package)
rm(list=ls()) #清空环境内变量 (Empty variables in the environment)
options(stringsAsFactors = F) #避免自动将字符串转换为R语言因子 (Avoid automatically converting strings to R language factors)
suppressMessages(library(GEOquery))
library(stringr)
library(ggplot2)
library(reshape2)
library(limma)
library("GEOquery")
setwd("E:/YQS-427-1WGCNA-Diff-ROC诊断/test/脑的数据集")
exp<- read.csv(file = "exp_count.csv",header = T,row.names = 1)#使用的数据集是GSE68719 (The dataset used is GSE68719)
group=read.csv(file = "cl.csv")
dim(exp)#74 17580
table(group$Type)#Control Diesease

44 29
#读取来之后就是count数据，先除去0多的。 (After reading, it is count data, first remove the excess zeros.)
#我们可以对这一部分在大多数样本中count为0或≤1的geneID进行过滤。关于筛选没有非常严格的标准， (We can filter out geneIDs with count of 0 or ≤1 in most samples. There is no very strict standard for screening,)
#保留在至少75%的样本中基因count数都大于1的gene。 (Retain genes with gene count greater than 1 in at least 75% of the samples.)
exp1 <- exp[apply(exp, 1, function(x) sum(x > 1) > 74*0.75), ]
dim(exp1)#16808 74
#把有小数点的去了 (Remove the ones with decimals)
#先对symbol去重 (First, remove duplicatesfrom symbol)
#计算行平均值，按降序排列 (Calculate row averages and sort them in descending order)
expr=exp1
index=order(rowMeans(expr[,-1]),decreasing = T)
#调整表达谱的基因顺序 (Adjust the gene order of the expression profile)
expr_ordered=expr[index,]
#对于有重复的基因，保留第一次出现的那个，即行平均值大的那个 (For duplicate genes, keep the one that appears first, i.e., the one with a larger row average)
keep=!duplicated(expr_ordered$symbol)
#得到最后处理之后的表达谱矩阵 (Get the final processed expression profile matrix)

exp1 <- round(exp1) #对表达量矩阵取整数 (Round the expression matrix)
#准备进行差异分析 (Prepare for differential analysis)
#count的数据需要用DEseq2 来做差异. (Count data needs DEseq2 for differential analysis)
Agroup_list<-c(rep("Control",44),rep("Disease",29)) %>% factor(., levels = c("Control", "Disease"), ordered = F)

library(DESeq2)
library(dplyr)
#创建表达量矩阵样本名和分组一一对应的数据框（colData）(Create a data frame (colData) with expression matrix sample names and groups one-to-one)
colData <- data.frame(row.names =colnames(exp1),
condition=Agroup_list)
head(colData)
colnames(exp1)
rownames(colData)
#用函数DESeqDataSetFromMatrix()构建DESeqDataSet对象（dds），包括存储输入值、计算中间值、标准化处理等 (Use the DESeqDataSetFromMatrix() function to build the DESeqDataSet object (dds), including storing input values, computing intermediate values, and standardization processing)
dds <- DESeqDataSetFromMatrix(
countData = exp1,
colData = colData,
design = ~ condition)

dds
#使用函数DESeq()进行差异分析流程：计算差异倍数、P值等 (Use the DESeq() function for differential analysis workflow: calculate differential multiples, P values, etc.)
dds <- DESeq(dds)
dds

#使用函数results()从dds中提取差异分析结果表：给出样本间的基本均值、log2FC、标准误差、检验统计量、p 值和调整后 p 值；(Use the results() function to extract the differential analysis results table from dds: give the basic averages between samples, log2FC, standard error, test statistics, p value, and adjusted p value;)
res <- results(dds,
contrast = c("condition",rev(levels(Agroup_list))))
#转化为数据框： (Convert to data frame:)
RES<- as.data.frame(res)
head(RES)
write.csv(x = RES,file = "ALL_diff.csv")
#保存差异分析数据： (Save differential analysis data:)
save(dds,res,RES,file = c("DESeq2.Rdata"))
#将数据进行标准化 (Standardize the data)
vsd <- vst(dds, blind=FALSE) #VST
assay(vsd)[1:4,1:4]
vsdd <- as.data.frame(assay(vsd))
#保存标准化数据： (Save standardized data:)

write.csv(x = vsdd,file = "Standardization_exp.csv")
#标记up和down (Mark up and down)

foldChange=0.5
padj=0.05

logFC <-RES$log2FoldChange
deg.padj <- RES$pvalue
data <- data.frame(logFC=logFC,padj=deg.padj)
data$mRNA.diffsig[(data$padj > 0.05|data$padj=="NA")|(data$logFC < foldChange)& data$logFC -foldChange] <- "no"
data$mRNA.diffsig[data$padj <= 0.05 & data$logFC >= foldChange] <- "up"
data$mRNA.diffsig[data$paddata$mRNA.diffsig[data$padj <= 0.05 & data$logFC <= -foldChange] <- "down"
x_lim <- max(logFC,-logFC)
rownames(data)<- rownames(RES)
head(data)
table(data$mRNA.diffsig)

data['GAS7',] #和原表达数据验证 DOWN (Verify with original expression data DOWN)
vsdd['GAS7',]

data['AOC1',] #和原表达数据验证UP (Verify with original expression data UP)
vsdd['AOC1',30:50]
GAS7 <- vsdd['GAS7',]
bar_mat <- t(GAS7)
group
bar_mat <- as.data.frame(bar_mat)
bar_mat <- cbind(group,bar_mat)
library(RColorBrewer)
library(ggpubr)
colnames(bar_mat) <- c("sample","sam","GAS7")
bar_mat$sam <- factor(bar_mat$sam,levels=c("Control","Diesease"))
color <- c("#5CB85C","#337AB7")
my_comparisons <- list(c("Control","Diesease"))

#开始设置循环#为了保证正常输出所以在目标基因列表下添加了其他基因 (Start setting loop #In order to ensure normal output, other genes have been added under the target gene list)
p <- ggboxplot(bar_mat, x = "sam", y = "GAS7",
color = "sam", palette = "jco") +
rotate_x_text(angle = 90) #将x轴肿瘤名称旋转90°展示 (Rotate the tumor names on the x-axis 90° for display)

#验证完成 (Verification completed)
table(data$mRNA.diffsig) #down1314, up1829 no18055
deg=data

导出所有的差异结果 (Export all differential results)
deg = na.omit(deg) ## 去掉数据中有NA的行或列 (Remove rows or columns with NA in the data)

#为deg数据框添加几列 (Add a few columns to the deg data frame)
#1.加probe_id列，把行名变成一列）(Add probe_id column and turn row names into a column)



library(dplyr)
deg <- mutate(deg,probe_id=rownames(deg))
head(deg)
diff=deg%>%select(logFC,padj,mRNA.diffsig,)%>%filter(mRNA.diffsig!="no")
dim(diff)
# write.csv(diff, "differ.csv")
#write.csv(deg, "differ_all_notice.csv")
#################差异绘图（plotting differ）################

library(dplyr)
library(ggplot2)
exp=vsdd
dat  = deg



pdf(file = "volcano.pdf",width=8,height=6)
dat  = deg

p <- ggplot(data = dat,
            aes(x = logFC,
                y = -log10(padj))) +
  geom_point(alpha=0.4, size=3.5,
             aes(color=mRNA.diffsig)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-0.5,0.5),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  theme_bw()
p<- p +theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20))
p
dev.off()
###热图（heatmap）



x=deg$logFC
names(x)=deg$probe_id
cg = deg$probe_id[deg$mRNA.diffsig !="no"]
n=exp[cg,]
dim(n)
library(pheatmap)
annotation_col=data.frame(group=Agroup_list)
rownames(annotation_col)=colnames(n)
library(ggplotify)
#列不聚类
pdf(file = "heatmap.pdf",width=8,height=5)
heatmap_plot <- as.ggplot(pheatmap(n,show_colnames =F,
                                   show_rownames = F,
                                   scale = "row",
                                   cluster_cols = F,
                                   annotation_col=annotation_col,
                                   main = "differ_heatmap"))
options(repr.plot.width = 3, repr.plot.height =2)#调整图大小
heatmap_plot
dev.off()



###################用交集基因做富集分析##################（Enrichment analysis with intersecting genes）
#利用clusterProfiler的内部函数enrichGo（）进行有参GO富集. (Use the internal function enrichGo() of clusterProfiler to perform reference-based GO enrichment)
library(org.Hs.eg.db) #人类的H19 (Human H19)
##clusterProfiler的 GO 富集 (GO enrichment of clusterProfiler)
library(clusterProfiler)
library(stringr)
library(BiocGenerics)
require(DOSE)
library(ggplot2)
columns(org.Hs.eg.db) #查看数据库中有的类型 (Check the types available in the database)
#准备文件 (Prepare files)

#一般来说一个基因可以对应多个GO terms，一个GO terms也可以对应多个基因。
# (In general, a gene can correspond to multiple GO terms, and a GO term can also correspond to multiple genes)
#先将SYMBOL进行ID转换#=bitr()可以进行ID转换 (First, convert the SYMBOL to ID using bitr() for ID conversion)

setwd("E:/YQS-427-1WGCNA-Diff-ROC诊断/test/脑的数据集/WCGNA")
x <- read.csv(file = "Diff_WGCNA_Endoplasmic_genes.csv")
dim(x) #44
x <- x$ID
eg = bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db") #把id转换成ENTREZID (Convert id to ENTREZID)
head(eg)
dim(eg)
dim(eg) #44个相关基因中仅有24个具有ENTREZID。则就用具有ENTREZID的进行富集.更改计算方法non，不改没结果）
#(Among the 44 related genes, only 24 have ENTREZID. Then use the ones with ENTREZID for enrichment. Change the calculation method to non, 
#without changing there will be no result)
ego<-  enrichGO(
  gene  = eg$ENTREZID,
  keyType = "ENTREZID",
  OrgDb   = org.Hs.eg.db,
  ont     = "all",
  pAdjustMethod = "none",
  pvalueCutoff  = 0.05,
  qvalueCutoff = 0.3,
  readable      = TRUE)
dim(ego)#1497
ego_BP<-  enrichGO(
  gene  = eg$ENTREZID,
  keyType = "ENTREZID",
  OrgDb   = org.Hs.eg.db,
  ont     = "BP",
  pAdjustMethod = "none",
  pvalueCutoff  = 0.05,
  qvalueCutoff = 0.3,
  readable      = TRUE)
dim(ego_BP)#1332
head(ego_BP)
ego_CC<-  enrichGO(
  gene  = eg$ENTREZID,
  keyType = "ENTREZID",
  OrgDb   = org.Hs.eg.db,
  ont     = "CC",
  pAdjustMethod = "none",
  pvalueCutoff  = 0.05,
  qvalueCutoff = 0.3,
  readable      = TRUE)
dim(ego_CC)#77
head(ego_CC)

ego_MF<-  enrichGO(
  gene  = eg$ENTREZID,
  keyType = "ENTREZID",
  OrgDb   = org.Hs.eg.db,
  ont     = "MF",
  pAdjustMethod = "none",
  pvalueCutoff  = 0.05,
  qvalueCutoff = 0.3,
  readable      = TRUE)
dim(ego_MF)#88
head(ego_MF)
#write.csv(x = ego,file = "downGO.csv") #储存GO信息表（save GO information to file）
#write.csv(x = ego_BP,file = "downego_BP.csv")
#write.csv(x = ego_CC,file = "downego_CC.csv")
#write.csv(x = ego_MF,file = "downego_MF.csv")


#绘图
pdf(file = "GO.dot.pdf",width=8,height=10)
dotplot(ego, orderBy = "GeneRatio", color = "p.adjust",
        showCategory = 5, size = NULL,     font.size = 9, title = "GO Analysis",
        split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scales = "free")
dev.off()
######进行分开绘图，再合并########(be first plotting single figure then merge )
p1<- dotplot(
  ego_BP,
  x = "GeneRatio",
  color = "p.adjust",
  title = "Top 10 of GO BP terms Enrichment",
  showCategory = 10,
  label_format = 30,
  font.size = 10
)
p2<- dotplot(
  ego_CC,
  x = "GeneRatio",
  color = "p.adjust",
  title = "Top 10 of GO CC terms Enrichment",
  showCategory = 10,
  label_format = 30,
  font.size = 10
)
p3<- dotplot(
  ego_MF,
  x = "GeneRatio",
  color = "p.adjust",
  title = "Top 10 of GO MF terms Enrichment",
  showCategory = 10,
  label_format = 30,
  font.size = 10
)
pdf(file = "downGO.dot.pdf",width=8,height=10)

library(patchwork)
p1/p2/p3
#KEGG分析
#KEGG分析
dev.off()
ekegg <- enrichKEGG(
  gene = eg$ENTREZID,
  keyType = "kegg",
  organism  = 'hsa',#hsa是表示人(hsa means human)
  pvalueCutoff  = 0.05,
  pAdjustMethod  = "none",
  qvalueCutoff = 1
)
head(ekegg)
dim(ekegg)#92
#绘图(plot)
pdf(file = "KEGG.dot.pdf",width=8,height=4)
dotplot(ekegg, orderBy='pvalue',decreasing=FALSE, showCategory = 10)
dev.off()
#write.csv(x = ekegg,file = "KEGG.csv")

##################进行ROC验证#################(ROC analysis)
library(pROC)
library(ROCR)
library(ggplot2)
library(tidyverse)

#进行ROC验证前先把基因上下调的给分出来 (Before performing ROC validation, separate the upregulated and downregulated genes)

setwd("E:/YQS-427-1WGCNA-Diff-ROC诊断/test/脑的数据集/WCGNA")
list = read.csv(file = "Diff_WGCNA_Endoplasmic_genes.csv")
list2 = diff[list$ID,]
dim(list2)
ROC_up = list2 %>% select(mRNA.diffsig) %>% filter(mRNA.diffsig == "up")
dim(ROC_up) #57个上调 (57 upregulated)
ROC_down = list2 %>% select(mRNA.diffsig) %>% filter(mRNA.diffsig == "down")
dim(ROC_down) #11个下调 (11 downregulated)
######先计算上调的基因的AUC值################# (First, calculate the AUC value of upregulated genes)
ROC_up1 <- vsdd[rownames(ROC_up),]
dim(vsdd)
dim(ROC_up1)
dim(group)
ROC_up1 = t(ROC_up1)
table(group$Type)
ROC_up1 = cbind(group, ROC_up1)
rt = ROC_up1
rt = rt[, -1]
#先计算看看有多少基因是大于0.8的 (First, calculate how many genes have an AUC value greater than 0.8)

计算type为结局变量的基因的相关信息 (Calculate the correlation information of genes with the outcome variable "type")
name = colnames(rt)
name = name[-1]

result <- rt %>% roc(response = Type, predictor = CLU,
levels = c("Control", "Diesease"), direction = "<") #指定计算poor状态的AUC (Calculate the AUC for poor condition)
result["auc"]

L = list()
for (i in 1:57) { result[[i]] <- rt %>% roc_(response = "Type", predictor = name[i],
levels = c("Control", "Diesease"), direction = "<") #指定计算Diesease状态的AUC (Calculate the AUC for Diesease condition)
L[[i]] = result[[i]]["auc"]

}
L[[1]]

用sapply将每一个元素变成向量，组成一个矩阵。下面两种写法等价 (Use sapply to convert each element into a vector and form a matrix. The following two methods are equivalent)
data.frame(t(sapply(L, c)))
T = data.frame(t(sapply(L, [))) #将刚从上面计算出来的AUC值转换成data.frame格式 (Convert the AUC value just calculated above to data.frame format)
colnames(T) = name
T = t(T)
colnames(T) = "AUC_value"
class(T)
T = data.frame(T)
dim(T)
#筛选出>0.8的,有15个基因是大于等于0.8的 (Filter out AUC values > 0.8, there are 15 genes with AUC values greater than or equal to 0.8)
AUC <- T %>% select(AUC_value) %>% filter(AUC_value >= 0.8)
dim(AUC)
AUC = data.frame(AUC)
class(AUC)
up_AUC = AUC
write.csv(x = T, file = "up_genes_AUC_great_than.csv")

########使用下调基因的来进行计算##### (Calculate using downregulated genes)
ROC_down1 <- vsdd[rownames(ROC_down),]
dim(vsdd)
dim(ROC_down1)
dim(group)
ROC_down1 = t(ROC_down1)
table(group$Type)
ROC_down1 = cbind(group, ROC_down1)
rt = ROC_down1
rt = rt[, -1]
#先计算看看有多少基因是大于0.8的 (First, calculate how many genes have an AUC value greater than 0.8)

计算type为结局变量的基因的相关信息 (Calculate the correlation information of genes with the outcome variable "type")
name = colnames(rt)
name = name[-1]

result <- rt %>% roc(response = Type, predictor = CLU,
levels = c("Control", "Diesease"), direction = "<") #指定计算poor状态的AUC (Calculate the AUC for poor condition)
result["auc"]

L = list()
for (i in 1:length(name)) { result[[i]] <- rt %>% roc_(response = "Type", predictor = name[i],
levels = c("Control", "Diesease"), direction = "<") #指定计算Diesease状态的AUC (Calculate the AUC for Diesease condition)
L[[i]] = result[[i]]["auc"]

}
L[[1]]

用sapply将每一个元素变成向量，组成一个矩阵。下面两种写法等价 (Use sapply to convert each element into a vector and form a matrix. The following two methods are equivalent)
data.frame(t(sapply(L, c)))
T = data.frame(t(sapply(L, [))) #将刚从上面计算出来的AUC值转换成data.frame格式 (Convert the AUC value just calculated above to data.frame format)
colnames(T) = name
T = t(T)
colnames(T) = "AUC_value"
class(T)
T = data.frame(T)
dim(T)
#筛选出>0.8的,有多少个基因是大于等于0.8的 (Filter out AUC values > 0.8, how many genes have AUC values greater than or equal to 0.8)
AUC <- T %>% select(AUC_value) %>% filter(AUC_value >= 0.8)
dim(AUC)
AUC = data.frame(AUC)
class(AUC)
down_AUC = AUC
write.csv(x = T, file = "down_genes_AUC_great_than.csv")



#####################画单图确认结果##################################
plot(result,
lwd=5, #曲线线宽 (Curve linewidth)
lty=1, #线条形状（0-6）(Line shape (0-6))
print.auc=TRUE, # 图像上输出AUC的值 (Output AUC value on the image)
print.auc.x=0.4, print.auc.y=0.5, # 设置AUC值坐标 (Set AUC value coordinates)
auc.polygon=TRUE, # 将ROC曲线下面积转化为多边形 (Convert the area under the ROC curve to a polygon)
auc.polygon.col="#fff7f7", # ROC曲线下填充色 (Fill color under the ROC curve)
grid=c(0.5, 0.2), # 设置两轴网格线的间隔 (Set the interval of grid lines for both axes)
grid.col=c("black", "black"), # 设置两轴间隔线条的颜色 (Set the color of the interval lines for both axes)
print.thres=TRUE, # 图像上输出最佳截断值 (Output the best cutoff value on the image)
main="ROC curves of CLU", # 添加图形标题 (Add graph title)
col="#FF2E63", # 设置ROC曲线颜色 (Set ROC curve color)
legacy.axes=TRUE)

setwd("E:/YQS-427-1WGCNA-Diff-ROC诊断/test/脑的数据集/WCGNA/output/up_output")

######################## 出down的图 (Plot downregulated genes)
list=rbind(up_AUC,down_AUC) #17个 (17 items)
rt=vsdd[rownames(down_AUC),]
dim(rt) #2
rt=t(rt)
rt=cbind(group,rt)
rt=rt[,-1]
dat <- pivot_longer(rt,-Type,names_to = "Index", #设置指标列列名 (Set index column name)
values_to = "expression") #设置数值列列名 (Set value column name)
head(dat)

A tibble: 6 x 3
#Type Index expression
#<chr> <chr> <dbl>
#1 Control CLU 17.0
#2 Control HSPA1A 11.3
#3 Control HSPB1 12.6
#4 Control VIM 12.7
#5 Control CD74 12.2
#6 Control G3BP2 11.9
#构建循环 (Construct loop)
setwd("E:/YQS-427-1WGCNA-Diff-ROC诊断/test/脑的数据集/WCGNA/output/down_output")
result <- unique(dat$Index);result
lapply(result, function(x){
file_me=paste(x,".pdf",sep = "")
pdf(file_me,width = 8,height = 5)
plot(roc(dat %>% filter(Index == x),response = Type, predictor =expression,
levels=c("Control", "Diesease"),direction = ">"),
lwd=5,
lty=1,
print.auc=TRUE,
print.auc.x=0.4, print.auc.y=0.5,
auc.polygon=TRUE,
auc.polygon.col="#fff7f7",
grid=c(0.5, 0.2),
grid.col=c("black", "black"),
print.thres=TRUE,
main= paste(x,"_ROC curves",sep = ""),
col="#FF2E63",
legacy.axes=TRUE)
dev.off()
})

出up的图 (Plot upregulated genes)
list=rbind(up_AUC,down_AUC) #17个 (17 items)
rt=vsdd[rownames(up_AUC),]
dim(rt) #15
rt=t(rt)
rt=cbind(group,rt)
rt=rt[,-1]
dat <- pivot_longer(rt,-Type,names_to = "Index", #设置指标列列名 (Set index column name)
values_to = "expression") #设置数值列列名 (Set value column name)
head(dat)

A tibble: 6 x 3
#Type Index expression
#<chr> <chr> <dbl>
#1 Control CLU 17.0
#2 Control HSPA1A 11.3
#3 Control HSPB1 12.6
#4 Control VIM 12.7
#5 Control CD74 12.2
#6 Control G3BP2 11.9
#构建循环 (Construct loop)



setwd("E:/YQS-427-1WGCNA-Diff-ROC诊断/test/脑的数据集/WCGNA/output/up_output")
result <- unique(dat$Index);result
lapply(result, function(x){
  file_me=paste(x,".pdf",sep = "")
  pdf(file_me,width = 8,height = 5)
  plot(roc(dat %>% filter(Index == x),response = Type, predictor =expression,
           levels=c("Control", "Diesease"),direction = "<"),
       lwd=5,
       lty=1,
       print.auc=TRUE,
       print.auc.x=0.4, print.auc.y=0.5,
       auc.polygon=TRUE,
       auc.polygon.col="#fff7f7",
       grid=c(0.5, 0.2),
       grid.col=c("black", "black"),
       print.thres=TRUE,
       main= paste(x,"_ROC curves",sep = ""),
       col="#FF2E63",
       legacy.axes=TRUE)
  dev.off()
})

###准备找验证集进行验证######## (Prepare to find and validate the dataset)

rm(list=ls()) #清空环境内变量 (Clear all variables in the environment)
options(stringsAsFactors = F) #避免自动将字符串转换为R语言因子 (Avoid automatic conversion of strings to R factors)
suppressMessages(library(GEOquery))
library(stringr)
library(ggplot2)
library(reshape2)
library(limma)
library("GEOquery")
setwd("E:/YQS-427-1WGCNA-Diff-ROC诊断/验证集合") (Set working directory)

#使用GEOquery()函数来下载GEO数据然后使用exprs()函数获取表达矩阵. (Use GEOquery() function to download GEO data and then use exprs() function to get the expression matrix)
gset = getGEO('GSE20168', destdir=".", AnnotGPL = F, getGPL = F) #设置后两个为T比较耗时且作用不大 (Setting the last two parameters to T is time-consuming and has little effect)
exp<-exprs(gset[[1]]) #exp即为表达矩阵 (exp is the expression matrix)

#进行重注释 (Perform re-annotation)

gpl = getGEO('GPL96',destdir = ".")#这种取出的对象是s4 (This object is of type s4)
id_pre = gpl@dataTable@table#使用@调用s4对象中某个的值 (Use @ to call a value in the s4 object)
colnames(id_pre)
head(id_pre)

#先都转换成GB_ACC (First, convert all to GB_ACC)
ids2 = id_pre[,c("ID","Gene Symbol")] #只需要"ID","GB_ACC"这两列 (Only need "ID", "GB_ACC" columns)
colnames(ids2) = c("probe_id","symbol") #给这两列重新取名 (Rename these two columns)

tail(sort(table(ids2$symbol)))#看每个基因用了些什么探针 (See which probes each gene used)
ids=as.data.frame(ids2)
table(rownames(exp) %in% ids$probe_id)#选取了exp中的行名和ids（s3）文件中的probe_id对应. table（）%in%（）的意思 (Selected the row names in exp and the probe_id in the ids(s3) file corresponding. table() %in%() means)
exp = exp[rownames(exp) %in% ids$probe_id,]#把没有对应名字的过滤 (Filter out those without corresponding names)
table(rownames(exp) %in% ids$probe_id)#对比探针，探针数只剩下26183的了 (Compare probes, only 26183 probes left)

ids=ids[match(rownames(exp),ids$probe_id),]#使用match函数把ids里的探针顺序和表达矩阵的统一. (Use match function to unify the probe order in ids with the expression matrix)
head(ids)
head(exp)#确定顺序对 (Confirm the order is correct)

tmp = by(exp,
ids$symbol,
function(x) rownames(x)[which.max(rowMeans(x))])
probes = as.character(tmp)
dim(exp)

exp = exp[rownames(exp) %in% probes,] # 过滤有多个探针的基因 (Filter out genes with multiple probes)
dim(exp)
head(tmp)
rownames(exp)=ids[match(rownames(exp),ids$probe_id),2]
head(exp)

#write.csv(x = exp,file = "exp.csv")
exp=read.csv(file = "exp.csv")
expr=exp

#计算行平均值，按降序排列 (Calculate row means and sort in descending order)
index=order(rowMeans(expr[,-1]),decreasing = T)

#调整表达谱的基因顺序 (Adjust the gene order of the expression profile)
expr_ordered=expr[index,]

#对于有重复的基因，保留第一次出现的那个，即行平均值大的那个 (For duplicate genes, keep the one that appears first, i.e., the one with the larger row mean)
keep=!duplicated(expr_ordered$X)

#得到最后处理之后的表达谱矩阵 (Obtain the final processed expression profile matrix)
expr_max=expr_ordered[keep,]
rownames(expr_max)=expr_max$X
expr_max=expr_max[,-1]
exp=expr_max
write.csv(x = exp,file = "exp.csv")

dim(exp)
group<- read.csv(file ="group.csv")
name=read.csv(file = "top5_hubgenes.csv")#基因HSPB1取别名HSPBAP1进行分析 (Analyze gene HSPB1 with alias HSPBAP1)

table(rownames(exp)%in%name$ID)
dim(name)

table(group$Type)
exp1=exp[name$ID,]

#数据处理 (Data processing)
bar_mat<-t(exp1)
table(group$Type)
group$Type=c(rep("Control",15),rep("Disease",14))
anno<-group
table(group$Type)

rownames(anno)=group$ID
anno$type2<-anno$Type
anno <- anno[rownames(bar_mat),]
bar_mat<-bar_mat[rownames(anno),]
bar_mat<-as.data.frame(bar_mat)
bar_mat$group=anno$Type

#绘制 (Drawing)
library(RColorBrewer)
library(ggpubr)
library(ggplot2)

因子水平 (Factor level)
group
bar_mat$group<-factor(bar_mat$group,levels=c("Control","Disease"))

颜色、分组比较设置 (Color and grouping comparison settings)
color <-c("#3399CC","#FF6600")
my_comparisons <- list(c("Control", "Disease"))

#循环 (Loop)

提取需要循环绘制的基因名 (Extract gene names for looping)
gc <- colnames(bar_mat)
dim(bar_mat)

#开始批量绘制 (Start batch drawing)



plist<-list()
for (i in 1:length(gc)){
  bar_tmp<-bar_mat[,c(gc[i],"group")]
  colnames(bar_tmp)<-c("Expression","group")
  pb1<-ggboxplot(bar_tmp,
                 x="group",
                 y="Expression",
                 color="group",
                 ylab = "Expression_of_gene(log2_transform)",
                 fill= NULL,
                 add = "jitter",
                 bxp.errorbar.width = 0.6,
                 width = 0.4,
                 size=0.01,
                 font.label = list(size=30),
                 palette = color)+
    theme(panel.background =element_blank())
  pb1<-pb1+theme(axis.line=element_line(colour="black"))+theme(axis.title.x = element_blank())
  pb1<-pb1+theme(axis.title.y = element_text(size = 15,vjust = 1,hjust = 1))+theme(axis.text.x = element_text(size = 15,angle = 45,vjust = 1,hjust = 1))
  pb1<-pb1+theme(axis.text.y = element_text(size = 15))+ggtitle(gc[i])+theme(plot.title = element_text(hjust = 0.5,size=15,face="bold"))
  #pb1<-pb1+theme(legend.position = "NA")#分组显示
  #pb1<-pb1+stat_compare_means(method="t.test",hide.ns = F,comparisons =my_comparisons,label="p.signif")#这个是加星星的。
  pb1<- pb1 + stat_compare_means(method = "t.test",comparisons =my_comparisons,aes(label = paste0(..method.., "\n", "p =", ..p.format..)))
  plist[[i]]<-pb1
}

#合并
library(cowplot)
pall<-plot_grid(plist[[1]],plist[[2]],plist[[3]],
                plist[[4]],plist[[5]],ncol=5)

pdf(file = "hubgenes_exp_human_box.pdf",width=18,height=5)
pall
dev.off()
#批量导出单图

dir()
gc
for (i in 1:5) {
  pdf(file = paste0(c(gc[i]),"_",".pdf"),width=10,height=5)
  p<- plist[[i]]
  print(p)
  dev.off()
}
#验证一下趋势
setwd("E:/YQS-427-1WGCNA-Diff-ROC诊断/test/脑的数据集")
diff=read.csv(file = "differ.csv")
rownames(diff)=diff$X
name
diff[name$Name,]
rownames(diff)
