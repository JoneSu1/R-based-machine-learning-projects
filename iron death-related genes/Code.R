library(tidyr)
library(dplyr)
library(limma)
library(ggpubr)
library(pheatmap)

#expressions file
GSE12288 <- read.table('GSE12288_series_matrix.txt',header = T,sep = '\t',check.names = F)
gpl96 <- read.table('F:/workspace/GPL/GPL96-57554.filt.txt',header = T,sep = '\t'.
                  check.names = F)
GSE12288$Symbol <- gpl96$`Gene Symbol`
GSE12288 <- GSE12288[,c(ncol(GSE12288),1:(ncol(GSE12288)-1))]
GSE12288 <- na.omit(GSE12288)
GSE12288 <- GSE12288[!GSE12288$Symbol=='',]
GSE12288 <- GSE12288[,-2]
newdata <- data.frame()
# Keep only one
for (i in 1:nrow(GSE12288)){
  #for (i in 1:2){
  if (grepl('///',GSE12288[i,1])){
    #print(i)
    a <- strsplit(GSE12288[i,1], '///', fixed=TRUE)[[1]]
    #print(a)
    newdata <-rbind(newdata,cbind(data.frame(Symbol=a[1]),GSE12288[i,2:ncol(GSE12288)]))
  } else {
    newdata <- rbind(newdata,GSE12288[i,])
  }
}
save(newdata,file='newdata.Rdata')
newdata$Symbol <- gsub(' ','',newdata$Symbol)
GSE12288.unique <- aggregate(. ~Symbol,newdata,mean)
rownames(GSE12288.unique ) <- GSE12288.unique$Symbol
sample <- read.table('GSE12288_sample.filt.txt',header = T,sep = '\t',check.names = F)

#Control 112,Case 110
GSE12288.final <- GSE12288.unique[,sample$Acc]
library(limma)
GSE12288.final=normalizeBetweenArrays(GSE12288.final)
sample <- read.table('GSE12288_sample.filt.txt',header = T,sep = '\t',check.names = F)
group_list <- factor(sample$Type,levels = c('control','case'))
boxplot(GSE12288.final,outline=FALSE, notch=T,col=group_list, las=2)
write.table(cbind(data.frame(Symbol=rownames(GSE12288.final),GSE12288.final)),file = 'GSE12288.final.txt',sep = '\t',row.names = F,quote = F)


#GSE20680#
library(GEOquery)
GSE20680 <- read.table('GSE20680_series_matrix.txt',header = T,sep = '\t',check.names = F)
gpl4133 <- Table(getGEO(filename='F:/workspace/GPL/GPL4133.annot.gz'))
gpl4133 <- gpl4133[,c(1,3)]
gpl4133 <- gpl4133[!gpl4133$`Gene symbol`=='',]
GSE20680 <- merge(GSE20680,gpl4133,by.x = 'ID_REF',by.y = 'ID')
GSE20680 <- GSE20680[,c(ncol(GSE20680),1:(ncol(GSE20680)-1))]
GSE20680 <- GSE20680[,-2]
colnames(GSE20680)[1] <- 'Symbol'
newdata2 <- data.frame()
# Keep only one
for (i in 1:nrow(GSE20680)){
  #for (i in 1:2){
  if (grepl('///',GSE20680[i,1])){
    #print(i)
    a <- strsplit(GSE20680[i,1], '///', fixed=TRUE)[[1]]
    #print(a)
    newdata2 <-rbind(newdata2,cbind(data.frame(Symbol=a[1]),GSE20680[i,2:ncol(GSE20680)]))
  } else {
    newdata2 <- rbind(newdata2,GSE20680[i,])
  }
}
save(newdata2,file='newdata2.Rdata')
newdata2$Symbol <- gsub(' ','',newdata2$Symbol)
newdata2[,2:ncol(newdata2)] <- apply(newdata2[,2:ncol(newdata2)],2,function(x){as.numeric(x)})
GSE20680.unique <- aggregate(. ~Symbol,newdata2,mean)
rownames(GSE20680.unique ) <- GSE20680.unique$Symbol
sample <- read.table('GSE20680_sample.filt.txt',header = T,sep = '\t',check.names = F)
#Control 52,Case 87
GSE20680.final <- GSE20680.unique[,sample$Acc]
library(limma)
GSE20680.final=normalizeBetweenArrays(GSE20680.final)
sample <- read.table('GSE20680_sample.filt.txt',header = T,sep = '\t',check.names = F)
group_list <- factor(sample$Type,levels = c('Control','Case'))
boxplot(GSE20680.final,outline=FALSE, notch=T,col=group_list, las=2)
write.table(cbind(data.frame(Symbol=rownames(GSE20680.final),GSE20680.final)),file = 'GSE20680.final.txt',sep = '\t',row.names = F,quote = F)


#GSE12288 Differential Expression ####
# sample <- read.table('GSE12288_sample.filt.txt',header = T,sep = '\t',check.names = F)
# group_list <- factor(sample$Type,levels = c('control','case'))
# pvalue < -0.05
# logFoldChange <- 0
# dat <- apply(GSE12288.final,2,function(x){log2(x+1)})
# design = model.matrix(~0+group_list, data=group_list)
# colnames(design) <- levels(group_list)
# rownames(design) <- colnames(dat)
# design
# # Construct the difference comparison matrix
# contrast.matrix <- makeContrasts(case-control, levels = design)
# fit <- lmFit(dat, design)
# fit2 <- contrasts.fit(fit, contrast.matrix)
# fit2 <- eBayes(fit2)
# options(digits = 4) # output the result of the decimal point to retain 4 digits
# allDiff=topTable(fit2,coef=1,number=Inf)
# allDiff <- na.omit(allDiff)
# write.table(cbind(Symbol=rownames(allDiff),allDiff),file="GSE12288.limmaOut.txt",sep="\t",row.names = F,quote = F)
# 
# diffSig = allDiff[(allDiff$P.Value < pvalue & (allDiff$logFC>logFoldChange | allDiff$logFC<(-logFoldChange))),]
# write.table(cbind(Symbol=rownames(diffSig),diffSig), file="GSE12288.diffSig.txt",sep="\t",row.names = F,quote = F)
# 
# diffUp = allDiff[(allDiff$P.Value < pvalue & (allDiff$logFC>logFoldChange)),]
# write.table(cbind(Symbol=rownames(diffUp),diffUp), file="GSE12288.up.txt",sep="\t",row.names = F,quote = F)
# 
# diffDown = allDiff[(allDiff$P.Value < pvalue & (allDiff$logFC<(-logFoldChange))),]
# write.table(cbind(Symbol=rownames(diffDown),diffDown), file="GSE12288.down.txt",sep="\t",row.names = F,quote = F)
# 
#GSE20680 Difference Expression ####
# sample <- read.table('GSE20680_sample.filt.txt',header = T,sep = '\t',check.names = F)
# group_list <- factor(sample$Type,levels = c('Control','Case'))
# pvalue < -0.05
# logFoldChange <- 0
# dat <- GSE20680.final
# design = model.matrix(~0+group_list, data=group_list)
# colnames(design) <- levels(group_list)
# rownames(design) <- colnames(dat)
# design
# # Construct the difference comparison matrix
# contrast.matrix <- makeContrasts(Case-Control, levels = design)
# fit <- lmFit(dat, design)
# fit2 <- contrasts.fit(fit, contrast.matrix)
# fit2 <- eBayes(fit2)
# options(digits = 4) # output the result of the decimal point to retain 4 digits
# allDiff=topTable(fit2,coef=1,number=Inf)
# allDiff <- na.omit(allDiff)
# write.table(cbind(Symbol=rownames(allDiff),allDiff),file="GSE20680.limmaOut.txt",sep="\t",row.names = F,quote = F)
# 
# diffSig = allDiff[(allDiff$P.Value < pvalue & (allDiff$logFC>logFoldChange | allDiff$logFC<(-logFoldChange))),]
# write.table(cbind(Symbol=rownames(diffSig),diffSig), file="GSE20680.diffSig.txt",sep="\t",row.names = F,quote = F)
# 
# diffUp = allDiff[(allDiff$P.Value < pvalue & (allDiff$logFC>logFoldChange)),]
# write.table(cbind(Symbol=rownames(diffUp),diffUp), file="GSE20680.up.txt",sep="\t",row.names = F,quote = F)
# 
# diffDown = allDiff[(allDiff$P.Value < pvalue & (allDiff$logFC<(-logFoldChange))),]
# write.table(cbind(Symbol=rownames(diffDown),diffDown), file="GSE20680.down.txt",sep="\t",row.names = F,quote = F)


#VolcanoMap
# pvalue < -0.05
# logFC <- 1
# allDiff <- read.table('limmaOut.txt',header = T,check.names = F,row.names = 1,sep = '\t')
# allDiff$Significant <- ifelse(allDiff$P.Value<pvalue & abs(allDiff$logFC)>= logFC.
# ifelse(allDiff$logFC> logFC,'up','down'),'no')
# mycol <- c("#3CB371", "#3D3D3D", "#FF4500")
# pdf(file="Fig1.MDD.DEGs.pdf",width=6,height=5.5)
# p <- ggplot(allDiff, aes(logFC, -log10(P.Value), colour= Significant))+
# geom_point(size=1.2)+theme_bw()+
# scale_color_manual(values = mycol,name='Significant')+
# labs(title="DEGs",x="Log2FC",y="-log10 (P.value)")+
# geom_hline(yintercept = -log10(pvalue),linetype=3,lwd = 1)+
# geom_vline(xintercept = c(-logFC, logFC), linetype=3,lwd = 1)+
# theme(axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15).
# legend.text = element_text(size = 15),text = element_text(size = 15).
# plot.title = element_text(hjust = 0.5))
# print(p)
# dev.off()

#t test for expression of iron death genes
fe <- read.csv('F:/workspace/Fe_gene.csv',header = T) %>% . $Symbol
fe <- unique(fe)
index1 <- rownames(GSE12288.final)[rownames(GSE12288.final) %in% fe]
index2 <- rownames(GSE20680.final)[rownames(GSE20680.final) %in% fe]
index2 <- index2[-grep('HBA1',index2)]

#GSE12288
hubexp <- GSE12288.final[index1,]
hubexp <- na.omit(hubexp)
hubexp <- data.frame(t(hubexp),check.names = F)
hubexp$Type <- factor(c(rep('Control',112),rep('Case',110)),levels = c('Control','Case'))
hubexp[,1:(ncol(hubexp)-1)] <- apply(hubexp[,1:(ncol(hubexp)-1)],2,function(x){log2(x+1)})
#t test
tdata <- data.frame()
for (i in colnames(hubexp)[1:(ncol(hubexp)-1)]{
  data = hubexp[,c(i,'Type')]
  colnames(data) <- c('Gene','Type')
  test <- t.test(Gene~Type,data)
  tdata <- rbind(tdata,data.frame(Gene = i,Pvalue = test$p.value))
}
hub <- tdata[tdata$Pvalue<0.05,]$Gene
hubexp <- hubexp[,c(hub,'Type')]
dbox <- gather(hubexp,FeGene,Expression,1:(ncol(hubexp)-1))
pdf('Fig1.GSE12288.Fe.pdf',width = 7.5,height = 4)
p <- ggboxplot(dbox,x='FeGene',y='Expression',fill = 'Type'.
               ylab = 'Relative Expression'.
               xlab = 'FeGene',palette = "jco".
               #add = 'jitter'.
               size =0.6)+ylim(c(0,11))+
  theme(axis.text.x=element_text(size=10,angle = 60,hjust=1).
        axis.text.y=element_text(size=13).
        axis.title = element_text(size=15).
        legend.text = element_text(size=15).
        legend.title = element_text(size=15).
        axis.line = element_line(size=1.2))+
  stat_compare_means(size = 4.5,aes(group=Type).
                     label = 'p.signif',method = 't.test',label.y = 10.5)
p + font('xlab',face = 'bold') + font('ylab',face = 'bold') +
  font('x.text',face = 'bold')+font('y.text',face = 'bold')+
  font('legend.title',face = 'bold')+font('legend.text',face = 'bold')
dev.off()


#GSE20680
# hubexp <- GSE20680.final[index2,]
# hubexp <- na.omit(hubexp)
# hubexp <- data.frame(t(hubexp),check.names = F)
# hubexp$Type <- c(rep('Control',52),rep('Case',87))
#t test
# tdata <- data.frame()
# for (i in colnames(hubexp)[1:(ncol(hubexp)-1)]){
# data = hubexp[,c(i,'Type')]
# colnames(data) <- c('Gene','Type')
# test <- t.test(Gene~Type,data)
# tdata <- rbind(tdata,data.frame(Gene = i,Pvalue = test$p.value))
# }
# hub <- tdata[tdata$Pvalue<0.05,]$Gene
# hubexp <- hubexp[,c(hub,'Type')]
# dbox <- gather(hubexp,FeGene,Expression,1:(ncol(hubexp)-1))
#pdf('GSE12288.fe.pdf',width = 7.5,height = 5)
# p <- ggboxplot(dbox,x='FeGene',y='Expression',fill = 'Type'.
# ylab = 'Relative Expression'.
# xlab = 'FeGene',palette = "jco".
#add ='jitter'.
# size =0.8)+
# theme(axis.text.x=element_text(size=10,angle = 60,hjust=1).
# axis.text.y=element_text(size=13).
# axis.title = element_text(size=15).
# legend.text = element_text(size=15).
# legend.title = element_text(size=15).
# axis.line = element_line(size=1.2))+
# stat_compare_means(size = 4.5,aes(group=Type).
# label = 'p.signif',method = 't.test')
# p + font('xlab',face = 'bold') + font('ylab',face = 'bold') +
# font('x.text',face = 'bold')+font('y.text',face = 'bold')+
# font('legend.title',face = 'bold') + font('legend.text',face = 'bold')
# dev.off()



#Heatmap
diff <- hub
diffexp <- GSE12288.final[diff,]
sample <- read.table('GSE12288_sample.filt.txt',header = T,sep = '\t',check.names = F)
annotation_col <- data.frame(Type = factor(sample$Type,levels = c("Control", "Case")))
rownames(annotation_col) <- colnames(diffexp)
color.key <- c("#3300CC", "#3399FF", "white", "#FF3333", "#CC0000")
pdf(file = 'Fig2.DEGs.MDD.heatmap.pdf',width = 7,height = 4.5) ## When actually making the map, screen the differential genes first, then draw the map with the differential genes; or draw the map with a large variation of Top several thousand genes
pheatmap(diffexp,cellwidth = 1.2
         ,cellheight = 14.
         method="pearson", # method to calculate correlation between gene or sample, optionally "pearson" (default), "kendall", or "spearman"
         scale="row", #scale for genes
         cluster_rows=T,#clustering for genes ## don't cluster the samples, otherwise the sample name order will change
         cluster_cols=F,#No clustering for sample
         color = colorRampPalette(color.key)(100),## the color field that characterizes the high and low expressions
         show_colnames=F,show_rownames =T.
         annotation_col = annotation_col.
         #treeheight_row = "0",treeheight_col = "0",#not draw tree
         border_color = "NA")
dev.off()

#correlation analysis ####
library(corrplot)
library(psych)
#Matrix as a whole calculates the correlation itself
diff <- hub
diffexp <- GSE12288.final[diff,]
rt <- data.frame(t(diffexp),check.names = F)
Cor <- corr.test(rt,method = "pearson")
M <- Cor$r
p.mat <- Cor$p
head(M)
col2 <- colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#FFFFFF". 
                           "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"))

#Relevance Heat Map
pdf('Fig3.corrplot.pdf',width = 7,height = 7)
corrplot(M, type="full", order="hclust". 
         p.mat = p.mat, sig.level = 0.05.
         insig = "blank",col=col2(200))
dev.off()


#Enrichment Analysis
library(stringr)
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

overlap <- read.table('hub.txt',header = T) %>% . $Symbol
gene <- mapIds(org.Hs.eg.db,keys = overlap,column = 'ENTREZID',keytype = 'SYMBOL',multiVals='filter')
save(gene,file='gene.Rdata')
kk <- enrichGO(gene = gene.
               OrgDb = org.Hs.eg.db. 
               pvalueCutoff =0.05. 
               qvalueCutoff = 0.2.
               ont="all".
               readable =T
)
k <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.2)  

write.table(kk,file="enrichGO.txt",quote=F,row.names = F,sep = '\t')     
write.table(k,file="enrichKEGG.txt",quote=F,row.names = F,sep = '\t') 

#Mapping
pdf(file="Fig5.enrichGO.pdf",width = 9.5,height = 7)
barplot(kk, drop = TRUE,showCategory = 15,split="ONTOLOGY") + 
  scale_y_discrete(labels=function(x) str_wrap(x, width=80))+
  facet_grid(ONTOLOGY~. , scale='free') 
dev.off()


#KEGG draws Bagua diagram
boxdata <- read.table('enrichKEGG.txt',header = T,sep='\t',check.names = F)
#boxdata <- boxdata[order(boxdata$LogP,decreasing = T),]
pdf('Fig6.enrichKEGG.pdf',width = 10,height = 4)
ggplot(boxdata,aes(x=Count,y=Description)) +
  geom_bar(stat='identity',aes(fill=-log10(p.adjust)))+
  scale_fill_gradient(low="Khaki1",high="Firebrick1")+
  theme_bw()+xlab('Count')+ylab('')+
  theme(axis.title = element_text(size=18).
        axis.text.x = element_text(size=14).
        axis.text.y = element_text(size=17).
        legend.text = element_text(size=14).
        legend.title = element_text(size=14))+
  theme(panel.grid.major = element_blank().
        panel.grid.minor = element_blank(),axis.line = element_line(color = "black"))
dev.off()



#LASSO
library(glmnet)
library(pROC)
tar <- read.table('hub.txt',header = T) %>% . $Symbol
rtt <- data.frame(t(GSE12288.final[tar,]),check.names = F)
rtt$Type <- c(rep('Control',112),rep('Case',110))
roc2 <- data.frame()
x=as.matrix(rtt[,c(1:(ncol(rtt)-1))])
y=data.matrix(rtt$Type)
roc2 <- data.frame()
for (i in 1:1000) {
  set.seed(402)
  fit <- glmnet(x, y, family = "binomial", type.measure="class")
  pdf(file = "Fig7.lasso.pdf", width = 9, height = 5)
  par(mfrow=c(1,2))
  plot(fit, xvar = "lambda", label = TRUE)
  #dev.off()
  cvfit <- cv.glmnet(x, y, family="binomial",type.measure = 'class')
  plot(cvfit)
  abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
  dev.off()
  coef <- coef(fit, s = cvfit$lambda.min)
  index <- which(as.matrix(coef)! = 0)
  actCoef <- coef[index]
  lassoGene=row.names(coef)[index]
  lassoGene
  write.table(lassoGene,file="lassoGene.txt",sep="\t",quote=F,row.names=F,col.names=F)
  out <- predict(cvfit,newx=x,s='lambda.min',type='class')
  group <- data.frame(Sample=rownames(rtt),Type=rtt$Type)
  rownames(group) <- group$Sample
  out <- data.frame(out,check.names = F)
  out$Type <- group[rownames(out),]$Type
  colnames(out) <- c('predict','Type')
  #out$predict <- sapply(out$predict,function(x){as.numeric(x)})
  out$predict <- factor(out$predict,labels = c(0,1))
  out$Type <- factor(out$Type,labels = c(0,1))
  pdf('Fig8.lasso.ROC.pdf',width = 6,height = 6)
  p <- plot.roc(as.numeric(out[,2]),as.numeric(out[,1]),ylim=c(0,1),xlim=c(1,0).
                smooth=F, #draw smooth curve
                main="", auc.polygon=T.
                auc.polygon.col="grey91".
                #print.thres="best", #write the threshold on the graph whose sum of sensitivity+ specificity is the largest
                col='MidnightBlue',# the color of the line
                lwd=2, #Thickness of the line
                legacy.axes=T,print.auc=T)
  # roc2 <- rbind(roc2,data.frame(LG = paste(lassoGene[-1],sep = '',collapse = ',').
  # AUC= as.numeric(p$auc),seed=i,n = length(lassoGene[-1]) ))
  dev.off()
}


#SVM####
library(tidyverse)
library(glmnet)
source('msvmRFE.R') #Folder comes with
library(VennDiagram)
library(sigFeature)
library(e1071)
library(caret)
library(randomForest)

load('rtt.Rdata')
input <- rtt[,c(ncol(rtt),1:(ncol(rtt)-1))]
The first column of data in #input is the sample type, 0 and 1
input$Type <- ifelse(input$Type=='Control',0,1)
svminfo <- data.frame()
for ( i in 1:10) {
  #set.seed(6666)
  set.seed(i)
  svmRFE(input, k = 10, halve.above = 20) #halve.above splits the data and assigns random numbers
  nfold = 10
  nrows = nrow(input)
  folds = rep(1:nfold, len=nrows)[sample(nrows)
  folds = lapply(1:nfold, function(x) which(folds == x))
  # results = lapply(folds, svmRFE.wrap, input, k=10, halve.above=50) # feature selection
  results = lapply(folds, svmRFE.wrap, input, k=10, halve.above=20) # feature selection
  top.features = WriteFeatures(results, input, save=F) # View main variables
  head(top.features)
  #Save the features found by SVM-REF to a file
  # write.csv(top.features, "feature_svm.csv",row.names = F)
  # Enable parallel as the backend for foreach parallel computation
  library(doParallel)
  detectCores()
  cl <- makeCluster(3)
  registerDoParallel(cl)
  # If you do not fill in the field, the default is to return a list.
  # featsweep <- foreach(x=1:180,.combine='rbind', .packages = c("e1071")) %dopar% FeatSweep.wrap(x,results, input)
  #180 is the number of genes
  featsweep <- foreach(x=1:18, .packages = c("e1071")) %dopar% FeatSweep.wrap(x,results, input)
  # Don't forget to end the parallel
  stopCluster(cl)
  # save(featsweep,file = "featsweep.RData")
  
  # Runtime depends mainly on the number of variables selected, it is better not to select too many variables for a normal computer
  # Pick the first 5 variables for SVM model building and experience
  
  
  # Select the first 300 variables for SVM model construction, and then import the already run results
  #featsweep = lapply(1:300, FeatSweep.wrap, results, input) #300 variables
  # save(featsweep,file = "featsweep1030.RData")
  # (load("featsweep1030.RData"))
  # (load("featsweep1530.RData"))
  
  # Drawing
  no.info = min(prop.table(table(input[,1])))
  errors = sapply( featsweep, function(x) ifelse(is.null(x), NA, x$error))
  
  #dev.new(width=4, height=4, bg='white')
  # pdf("Fig13.B_svm-error.pdf", width = 5, height = 5)
  PlotErrors(errors, no.info=no.info) #View error rate
  dev.off()
  
  # dev.new(width=4, height=4, bg='white')
  # pdf("Fig14.B_svm-accuracy.pdf",width = 5,height = 5)
  Plotaccuracy(1-errors,no.info=no.info) #View accuracy
  dev.off()
  
  # The position where the red circle is located in the graph, i.e. the lowest error rate
  a <- which.min(errors) 
  
  # Compare the feature variables identified by lasso and SVM-REF method one and draw Venn diagram
  myoverlap <- intersect(lassoGene[-1], top.features[1:which.min(errors), "FeatureName"]) # intersect
  b <- paste(myoverlap,sep = '',collapse = ',')
  summary(lassoGene[-1]%in%top.features[1:which.min(errors), "FeatureName"]) 
  
  #ROC
  model <- svm(Type~. ,data = input, kernel = "linear")
  out <- predict(model,newx= input,type='class')
  group <- data.frame(Sample=rownames(input),Type=input$Type)
  rownames(group) <- group$Sample
  out <- data.frame(out,check.names = F)
  out$Type <- group[rownames(out),]$Type
  colnames(out) <- c('predict','Type')
  #out$predict <- sapply(out$predict,function(x){as.numeric(x)})
  # out$predict <- factor(out$predict,labels = c(0,1))
  # out$Type <- factor(out$Type,labels = c(0,1))
  # pdf('Fig15.SVM.ROC.pdf',width = 6,height = 6)
  p <- plot.roc(as.numeric(out[,2]),as.numeric(out[,1]),ylim=c(0,1),xlim=c(1,0).
                smooth=F, #draw smooth curve
                main="ROC of SVM", auc.polygon=T.
                auc.polygon.col="grey91".
                #print.thres="best", #write the threshold on the graph whose sum of sensitivity+ specificity is the largest
                col='MidnightBlue',# the color of the line
                lwd=2, #Thickness of the line
                legacy.axes=T,print.auc=T)
  roc2 <- as.numeric(p$auc)
  dev.off()
  svminfo <- rbind(svminfo,data.frame(seed=i,error = a,overlap = b,ROC=roc2))
}


#Single gene ROC, drawn on a graph
hub <- lassoGene[-1]
hubexp <- GSE20163.final[hub,]
hubexp <- na.omit(hubexp)
hubexp <- data.frame(t(hubexp),check.names = F)
hubexp$Type <- c(rep('Control',9),rep('PD',8))
rocdata <- data.frame(Sample = rownames(hubexp).
                      exp = hubexp$SLC18A2.
                      Type = hubexp$Type) 
#mycol <- brewer.pal(10,'Set3')
mycol <- brewer.pal(10,'Set2')
pdf('Fig19.markerROC.pdf', width = 7, height = 7)
x <- plot.roc(rocdata[,3],rocdata[,2],ylim=c(0,1),xlim=c(1,0).
              smooth=F, #draw smooth curve
              main="".
              #print.thres="best", #write the threshold on the graph whose sum of sensitivity+ specificity is the largest
              col=mycol[1],# color of the line
              lwd=2, #Thickness of the line
              legacy.axes=T)
j=1
auc.test <- paste0(colnames(hubexp)[1],' AUC : ',format(as.numeric(x$auc),digits=3))
for (i in colnames(hubexp[,2:(ncol(hubexp)-1)])){
  j=j+1
  rocdata <- data.frame(Sample = rownames(hubexp).
                        exp = hubexp[,i].
                        Type = hubexp$Type) 
  x <- plot.roc(rocdata[,3],rocdata[,2],ylim=c(0,1),xlim=c(1,0).
                smooth=F, #draw smooth curve
                main="".
                #print.thres="best", #write the threshold on the graph whose sum of sensitivity+ specificity is the largest
                col=mycol[j],# the color of the line
                lwd=2, #Thickness of the line
                legacy.axes=T,add=T)
  auc.test <- c(auc.test.
                paste0(i,' AUC : ',format(as.numeric(x$auc),digits=3)))
}
legend(0.6,0.2, auc.test,lwd=2,bty="n",col=mycol,cex = 1.3)
dev.off()



Correlation of #gene expression with coronary heart disease index
# All samples
# library(gmodels)
# sample <- read.table('GSE12288_sample.filt.txt',header = T,sep = '\t',check.names = F)
# cadi <- sample$CADi
# cadexp <- GSE12288.final[hub,]
# cadexp <- apply(cadexp,2,function(x){log2(x+1)})
# box2 <- lapply(hub,function(g){
# rt <- data.frame(Relative_Expression=as.numeric(cadexp[g,]).
# CADi=as.numeric(cadi))
# ggplot(data=rt, aes(x=Relative_Expression, y=CADi))+
# geom_point(color="blue")+ggtitle(g)+
# stat_smooth(method="lm",se=FALSE,color='red')+
# stat_cor(data=rt, method = "pearson")+theme_bw()
# 
# }) 
# pdf('lm.pdf', width = 10, height = 12)
# plot_grid(plotlist = box2,nrow = 5)
# dev.off() 

#disease sample
library(gmodels)
sample <- read.table('GSE12288_sample.filt.txt',header = T,sep = '\t',check.names = F)
cadi <- sample$CADi[113:nrow(sample)]
cadexp <- GSE12288.final[hub,113:nrow(sample)]
# cadexp <- apply(cadexp,2,function(x){log2(x+1)})
box2 <- lapply(hub,function(g){
  rt <- data.frame(Expression=as.numeric(cadexp[g,]).
                   CADi=as.numeric(cadi))
  ggplot(data=rt, aes(x=Expression, y=CADi))+
    geom_point(color="blue")+ggtitle(g)+
    stat_smooth(method="lm",se=FALSE,color='red')+
    stat_cor(data=rt, method = "pearson")+theme_bw()
  
}) 
pdf('lm1.pdf', width = 10, height = 12)
plot_grid(plotlist = box2,nrow = 5)
dev.off()

########################CHIESE NOTICE##########
library(tidyr)
library(dplyr)
library(limma)
library(ggpubr)
library(pheatmap)

#表达量文件
GSE12288 <- read.table('GSE12288_series_matrix.txt',header = T,sep = '\t',check.names = F)
gpl96 <- read.table('F:/workspace/GPL/GPL96-57554.filt.txt',header = T,sep = '\t',
                  check.names = F)
GSE12288$Symbol <- gpl96$`Gene Symbol`
GSE12288 <- GSE12288[,c(ncol(GSE12288),1:(ncol(GSE12288)-1))]
GSE12288 <- na.omit(GSE12288)
GSE12288 <- GSE12288[!GSE12288$Symbol=='',]
GSE12288 <- GSE12288[,-2]
newdata <- data.frame()
#只保留一个
for (i in 1:nrow(GSE12288)){
  #for (i in 1:2){
  if (grepl('///',GSE12288[i,1])){
    #print(i)
    a <- strsplit(GSE12288[i,1], '///', fixed=TRUE)[[1]]
    #print(a)
    newdata <-rbind(newdata,cbind(data.frame(Symbol=a[1]),GSE12288[i,2:ncol(GSE12288)]))
  } else {
    newdata <- rbind(newdata,GSE12288[i,])
  }
}
save(newdata,file='newdata.Rdata')
newdata$Symbol <- gsub(' ','',newdata$Symbol)
GSE12288.unique <- aggregate(.~Symbol,newdata,mean)
rownames(GSE12288.unique ) <- GSE12288.unique$Symbol
sample <- read.table('GSE12288_sample.filt.txt',header = T,sep = '\t',check.names = F)

#Control 112,Case 110
GSE12288.final <- GSE12288.unique[,sample$Acc]
library(limma)
GSE12288.final=normalizeBetweenArrays(GSE12288.final)
sample <- read.table('GSE12288_sample.filt.txt',header = T,sep = '\t',check.names = F)
group_list <- factor(sample$Type,levels = c('control','case'))
boxplot(GSE12288.final,outline=FALSE, notch=T,col=group_list, las=2)
write.table(cbind(data.frame(Symbol=rownames(GSE12288.final),GSE12288.final)),file = 'GSE12288.final.txt',sep = '\t',row.names = F,quote = F)


#GSE20680#
library(GEOquery)
GSE20680 <- read.table('GSE20680_series_matrix.txt',header = T,sep = '\t',check.names = F)
gpl4133 <- Table(getGEO(filename='F:/workspace/GPL/GPL4133.annot.gz'))
gpl4133 <- gpl4133[,c(1,3)]
gpl4133 <- gpl4133[!gpl4133$`Gene symbol`=='',]
GSE20680 <- merge(GSE20680,gpl4133,by.x = 'ID_REF',by.y = 'ID')
GSE20680 <- GSE20680[,c(ncol(GSE20680),1:(ncol(GSE20680)-1))]
GSE20680 <- GSE20680[,-2]
colnames(GSE20680)[1] <- 'Symbol'
newdata2 <- data.frame()
#只保留一个
for (i in 1:nrow(GSE20680)){
  #for (i in 1:2){
  if (grepl('///',GSE20680[i,1])){
    #print(i)
    a <- strsplit(GSE20680[i,1], '///', fixed=TRUE)[[1]]
    #print(a)
    newdata2 <-rbind(newdata2,cbind(data.frame(Symbol=a[1]),GSE20680[i,2:ncol(GSE20680)]))
  } else {
    newdata2 <- rbind(newdata2,GSE20680[i,])
  }
}
save(newdata2,file='newdata2.Rdata')
newdata2$Symbol <- gsub(' ','',newdata2$Symbol)
newdata2[,2:ncol(newdata2)] <- apply(newdata2[,2:ncol(newdata2)],2,function(x){as.numeric(x)})
GSE20680.unique <- aggregate(.~Symbol,newdata2,mean)
rownames(GSE20680.unique ) <- GSE20680.unique$Symbol
sample <- read.table('GSE20680_sample.filt.txt',header = T,sep = '\t',check.names = F)
#Control 52,Case 87
GSE20680.final <- GSE20680.unique[,sample$Acc]
library(limma)
GSE20680.final=normalizeBetweenArrays(GSE20680.final)
sample <- read.table('GSE20680_sample.filt.txt',header = T,sep = '\t',check.names = F)
group_list <- factor(sample$Type,levels = c('Control','Case'))
boxplot(GSE20680.final,outline=FALSE, notch=T,col=group_list, las=2)
write.table(cbind(data.frame(Symbol=rownames(GSE20680.final),GSE20680.final)),file = 'GSE20680.final.txt',sep = '\t',row.names = F,quote = F)


# #GSE12288差异表达####
# sample <- read.table('GSE12288_sample.filt.txt',header = T,sep = '\t',check.names = F)
# group_list <- factor(sample$Type,levels = c('control','case'))
# pvalue <-0.05
# logFoldChange <- 0
# dat <- apply(GSE12288.final,2,function(x){log2(x+1)})
# design = model.matrix(~0+group_list, data=group_list)
# colnames(design) <- levels(group_list)
# rownames(design) <- colnames(dat)
# design
# # 构建差异比较矩阵
# contrast.matrix <- makeContrasts(case-control, levels = design)
# fit <- lmFit(dat, design)
# fit2 <- contrasts.fit(fit, contrast.matrix)
# fit2 <- eBayes(fit2)
# options(digits = 4)#输出结果的小数点后保留4位
# allDiff=topTable(fit2,coef=1,number=Inf)
# allDiff <- na.omit(allDiff)
# write.table(cbind(Symbol=rownames(allDiff),allDiff),file="GSE12288.limmaOut.txt",sep="\t",row.names = F,quote = F)
# 
# diffSig = allDiff[(allDiff$P.Value < pvalue & (allDiff$logFC>logFoldChange | allDiff$logFC<(-logFoldChange))),]
# write.table(cbind(Symbol=rownames(diffSig),diffSig), file="GSE12288.diffSig.txt",sep="\t",row.names = F,quote = F)
# 
# diffUp = allDiff[(allDiff$P.Value < pvalue & (allDiff$logFC>logFoldChange)),]
# write.table(cbind(Symbol=rownames(diffUp),diffUp), file="GSE12288.up.txt",sep="\t",row.names = F,quote = F)
# 
# diffDown = allDiff[(allDiff$P.Value < pvalue & (allDiff$logFC<(-logFoldChange))),]
# write.table(cbind(Symbol=rownames(diffDown),diffDown), file="GSE12288.down.txt",sep="\t",row.names = F,quote = F)
# 
# #GSE20680差异表达####
# sample <- read.table('GSE20680_sample.filt.txt',header = T,sep = '\t',check.names = F)
# group_list <- factor(sample$Type,levels = c('Control','Case'))
# pvalue <-0.05
# logFoldChange <- 0
# dat <- GSE20680.final
# design = model.matrix(~0+group_list, data=group_list)
# colnames(design) <- levels(group_list)
# rownames(design) <- colnames(dat)
# design
# # 构建差异比较矩阵
# contrast.matrix <- makeContrasts(Case-Control, levels = design)
# fit <- lmFit(dat, design)
# fit2 <- contrasts.fit(fit, contrast.matrix)
# fit2 <- eBayes(fit2)
# options(digits = 4)#输出结果的小数点后保留4位
# allDiff=topTable(fit2,coef=1,number=Inf)
# allDiff <- na.omit(allDiff)
# write.table(cbind(Symbol=rownames(allDiff),allDiff),file="GSE20680.limmaOut.txt",sep="\t",row.names = F,quote = F)
# 
# diffSig = allDiff[(allDiff$P.Value < pvalue & (allDiff$logFC>logFoldChange | allDiff$logFC<(-logFoldChange))),]
# write.table(cbind(Symbol=rownames(diffSig),diffSig), file="GSE20680.diffSig.txt",sep="\t",row.names = F,quote = F)
# 
# diffUp = allDiff[(allDiff$P.Value < pvalue & (allDiff$logFC>logFoldChange)),]
# write.table(cbind(Symbol=rownames(diffUp),diffUp), file="GSE20680.up.txt",sep="\t",row.names = F,quote = F)
# 
# diffDown = allDiff[(allDiff$P.Value < pvalue & (allDiff$logFC<(-logFoldChange))),]
# write.table(cbind(Symbol=rownames(diffDown),diffDown), file="GSE20680.down.txt",sep="\t",row.names = F,quote = F)


# #火山图
# pvalue <-0.05
# logFC <- 1
# allDiff <- read.table('limmaOut.txt',header = T,check.names = F,row.names = 1,sep = '\t')
# allDiff$Significant <- ifelse(allDiff$P.Value<pvalue & abs(allDiff$logFC)>= logFC,
#                               ifelse(allDiff$logFC> logFC,'up','down'),'no')
# mycol <- c("#3CB371","#3D3D3D","#FF4500")
# pdf(file="Fig1.MDD.DEGs.pdf",width=6,height=5.5)
# p <- ggplot(allDiff, aes(logFC, -log10(P.Value), colour= Significant))+
#   geom_point(size=1.2)+theme_bw()+
#   scale_color_manual(values = mycol,name='Significant')+
#   labs(title="DEGs",x="Log2FC",y="-log10 (P.value)")+
#   geom_hline(yintercept = -log10(pvalue),linetype=3,lwd = 1)+
#   geom_vline(xintercept = c(-logFC, logFC), linetype=3,lwd = 1)+
#   theme(axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15),
#         legend.text = element_text(size = 15),text = element_text(size = 15),
#         plot.title = element_text(hjust = 0.5))
# print(p)
# dev.off()

#t检验铁死亡基因的表达
fe <- read.csv('F:/workspace/Fe_gene.csv',header = T) %>% .$Symbol
fe <- unique(fe)
index1 <- rownames(GSE12288.final)[rownames(GSE12288.final) %in% fe]
index2 <- rownames(GSE20680.final)[rownames(GSE20680.final) %in% fe]
index2 <- index2[-grep('HBA1',index2)]

#GSE12288
hubexp <- GSE12288.final[index1,]
hubexp <- na.omit(hubexp)
hubexp <- data.frame(t(hubexp),check.names = F)
hubexp$Type <- factor(c(rep('Control',112),rep('Case',110)),levels = c('Control','Case'))
hubexp[,1:(ncol(hubexp)-1)] <- apply(hubexp[,1:(ncol(hubexp)-1)],2,function(x){log2(x+1)})
#t检验
tdata <- data.frame()
for (i in colnames(hubexp)[1:(ncol(hubexp)-1)]){
  data = hubexp[,c(i,'Type')]
  colnames(data) <- c('Gene','Type')
  test <- t.test(Gene~Type,data)
  tdata <- rbind(tdata,data.frame(Gene = i,Pvalue = test$p.value))
}
hub <- tdata[tdata$Pvalue<0.05,]$Gene
hubexp <- hubexp[,c(hub,'Type')]
dbox <- gather(hubexp,FeGene,Expression,1:(ncol(hubexp)-1))
pdf('Fig1.GSE12288.Fe.pdf',width = 7.5,height = 4)
p <- ggboxplot(dbox,x='FeGene',y='Expression',fill ='Type',
               ylab = 'Relative Expression',
               xlab ='FeGene',palette = "jco",
               #add ='jitter',
               size =0.6)+ylim(c(0,11))+
  theme(axis.text.x=element_text(size=10,angle = 60,hjust=1),
        axis.text.y=element_text(size=13),
        axis.title = element_text(size=15),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15),
        axis.line = element_line(size=1.2))+
  stat_compare_means(size = 4.5,aes(group=Type),
                     label = 'p.signif',method = 't.test',label.y = 10.5)
p + font('xlab',face = 'bold')+font('ylab',face = 'bold')+
  font('x.text',face = 'bold')+font('y.text',face = 'bold')+
  font('legend.title',face = 'bold')+font('legend.text',face = 'bold')
dev.off()


# #GSE20680
# hubexp <- GSE20680.final[index2,]
# hubexp <- na.omit(hubexp)
# hubexp <- data.frame(t(hubexp),check.names = F)
# hubexp$Type <- c(rep('Control',52),rep('Case',87))
# #t检验
# tdata <- data.frame()
# for (i in colnames(hubexp)[1:(ncol(hubexp)-1)]){
#   data = hubexp[,c(i,'Type')]
#   colnames(data) <- c('Gene','Type')
#   test <- t.test(Gene~Type,data)
#   tdata <- rbind(tdata,data.frame(Gene = i,Pvalue = test$p.value))
# }
# hub <- tdata[tdata$Pvalue<0.05,]$Gene
# hubexp <- hubexp[,c(hub,'Type')]
# dbox <- gather(hubexp,FeGene,Expression,1:(ncol(hubexp)-1))
# #pdf('GSE12288.fe.pdf',width = 7.5,height = 5)
# p <- ggboxplot(dbox,x='FeGene',y='Expression',fill ='Type',
#                ylab = 'Relative Expression',
#                xlab ='FeGene',palette = "jco",
#                #add ='jitter',
#                size =0.8)+
#   theme(axis.text.x=element_text(size=10,angle = 60,hjust=1),
#         axis.text.y=element_text(size=13),
#         axis.title = element_text(size=15),
#         legend.text = element_text(size=15),
#         legend.title = element_text(size=15),
#         axis.line = element_line(size=1.2))+
#   stat_compare_means(size = 4.5,aes(group=Type),
#                      label = 'p.signif',method = 't.test')
# p + font('xlab',face = 'bold')+font('ylab',face = 'bold')+
#   font('x.text',face = 'bold')+font('y.text',face = 'bold')+
#   font('legend.title',face = 'bold')+font('legend.text',face = 'bold')
# dev.off()



#热图
diff <- hub
diffexp <- GSE12288.final[diff,]
sample <- read.table('GSE12288_sample.filt.txt',header = T,sep = '\t',check.names = F)
annotation_col <- data.frame(Type = factor(sample$Type,levels = c("Control","Case")))
rownames(annotation_col) <- colnames(diffexp)
color.key <- c("#3300CC","#3399FF","white","#FF3333","#CC0000")
pdf(file = 'Fig2.DEGs.MDD.heatmap.pdf',width =7,height = 4.5)  ##实际作图时，先筛差异基因，再用差异基因画图；或者用变化大的Top几千个基因画图
pheatmap(diffexp,cellwidth = 1.2
         ,cellheight = 14,
         method="pearson", #计算gene或sample之间的相关性的方法，可选"pearson" (default), "kendall", or "spearman"
         scale="row", #为基因做scale
         cluster_rows=T,#为基因做聚类##不要给样本做聚类，否则样本名顺序会变
         cluster_cols=F,#不为sample做聚类
         color = colorRampPalette(color.key)(100),###表征表达量高低的颜色域
         show_colnames=F,show_rownames =T,
         annotation_col = annotation_col,
         #treeheight_row = "0",treeheight_col = "0",#不画树
         border_color = "NA")
dev.off()

#相关性分析####
library(corrplot)
library(psych)
#矩阵整体自身计算相关性
diff <- hub
diffexp <- GSE12288.final[diff,]
rt <- data.frame(t(diffexp),check.names = F)
Cor <- corr.test(rt,method = "pearson")
M <- Cor$r
p.mat <- Cor$p
head(M)
col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3", "#92C5DE","#D1E5F0", "#FFFFFF", 
                           "#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))

#相关性热图
pdf('Fig3.corrplot.pdf',width = 7,height = 7)
corrplot(M, type="full", order="hclust", 
         p.mat = p.mat, sig.level = 0.05,
         insig = "blank",col=col2(200))
dev.off()


#富集分析
library(stringr)
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

overlap <- read.table('hub.txt',header = T) %>% .$Symbol
gene <- mapIds(org.Hs.eg.db,keys = overlap,column = 'ENTREZID',keytype = 'SYMBOL',multiVals='filter')
save(gene,file='gene.Rdata')
kk <- enrichGO(gene = gene,
               OrgDb = org.Hs.eg.db, 
               pvalueCutoff =0.05, 
               qvalueCutoff = 0.2,
               ont="all",
               readable =T
)
k <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.2)  

write.table(kk,file="enrichGO.txt",quote=F,row.names = F,sep = '\t')               
write.table(k,file="enrichKEGG.txt",quote=F,row.names = F,sep = '\t') 

#作图
pdf(file="Fig5.enrichGO.pdf",width = 9.5,height = 7)
barplot(kk, drop = TRUE,showCategory = 15,split="ONTOLOGY") + 
  scale_y_discrete(labels=function(x) str_wrap(x, width=80))+
  facet_grid(ONTOLOGY~., scale='free') 
dev.off()


#KEGG画八卦图
boxdata <- read.table('enrichKEGG.txt',header = T,sep='\t',check.names = F)
#boxdata <- boxdata[order(boxdata$LogP,decreasing = T),]
pdf('Fig6.enrichKEGG.pdf',width = 10,height = 4)
ggplot(boxdata,aes(x=Count,y=Description)) +
  geom_bar(stat='identity',aes(fill=-log10(p.adjust)))+
  scale_fill_gradient(low="Khaki1",high="Firebrick1")+
  theme_bw()+xlab('Count')+ylab('')+
  theme(axis.title = element_text(size=18),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=17),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))
dev.off()



#LASSO
library(glmnet)
library(pROC)
tar <- read.table('hub.txt',header = T) %>% .$Symbol
rtt <- data.frame(t(GSE12288.final[tar,]),check.names = F)
rtt$Type <- c(rep('Control',112),rep('Case',110))
roc2 <- data.frame()
x=as.matrix(rtt[,c(1:(ncol(rtt)-1))])
y=data.matrix(rtt$Type)
roc2 <- data.frame()
for (i in 1:1000) {
  set.seed(402)
  fit <- glmnet(x, y, family = "binomial", type.measure="class")
  pdf(file = "Fig7.lasso.pdf",width = 9,height = 5)
  par(mfrow=c(1,2))
  plot(fit, xvar = "lambda", label = TRUE)
  #dev.off()
  cvfit <- cv.glmnet(x, y, family="binomial",type.measure = 'class')
  plot(cvfit)
  abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
  dev.off()
  coef <- coef(fit, s = cvfit$lambda.min)
  index <- which(as.matrix(coef)!= 0)
  actCoef <- coef[index]
  lassoGene=row.names(coef)[index]
  lassoGene
  write.table(lassoGene,file="lassoGene.txt",sep="\t",quote=F,row.names=F,col.names=F)
  out <- predict(cvfit,newx=x,s='lambda.min',type='class')
  group <- data.frame(Sample=rownames(rtt),Type=rtt$Type)
  rownames(group) <- group$Sample
  out <- data.frame(out,check.names = F)
  out$Type <- group[rownames(out),]$Type
  colnames(out) <- c('predict','Type')
  #out$predict <- sapply(out$predict,function(x){as.numeric(x)})
  out$predict <- factor(out$predict,labels = c(0,1))
  out$Type <- factor(out$Type,labels = c(0,1))
  pdf('Fig8.lasso.ROC.pdf',width = 6,height = 6)
  p <- plot.roc(as.numeric(out[,2]),as.numeric(out[,1]),ylim=c(0,1),xlim=c(1,0),
                smooth=F, #绘制平滑曲线
                main="", auc.polygon=T,
                auc.polygon.col="grey91",
                #print.thres="best", #把阈值写在图上，其sensitivity+ specificity之和最大
                col='MidnightBlue',#线的颜色
                lwd=2, #线的粗细
                legacy.axes=T,print.auc=T)
  # roc2 <- rbind(roc2,data.frame(LG = paste(lassoGene[-1],sep = '',collapse = ','),
  #                               AUC= as.numeric(p$auc),seed=i,n = length(lassoGene[-1]) ))
  dev.off()
}


#SVM####
library(tidyverse)
library(glmnet)
source('msvmRFE.R')   #文件夹内自带
library(VennDiagram)
library(sigFeature)
library(e1071)
library(caret)
library(randomForest)

load('rtt.Rdata')
input <- rtt[,c(ncol(rtt),1:(ncol(rtt)-1))]
#input第一列数据是样本类型，0和1
input$Type <- ifelse(input$Type=='Control',0,1)
svminfo <- data.frame()
for ( i in 1:10) {
  #set.seed(6666)
  set.seed(i)
  svmRFE(input, k = 10, halve.above = 20) #halve.above分割数据，分配随机数
  nfold = 10
  nrows = nrow(input)
  folds = rep(1:nfold, len=nrows)[sample(nrows)]
  folds = lapply(1:nfold, function(x) which(folds == x))
  # results = lapply(folds, svmRFE.wrap, input, k=10, halve.above=50) #特征选择
  results = lapply(folds, svmRFE.wrap, input, k=10, halve.above=20) #特征选择
  top.features = WriteFeatures(results, input, save=F) #查看主要变量
  head(top.features)
  #把SVM-REF找到的特征保存到文件
  # write.csv(top.features,"feature_svm.csv",row.names = F)
  # 启用parallel作为foreach并行计算的后端
  library(doParallel)
  detectCores()
  cl <- makeCluster(3)
  registerDoParallel(cl)
  # 并行计算方式 “.combine"表示合并方式，rbind就是按行成矩阵（data.frame),如果不填，默认返回一个list。
  # featsweep <- foreach(x=1:180,.combine='rbind', .packages = c("e1071")) %dopar% FeatSweep.wrap(x,results, input)
  #180是基因数
  featsweep <- foreach(x=1:18, .packages = c("e1071")) %dopar% FeatSweep.wrap(x,results, input)
  #别忘了结束并行
  stopCluster(cl)
  # save(featsweep,file = "featsweep.RData")
  
  # 运行时间主要取决于选择变量的个数，一般的电脑还是不要选择太多变量
  # 选前5个变量进行SVM模型构建，体验一下
  
  
  # 选前300个变量进行SVM模型构建，然后导入已经运行好的结果
  #featsweep = lapply(1:300, FeatSweep.wrap, results, input) #300个变量
  # save(featsweep,file = "featsweep1030.RData")
  # (load("featsweep1030.RData"))
  # (load("featsweep1530.RData"))
  
  # 画图
  no.info = min(prop.table(table(input[,1])))
  errors = sapply( featsweep, function(x) ifelse(is.null(x), NA, x$error))
  
  #dev.new(width=4, height=4, bg='white')
  # pdf("Fig13.B_svm-error.pdf",width = 5,height = 5)
  PlotErrors(errors, no.info=no.info) #查看错误率
  dev.off()
  
  # dev.new(width=4, height=4, bg='white')
  # pdf("Fig14.B_svm-accuracy.pdf",width = 5,height = 5)
  Plotaccuracy(1-errors,no.info=no.info) #查看准确率
  dev.off()
  
  # 图中红色圆圈所在的位置，即错误率最低点
  a <- which.min(errors) 
  
  # 比较lasso和SVM-REF方法一找出的特征变量，画Venn图
  myoverlap <- intersect(lassoGene[-1], top.features[1:which.min(errors), "FeatureName"]) #交集
  b <- paste(myoverlap,sep = '',collapse = ',')
  summary(lassoGene[-1]%in%top.features[1:which.min(errors), "FeatureName"]) 
  
  #ROC
  model <- svm(Type~.,data = input, kernel = "linear")
  out <- predict(model,newx= input,type='class')
  group <- data.frame(Sample=rownames(input),Type=input$Type)
  rownames(group) <- group$Sample
  out <- data.frame(out,check.names = F)
  out$Type <- group[rownames(out),]$Type
  colnames(out) <- c('predict','Type')
  #out$predict <- sapply(out$predict,function(x){as.numeric(x)})
  # out$predict <- factor(out$predict,labels = c(0,1))
  # out$Type <- factor(out$Type,labels = c(0,1))
  # pdf('Fig15.SVM.ROC.pdf',width = 6,height = 6)
  p <- plot.roc(as.numeric(out[,2]),as.numeric(out[,1]),ylim=c(0,1),xlim=c(1,0),
                smooth=F, #绘制平滑曲线
                main="ROC of SVM", auc.polygon=T,
                auc.polygon.col="grey91",
                #print.thres="best", #把阈值写在图上，其sensitivity+ specificity之和最大
                col='MidnightBlue',#线的颜色
                lwd=2, #线的粗细
                legacy.axes=T,print.auc=T)
  roc2 <- as.numeric(p$auc)
  dev.off()
  svminfo <- rbind(svminfo,data.frame(seed=i,error = a,overlap = b,ROC=roc2))
}


#单基因ROC，画在一幅图上
hub <- lassoGene[-1]
hubexp <- GSE20163.final[hub,]
hubexp <- na.omit(hubexp)
hubexp <- data.frame(t(hubexp),check.names = F)
hubexp$Type <- c(rep('Control',9),rep('PD',8))
rocdata <- data.frame(Sample = rownames(hubexp),
                      exp = hubexp$SLC18A2,
                      Type = hubexp$Type) 
#mycol <- brewer.pal(10,'Set3')
mycol <- brewer.pal(10,'Set2')
pdf('Fig19.markerROC.pdf',width = 7,height = 7)
x <- plot.roc(rocdata[,3],rocdata[,2],ylim=c(0,1),xlim=c(1,0),
              smooth=F, #绘制平滑曲线
              main="",
              #print.thres="best", #把阈值写在图上，其sensitivity+ specificity之和最大
              col=mycol[1],#线的颜色
              lwd=2, #线的粗细
              legacy.axes=T)
j=1
auc.test <- paste0(colnames(hubexp)[1],' AUC : ',format(as.numeric(x$auc),digits=3))
for (i in colnames(hubexp[,2:(ncol(hubexp)-1)])){
  j=j+1
  rocdata <- data.frame(Sample = rownames(hubexp),
                        exp = hubexp[,i],
                        Type = hubexp$Type) 
  x <- plot.roc(rocdata[,3],rocdata[,2],ylim=c(0,1),xlim=c(1,0),
                smooth=F, #绘制平滑曲线
                main="",
                #print.thres="best", #把阈值写在图上，其sensitivity+ specificity之和最大
                col=mycol[j],#线的颜色
                lwd=2, #线的粗细
                legacy.axes=T,add=T)
  auc.test <- c(auc.test,
                paste0(i,' AUC : ',format(as.numeric(x$auc),digits=3)))
}
legend(0.6,0.2, auc.test,lwd=2,bty="n",col=mycol,cex = 1.3)
dev.off()



#基因表达量与冠心病指数的相关性
# #全部样本
# library(gmodels)
# sample <- read.table('GSE12288_sample.filt.txt',header = T,sep = '\t',check.names = F)
# cadi <- sample$CADi
# cadexp <- GSE12288.final[hub,]
# cadexp <- apply(cadexp,2,function(x){log2(x+1)})
# box2 <- lapply(hub,function(g){
#   rt <- data.frame(Relative_Expression=as.numeric(cadexp[g,]),
#                    CADi=as.numeric(cadi))
#   ggplot(data=rt, aes(x=Relative_Expression, y=CADi))+
#     geom_point(color="blue")+ggtitle(g)+
#     stat_smooth(method="lm",se=FALSE,color='red')+
#     stat_cor(data=rt, method = "pearson")+theme_bw()
#   
# }) 
# pdf('lm.pdf',width = 10,height = 12)
# plot_grid(plotlist = box2,nrow = 5)
# dev.off() 

#疾病样本
library(gmodels)
sample <- read.table('GSE12288_sample.filt.txt',header = T,sep = '\t',check.names = F)
cadi <- sample$CADi[113:nrow(sample)]
cadexp <- GSE12288.final[hub,113:nrow(sample)]
# cadexp <- apply(cadexp,2,function(x){log2(x+1)})
box2 <- lapply(hub,function(g){
  rt <- data.frame(Expression=as.numeric(cadexp[g,]),
                   CADi=as.numeric(cadi))
  ggplot(data=rt, aes(x=Expression, y=CADi))+
    geom_point(color="blue")+ggtitle(g)+
    stat_smooth(method="lm",se=FALSE,color='red')+
    stat_cor(data=rt, method = "pearson")+theme_bw()
  
}) 
pdf('lm1.pdf',width = 10,height = 12)
plot_grid(plotlist = box2,nrow = 5)
dev.off() 
