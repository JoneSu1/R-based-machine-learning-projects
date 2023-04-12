library(tidyr)
library(dplyr)
library(survival)
library(survminer)
library(stringr)
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library(GSVA)
library(estimate)
library(ggplot2)
library(grid)
library(gridExtra)

# TCGA_BRCA.raw <- read.delim2('mRNAmatrix.txt',header = T,sep = '\t',check.names = F)
# save(TCGA_BRCA.raw,file='TCGA_BRCA.Rdata')
# ids <- read.csv('F:\\workspace\\生信资料\\mart_export.csv',header=T,check.names = F)
# ids.filt <- ids[,c(1,5,6)]
# ids.filt.sno <- ids.filt[ids.filt$`Gene type`=='snoRNA',]
# TCGA_BRCA.sno <- merge(ids.filt.sno,TCGA_BRCA.raw,by.y='id',by.x='Gene stable ID')
# TCGA_BRCA.sno <- TCGA_BRCA.sno[,-c(1,2)]
# TCGA_BRCA.sno[,2:ncol(TCGA_BRCA.sno)] <- apply(TCGA_BRCA.sno[,2:ncol(TCGA_BRCA.sno)],
#                                                2,function(x){as.numeric(x)})
# TCGA_BRCA.sno.unique <- aggregate(.~Symbol,TCGA_BRCA.sno,mean)
# write.table(TCGA_BRCA.sno.unique,'TCGA_BRCA.sno.unique.txt',
#             sep = '\t',row.names = F,quote = F)
# rownames(TCGA_BRCA.sno.unique) <- TCGA_BRCA.sno.unique$Symbol
# TCGA_BRCA.sno.unique <- TCGA_BRCA.sno.unique[,-1]
# TCGA_BRCA.sno.unique.log <- apply(TCGA_BRCA.sno.unique,2,function(x){log2(x+1)})

#SnoRNAic数据库直接处理好了癌症+正常的snoRNA表达矩阵####
#不能去重，去重后差异太小，这种特殊的snoRNA直接用原始ID
#读取snoRNA
#snoexp <- read.table('TCGA_BRCA.sno.cancer.txt',header = T,sep = '\t',check.names = F)
snoexp1 <- read.table('TCGA_BRCA.sno.cancer1.txt',header = T,sep = '\t',check.names = F)
#排序
sample <- read.table('sample_order.txt',header = T) %>% .$sample
#snoexp <- snoexp[,c('ENSEMBL','SYMBOL',sample)]
snoexp1 <- snoexp1[,c('snoRNA',sample)]
rownames(snoexp1) <- snoexp1$snoRNA
#TCGA_BRCA.sno.unique <- aggregate(.~Symbol,snoexp1,mean)
#write.table(TCGA_BRCA.sno.unique,'TCGA_BRCA.sno.unique.txt',sep = '\t',row.names = F,quote = F)
#dat <- TCGA_BRCA.sno.unique.log[,colSums(TCGA_BRCA.sno.unique.log)>1]
rownames(TCGA_BRCA.sno.unique) <- TCGA_BRCA.sno.unique$Symbol
TCGA_BRCA.sno.unique <- snoexp1[,-1]
TCGA_BRCA.sno.unique.log <- apply(TCGA_BRCA.sno.unique,2,function(x){log2(x+1)})
dat <- TCGA_BRCA.sno.unique.log
group_list=c(rep('Normal',104),rep('Tumor',1077))
## 强制限定顺序
group_list <- factor(group_list,levels = c("Normal","Tumor"))
pvalue <-0.05
logFoldChange <- 1
design = model.matrix(~0+group_list, data=group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(dat)
design
# 构建差异比较矩阵
contrast.matrix <- makeContrasts(Tumor-Normal, levels = design)
fit <- lmFit(dat, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
options(digits = 4)#输出结果的小数点后保留4位
allDiff=topTable(fit2,coef=1,number=Inf)
allDiff <- na.omit(allDiff)
write.table(cbind(Symbol=rownames(allDiff),allDiff),file="limmaOut.txt",sep="\t",row.names = F,quote = F)

diffSig = allDiff[(allDiff$P.Value < pvalue & (allDiff$logFC>logFoldChange | allDiff$logFC<(-logFoldChange))),]
#diffSig = allDiff[(allDiff$adj.P.Val < pvalue & (allDiff$logFC>logFoldChange | allDiff$logFC<(-logFoldChange))),]
write.table(cbind(Symbol=rownames(diffSig),diffSig), file="diffSig.txt",sep="\t",row.names = F,quote = F)

diffUp = allDiff[(allDiff$P.Value < pvalue & (allDiff$logFC>logFoldChange)),]
write.table(cbind(Symbol=rownames(diffUp),diffUp), file="up.txt",sep="\t",row.names = F,quote = F)

diffDown = allDiff[(allDiff$P.Value < pvalue & (allDiff$logFC<(-logFoldChange))),]
write.table(cbind(Symbol=rownames(diffDown),diffDown), file="down.txt",sep="\t",row.names = F,quote = F)

#火山图
pvalue <-0.05
logFC <- 1
allDiff <- read.table('limmaOut.txt',header = T,check.names = F,row.names = 1,sep = '\t')
allDiff$Significant <- ifelse(allDiff$P.Value<pvalue & abs(allDiff$logFC)>= logFC,
                              ifelse(allDiff$logFC> logFC,'up','down'),'no')
mycol <- c("#3CB371","#3D3D3D","#FF4500")
pdf(file="Fig1A.DESNOs.pdf",width=6,height=5.5)
p <- ggplot(allDiff, aes(logFC, -log10(P.Value), colour= Significant))+
  geom_point(size=1.2)+theme_bw()+
  scale_color_manual(values = mycol,name='Significant')+
  labs(title="DEGs",x="Log2FC",y="-log10 (P.value)")+
  geom_hline(yintercept = -log10(pvalue),linetype=3,lwd = 1)+
  geom_vline(xintercept = c(-logFC, logFC), linetype=3,lwd = 1)+
  theme(axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15),
        legend.text = element_text(size = 15),text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5))
print(p)
dev.off()

#热图
name  <- read.table('diffSig.symbol.txt',header = T,sep = '\t',check.names = F) %>% .$Name
diff <- read.table('diffSig.symbol.txt',header = T,sep = '\t',check.names = F) %>% .$Symbol
diffexp <- TCGA_BRCA.sno.unique.log[diff,]
annotation_col <- data.frame(Type = factor(group_list,levels = c("Normal","Tumor")))
rownames(annotation_col) <- colnames(diffexp)
rownames(diffexp) <- name
color.key <- c("#3300CC","#3399FF","white","#FF3333","#CC0000")
pdf(file = 'Fig1B.DEGsheatmap.pdf',width = 6,height =4)  ##实际作图时，先筛差异基因，再用差异基因画图；或者用变化大的Top几千个基因画图
pheatmap(diffexp,cellwidth = 0.2,cellheight = 1.7,
         method="spearman", #计算gene或sample之间的相关性的方法，可选"pearson" (default), "kendall", or "spearman"
         scale="row", #为基因做scale
         cluster_rows=T,#为基因做聚类##不要给样本做聚类，否则样本名顺序会变
         cluster_cols=F,#不为sample做聚类
         color = colorRampPalette(color.key)(50),###表征表达量高低的颜色域
         show_colnames=F,show_rownames =F,
         annotation_col = annotation_col,
         #treeheight_row = "0",treeheight_col = "0",#不画树
         border_color = "NA")
dev.off()


#构建风险模型####
#肿瘤样本！！！
tarexp <- diffexp[,105:ncol(diffexp)]
tarexp <- data.frame(t(tarexp),check.names = F)
tarexp$Patient <- substr(rownames(tarexp),1,12)
clinical <- read.table('clinical.txt',header = T,sep = '\t',check.names = F)
OSinput <- merge(clinical,tarexp,by = 'Patient')
write.table(OSinput,'OSinput.txt',sep = '\t',quote = F,row.names = F)
model.info <- data.frame()
rtt <- OSinput[,c(1,4,9,2,3,5,6,7,8,10:138)]
rtt$Type <- 'Tumor'
#seed =933
for (seed in c(1:1000)){
  seed <- 933
  set.seed(seed)
  index <- caret::createDataPartition(rtt[,'Type'], p =0.7)
  train <- rtt[index$Resample1,]
  test<- rtt[-index$Resample1,]
  write.table(train,'train.coxinput.txt',row.names = F,col.names = T,quote=F,sep = '\t')
  write.table(test,'test.coxinput.txt',row.names = F,col.names = T,quote=F,sep = '\t')
  coxinput <- train
  coxinput <- coxinput[,-ncol(coxinput)]
  
  #KM+单因素####
  print('start unicox')
  KM <- data.frame()
  for(i in colnames(coxinput[,10:ncol(coxinput)])){
    survival <- coxinput[,c('Patient','fustate','futime',i)]
    survival$exp <- ifelse(survival[,4]>median(survival[,4]),
                           'High_Expression','Low_Expression')
    diff <- survdiff(Surv(futime,fustate) ~ exp, data = survival)
    pValue=format((1-pchisq(diff$chisq,df=1)),digits=5)
    KM <- rbind(KM,data.frame(Gene=i,pvalue=pValue))
  }
    write.table(KM,file = 'KM.txt',sep = '\t',quote = F,row.names = F)
  km <- KM[KM$pvalue<0.05,]$Gene
  outTab=data.frame()
  if (length(km)>0){
    for(i in km){
      cox <- coxph(Surv(futime, fustate) ~ coxinput[,i], data = coxinput)
      coxSummary = summary(cox)
      coxP=coxSummary$coefficients[,"Pr(>|z|)"]
      outTab=rbind(outTab,
                   cbind(id=i,
                         HR=coxSummary$conf.int[,"exp(coef)"],
                         HR.95L=coxSummary$conf.int[,"lower .95"],
                         HR.95H=coxSummary$conf.int[,"upper .95"],
                         pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
      )
    }
    write.table(outTab,file="uniCox.xls",sep="\t",row.names=F,quote=F)
    outTab[,c(2:ncol(outTab))] <- apply(outTab[,c(2:ncol(outTab))],2,function(x){as.numeric(x)})
    outTab <- outTab[order(outTab$pvalue),]
    tabletext <- outTab
    tabletext[,2:ncol(tabletext)] <- apply(tabletext[,2:ncol(tabletext)],2,
                                           function(x){format(as.numeric(x),digits=4)})
    tabletext <- rbind(c("Variable","HR","lower 95%CI","upper 95%CI","pvalue"),tabletext)
    pdf('Fig2A.unicox.pdf',width = 10,height =6)
    forestplot(labeltext=tabletext,
               mean=c(NA,outTab$HR),#HR值
               lower=c(NA,outTab$HR.95L), #95%置信区间下限
               upper=c(NA,outTab$HR.95H),#95%置信区间上限
               graph.pos=6,#图在表中的列位置
               graphwidth = unit(.25,"npc"),#图在表中的宽度比
               fn.ci_norm="fpDrawDiamondCI",#box类型选择钻石
               col=fpColors(box="#00A896", lines="#02C39A", zero = "black"),#box颜色
               boxsize=0.4,#box大小固定
               lwd.ci=1.6,
               ci.vertices.height = 0.2,ci.vertices=T,#不显示区间
               zero=1,#zero线横坐标
               lwd.zero=2,#zero线宽
               xticks = c(0.6,1,1.8),#横坐标刻度根据需要可随意设置
               lwd.xaxis=2,
               xlab="HR",
               hrzl_lines=list("1" = gpar(lwd=2, col="black"),#第一行顶部加黑实线
                               "2" = gpar(lwd=1, col="grey50", lty=2),#第二行顶部加灰色虚线
                               "8" = gpar(lwd=2, col="black")),#最后一行底部加黑线，""中数字为nrow(tabletext) + 1
               txt_gp=fpTxtGp(label=gpar(cex=1.2),#各种字体大小设置
                              ticks=gpar(cex=0.85),
                              xlab=gpar(cex=1.2),
                              title=gpar(cex=1.2)),
               lineheight = unit(1.5,"cm"),#固定行高
               colgap = unit(0.5,"cm"),
               mar=unit(rep(1.5, times = 4), "cm"),
               new_page = F
    )
    dev.off()
    
    #多因素####
    print('start multicox')
    #outTab$pvalue <- sapply(outTab$pvalue,function(x){as.numeric(x)})
    outTab[,c(2:ncol(outTab))] <- apply(outTab[,c(2:ncol(outTab))],2,function(x){as.numeric(x)})
    unicox <- outTab[outTab$pvalue <= 0.05,]$id
    if (length(unicox)>1) {
      multi <- coxinput[,c('Patient','futime','fustate',unicox)]
      multi1 <- multi[,-1]
      rownames(multi1) <- multi$Patient
      multiCox=coxph(Surv(futime, fustate) ~ ., data = multi1)
      multiCox=step(multiCox,direction = "both")#根据AIC值过滤，向前向后方法进一步优化
      multiCoxSum=summary(multiCox)
      outTab1=data.frame()
      outTab1=cbind(
        coef=multiCoxSum$coefficients[,"coef"],
        HR=multiCoxSum$conf.int[,"exp(coef)"],
        HR.95L=multiCoxSum$conf.int[,"lower .95"],
        HR.95H=multiCoxSum$conf.int[,"upper .95"],
        pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
      if (length(row.names(outTab1)) >0){
        outTab1=cbind(id=row.names(outTab1),outTab1)
        write.table(outTab1,file="multiCox.xls",sep="\t",row.names=F,quote=F)
        outTab1 <- data.frame(outTab1)
        outTab1[,c(2:ncol(outTab1))] <- apply(outTab1[,c(2:ncol(outTab1))],2,function(x){as.numeric(x)})
        modelgene <- outTab1$id
        outTab1 <- outTab1[order(outTab1$pvalue),]
        tabletext <- outTab1
        tabletext[,2:ncol(tabletext)] <- apply(tabletext[,2:ncol(tabletext)],2,
                                               function(x){format(as.numeric(x),digits=4)})
        tabletext <- rbind(c("Variable",'coef',"HR","lower 95%CI","upper 95%CI","pvalue"),tabletext)
        pdf('Fig2B.multicox.pdf',width = 10,height =4)
        forestplot(labeltext=tabletext,
                   mean=c(NA,outTab1$HR),#HR值
                   lower=c(NA,outTab1$HR.95L), #95%置信区间下限
                   upper=c(NA,outTab1$HR.95H),#95%置信区间上限
                   graph.pos=7,#图在表中的列位置
                   graphwidth = unit(.25,"npc"),#图在表中的宽度比
                   fn.ci_norm="fpDrawDiamondCI",#box类型选择钻石
                   col=fpColors(box="#00A896", lines="#02C39A", zero = "black"),#box颜色
                   boxsize=0.3,#box大小固定
                   lwd.ci=1.6,
                   ci.vertices.height = 0.2,ci.vertices=T,#不显示区间
                   zero=1,#zero线横坐标
                   lwd.zero=2,#zero线宽
                   xticks = c(0.5,2),#横坐标刻度根据需要可随意设置
                   lwd.xaxis=2,
                   xlab="HR",
                   hrzl_lines=list("1" = gpar(lwd=2, col="black"),#第一行顶部加黑实线
                                   "2" = gpar(lwd=1, col="grey50", lty=2),#第二行顶部加灰色虚线
                                   "4" = gpar(lwd=2, col="black")),#最后一行底部加黑线，""中数字为nrow(tabletext) + 1
                   txt_gp=fpTxtGp(label=gpar(cex=1.2),#各种字体大小设置
                                  ticks=gpar(cex=0.85),
                                  xlab=gpar(cex=1.2),
                                  title=gpar(cex=1.2)),
                   lineheight = unit(2.5,"cm"),#固定行高
                   colgap = unit(0.5,"cm"),
                   mar=unit(rep(1.5, times = 4), "cm"),
                   new_page = F
        )
        dev.off()
        
        #预测训练集风险评分####
        riskScore=predict(multiCox,type="risk")
        coxGene=rownames(multiCoxSum$coefficients)
        coxGene=gsub("`","",coxGene)
        outCol=c('Patient',"futime","fustate",coxGene)
        write.table(cbind(multi[,outCol],riskScore),
                    file="risk.txt",
                    sep="\t",
                    quote=F,
                    row.names=F)
        rt <- read.table('risk.txt',header = T,sep = '\t',check.names = F)
        rt$risk <- ifelse(rt$riskScore>=median(rt$riskScore),'High_risk','Low_risk')
        rt <- rt[order(rt$risk,decreasing = T),]
        write.table(rt, file="risk.txt",sep="\t",quote=F,row.names=F)
        # #画风险曲线####
        riskdata <- read.table('risk.txt',header = T,check.names = T,sep = '\t')
        risk=riskdata[order(riskdata$riskScore),]
        median <- median(riskdata$riskScore)
        riskClass=riskdata[,"risk"]
        lowLength=length(riskClass[riskClass=="Low_risk"])
        highLength=length(riskClass[riskClass=="High_risk"])
        line=risk[,"riskScore"]
        #line[line>10]=10

        #pdf('train.riskscore+state.pdf',width = 10,height = 7)
        par(mar=c(6,4,2,3),mfrow=c(2,1))
        plot(line,
             type="p",
             pch=20,cex=0.7,
             xlab="Patients (increasing risk socre)",
             ylab="Risk score",
             col=c(rep("blue",lowLength),
                   rep("red",highLength)),main='Training set')
        legend('topleft',legend = c('Low_risk','High_risk'),
               col=c('blue','red'),pch = c(20,20),bty='n')
        abline(h=median,v=lowLength,lty=2)

        #生存状态
        color=as.vector(riskdata$fustate)
        color[color==1]="red"
        color[color==0]="blue"
        plot(riskdata$futime,
             pch=19,cex=0.7,
             xlab="Patients (increasing risk socre)",
             ylab="Survival time (days)",
             col=color,main='Training set')
        abline(v=lowLength,lty=2)
        legend('topright',legend = c('Live','Dead'),
               col=c('blue','red'),pch = c(20,20),bty='n')
        dev.off()
        #热图####
        rt <- read.table('risk.txt',header = T,sep = '\t',check.names = F)
        rownames(rt) <- rt$id
        rt <- rt[,c(4:(3+length(modelgene)),ncol(rt))]
        data <- rt[,1:(ncol(rt)-1)] %>% t(.)
        colnames(data) <- rownames(rt)
        #color.key <- c("#3300CC","#3399FF","white","#FF3333","#CC0000")
        color.key <- c("blue",'DodgerBlue1','Snow1','Firebrick1',"red")
        annotation_col = data.frame(Type = factor(rt$risk,levels = c('Low_risk','High_risk')))  ##分组，legend名称：Type,这里可以调整显示顺序
        rownames(annotation_col) <- rownames(rt)
        # ann_colors <- list(Type=c(Low_risk='DodgerBlue',
        #                    High_risk='GreenYellow'))
        pdf('Fig1G.train.risk.heatmap.pdf',width =8,height = 3)
        pheatmap(data,cellwidth = 0.4, cellheight = 60, fontsize = 12,
                 clustering_method="complete", #计算gene或sample之间的相关性的方法，可选"pearson" (default), "kendall", or "spearman"
                 scale="row", #为基因做scale
                 cluster_rows=T,#为基因做聚类##不要给样本做聚类，否则样本名顺序会变
                 cluster_cols=F,#不为sample做聚类
                 color = colorRampPalette(color.key)(100),###表征表达量高低的颜色域
                 show_colnames=F,show_rownames =T,
                 annotation_col = annotation_col,
                 #annotation_colors = ann_colors,##
                 #treeheight_row = "0",treeheight_col = "0",#不画树
                 border_color = "NA",main = 'Training Set')
        dev.off()
        
        #预测验证集风险评分####
        rtttt <- read.table('test.coxinput.txt',header = T,check.names = F,sep = '\t')
        rownames(rtttt) <- rtttt$Patient
        rtttt1 <- rtttt[,-ncol(rtttt)]
        riskScore=predict(multiCox,type="risk",newdata=rtttt1)
        coxGene=rownames(multiCoxSum$coefficients)
        coxGene=gsub("`","",coxGene)
        outCol=c("Patient","futime","fustate",coxGene)
        write.table(cbind(rtttt1[,outCol],riskScore),
                    file="test.risk.txt",
                    sep="\t",
                    quote=F,
                    row.names=F)
        rt <- read.table('test.risk.txt',header = T,sep = '\t',check.names = F)
        rt$risk <- ifelse(rt$riskScore>=median(rt$riskScore),'High_risk','Low_risk')
        rt <- rt[order(rt$risk,decreasing = T),]
        write.table(rt, file="test.risk.txt",sep="\t",quote=F,row.names=F)
        
        # #画验证集风险曲线#####
        riskdata <- read.table('test.risk.txt',header = T,check.names = T,sep = '\t')
        risk=riskdata[order(riskdata$riskScore),]
        median <- median(riskdata$riskScore)
        riskClass=riskdata[,"risk"]
        lowLength=length(riskClass[riskClass=="Low_risk"])
        highLength=length(riskClass[riskClass=="High_risk"])
        line=risk[,"riskScore"]
        #line[line>10]=10

        #pdf('test.riskscore+state.pdf',width = 8,height = 7)
        set.seed(123)
        par(mar=c(6,4,2,3),mfrow=c(2,1))
        plot(line,
             type="p",
             pch=20,cex=0.7,
             xlab="Patients (increasing risk socre)",
             ylab="Risk score",
             col=c(rep("blue",lowLength),
                   rep("red",highLength)),main='Testing set')
        legend('topleft',legend = c('Low_risk','High_risk'),
               col=c('blue','red'),pch = c(20,20),bty='n')
        abline(h=median,v=lowLength,lty=2)

        #生存状态
        color=as.vector(riskdata$fustate)
        color[color==1]="red"
        color[color==0]="blue"
        plot(riskdata$futime,
             pch=19,cex=0.7,
             xlab="Patients (increasing risk socre)",
             ylab="Survival time (days)",
             col=color,main='Testing set')
        abline(v=lowLength,lty=2)
        legend('topright',legend = c('Live','Dead'),
               col=c('blue','red'),pch = c(20,20),bty='n')
        dev.off()
        #验证集热图#####
        rt <- read.table('test.risk.txt',header = T,sep = '\t',check.names = F)
        rownames(rt) <- rt$id
        rt <- rt[,c(4:(3+length(modelgene)),ncol(rt))]
        data <- rt[,1:(ncol(rt)-1)] %>% t(.)
        colnames(data) <- rownames(rt)
        #color.key <- c("#3300CC","#3399FF","white","#FF3333","#CC0000")
        color.key <- c("blue",'DodgerBlue1','Snow1','Firebrick1',"red")
        annotation_col = data.frame(Type = factor(rt$risk,levels = c('Low_risk','High_risk')))  ##分组，legend名称：Type,这里可以调整显示顺序
        rownames(annotation_col) <- rownames(rt)
        # ann_colors <- list(Type=c(Low_risk='DodgerBlue',
        #                    High_risk='GreenYellow'))
        pdf('Fig1H.test.risk.heatmap.pdf',width =8,height = 3)
        pheatmap(data,cellwidth = 0.8, cellheight =60, fontsize = 12,
                 clustering_method="complete", #计算gene或sample之间的相关性的方法，可选"pearson" (default), "kendall", or "spearman"
                 scale="row", #为基因做scale
                 cluster_rows=T,#为基因做聚类##不要给样本做聚类，否则样本名顺序会变
                 cluster_cols=F,#不为sample做聚类
                 color = colorRampPalette(color.key)(100),###表征表达量高低的颜色域
                 show_colnames=F,show_rownames =T,
                 annotation_col = annotation_col,
                 #annotation_colors = ann_colors,##
                 #treeheight_row = "0",treeheight_col = "0",#不画树
                 border_color = "NA",main = 'Testing Set')
        dev.off()
        
        #训练集survival分析####
        print('start risk survival')
        survival <- read.table('risk.txt',sep = '\t',check.names = F,header = T)
        fit <- survfit(Surv(futime, fustate) ~ risk, data = survival)
        p <- ggsurvplot(fit, data = survival ,
                        ggtheme = theme_bw(), #想要网格就运行这行
                        conf.int = F, #不画置信区间，想画置信区间就把F改成T
                        #conf.int.style = "step",#置信区间的类型，还可改为ribbon
                        risk.table = TRUE, # 绘制累计风险曲线
                        #tables.theme = theme_void(),
                        tables.height = 0.25,
                        censor = T, #不显示观察值所在的位置
                        palette = c("Red3","DodgerBlue4"), #线的颜色对应高、低(自定义调色板)
                        xlab = "Time (Days)",
                        legend.title = "",#基因名写在图例题目的位置
                        font.legend = 15,#图例的字体大小
                        #font.title = 12,font.x = 10,font.y = 10,#设置其他字体大小
                        legend=c(0.8,0.8),
                        legend.labs = c('High_risk','Low_risk'),pval = TRUE)
        #在图例上标出高低分界点的表达量，和组内sample数量
        #legend.labs=c(paste0(">",round(svsort[fit$n[2],i],2),"(",fit$n[1],")"),
        #paste0("<",round(svsort[fit$n[2],i],2),"(",fit$n[2],")")),
        #标出pvalue、HR、95% CI
        #太小的p value标为p < 0.001)
        p$plot <- p$plot+ggplot2::annotate('text',x=4800,y=0.6,label=paste0('Low_risk (n = ',fit$n[1],')'),size=5)+
          ggplot2::annotate('text',x=1800,y=0.45,label=paste0('High_risk (n = ',fit$n[2],')'),size=5)
        pdf(file = "Fig1I.train.survival.pdf",width = 5.5,height = 6)
        print(p)
        dev.off()
        
        
        #验证集survival分析####
        survival <- read.table('test.risk.txt',sep = '\t',check.names = F,header = T)
        fit <- survfit(Surv(futime, fustate) ~ risk, data = survival)
        p <- ggsurvplot(fit, data = survival ,
                        ggtheme = theme_bw(), #想要网格就运行这行
                        conf.int = F, #不画置信区间，想画置信区间就把F改成T
                        #conf.int.style = "step",#置信区间的类型，还可改为ribbon
                        risk.table = TRUE, # 绘制累计风险曲线
                        #tables.theme = theme_void(),
                        tables.height = 0.25,
                        censor = T, #不显示观察值所在的位置
                        palette = c("Red3","DodgerBlue4"), #线的颜色对应高、低(自定义调色板)
                        xlab = "Time (Days)",
                        legend.title = "",#基因名写在图例题目的位置
                        font.legend = 15,#图例的字体大小
                        #font.title = 12,font.x = 10,font.y = 10,#设置其他字体大小
                        legend=c(0.8,0.8),
                        legend.labs = c('High_risk','Low_risk'),pval = TRUE)
        #在图例上标出高低分界点的表达量，和组内sample数量
        #legend.labs=c(paste0(">",round(svsort[fit$n[2],i],2),"(",fit$n[1],")"),
        #paste0("<",round(svsort[fit$n[2],i],2),"(",fit$n[2],")")),
        #标出pvalue、HR、95% CI
        #太小的p value标为p < 0.001)
        p$plot <- p$plot+ggplot2::annotate('text',x=5800,y=0.6,label=paste0('Low_risk (n = ',fit$n[1],')'),size=5)+
          ggplot2::annotate('text',x=2400,y=0.45,label=paste0('High_risk (n = ',fit$n[2],')'),size=5)
        pdf(file = "Fig1J.test.survival.pdf",width = 5.5,height = 6)
        print(p)
        dev.off()
        
        #survivalROC验证####
        print('start survivalROC')
        #训练集
        riskdata <- read.table('risk.txt',header = T,sep = '\t',check.names = F)
        mycol <- c('OrangeRed','Gold','GreenYellow','LightSlateBlue','DarkTurquoise')
        roc <- survivalROC(Stime = riskdata$futime,
                           status = riskdata$fustate,
                           marker = riskdata$riskScore,
                           predict.time = 365,method='KM')
        aucText=paste0("AUC for 1-Year Survival: ",sprintf("%.3f",roc$AUC))
        #pdf('Fig1K.train.ROC.pdf',width = 6.5,height = 7)
        plot(roc$FP,roc$TP, type="l",lwd = 3,xlim=c(0,1), ylim=c(0,1.1),
             xlab = "False Positive rate",
             ylab = "True Positive rate", col="OrangeRed",
             main = 'SurvivalROC for Training Set')
        abline(0,1)
        for (i in c(2:5)){
          roc1 <- survivalROC(Stime = riskdata$futime,
                              status = riskdata$fustate,
                              marker = riskdata$riskScore,
                              predict.time = 365*i,method='KM')
          aucText=c(aucText,paste0("AUC for ",i,"-Year Survival:",sprintf("%.3f",roc1$AUC)))
          lines(roc1$FP, roc1$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=mycol[i],lwd = 3)
        }
        legend("bottomright", aucText,lwd=2,bty="n",col=mycol[1:5])
        dev.off()
        AUC1 <- aucText
        
        ###验证集
        riskdata <- read.table('test.risk.txt',header = T,sep = '\t',check.names = F)
        mycol <- c('OrangeRed','Gold','GreenYellow','LightSlateBlue','DarkTurquoise')
        roc <- survivalROC(Stime = riskdata$futime,
                           status = riskdata$fustate,
                           marker = riskdata$riskScore,
                           predict.time = 365,method='KM')
        aucText=paste0("AUC for 1-Year Survival: ",sprintf("%.3f",roc$AUC))
        #pdf('Fig1L.test.ROC.pdf',width = 6.5,height = 7)
        plot(roc$FP,roc$TP, type="l",lwd = 3,xlim=c(0,1), ylim=c(0,1.1),
             xlab = "False Positive rate",
             ylab = "True Positive rate", col="OrangeRed",
             main = 'SurvivalROC for Testing Set')
        abline(0,1)
        for (i in c(2:5)){
          roc1 <- survivalROC(Stime = riskdata$futime,
                              status = riskdata$fustate,
                              marker = riskdata$riskScore,
                              predict.time = 365*i,method='KM')
          lines(roc1$FP, roc1$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=mycol[i],lwd = 3)
          aucText=c(aucText,paste0("AUC for ",i,"-Year Survival:",sprintf("%.3f",roc1$AUC)))
        }
        legend("bottomright", aucText,lwd=2,bty="n",col=mycol[1:5])
        dev.off()
        AUC2 <- aucText
        
        #独立预后,用所有的样本做####
        print('start indepedent')
        #变成数字，Patient用OSinput输出补全
        input <- read.table('survivalInput2.txt',header = T,sep = '\t',check.names = F,fill = NA)
        rownames(input) <- input$Patient
        cox.train <- read.table('train.coxinput.txt',header = T,sep = '\t',check.names = F) %>% .$Patient
        cox.vali <- read.table('test.coxinput.txt',header = T,sep = '\t',check.names = F) %>% .$Patient
        train <- input[cox.train,]
        vali <- input[cox.vali,] 
        risk.train <-read.table('risk.txt',header = T,sep = '\t',check.names = F) %>% .[,c('Patient','riskScore')]
        risk.vali <-read.table('test.risk.txt',header = T,sep = '\t',check.names = F) %>% .[,c('Patient','riskScore')]
        train <- merge(risk.train,train,by='Patient',sort=F)
        train <- train[,c(1,3,4,2,5:10)]
        vali <- merge(risk.vali,vali,by='Patient',sort=F)
        vali <- vali[,c(1,3,4,2,5:10)]
        ##训练集，单因素Cox检验
        #rtttt = train[,-c(1,2)]
        rtttt = rbind(train,vali)
        outTab2 = data.frame()
        for(i in colnames(rtttt[,4:ncol(rtttt)])){
          cox <- coxph(Surv(futime, fustate) ~ rtttt[,i], data = rtttt)
          coxSummary = summary(cox)
          coxP=coxSummary$coefficients[,"Pr(>|z|)"]
          outTab2=rbind(outTab2,
                        cbind(id=i,
                              HR=coxSummary$conf.int[,"exp(coef)"],
                              HR.95L=coxSummary$conf.int[,"lower .95"],
                              HR.95H=coxSummary$conf.int[,"upper .95"],
                              pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
          )
        }
        #outTab2$pvalue <- sapply(outTab2$pvalue,function(x){as.numeric(x)})
        outTab2[,c(2:5)] <- apply(outTab2[,c(2:5)],2,function(x){as.numeric(x)})
        outTab2 <- outTab2[order(outTab2$pvalue),]
        write.table(outTab2,file="indep.uniCox.xls",sep="\t",row.names=F,quote=F)
        tabletext <- outTab2
        tabletext[,2:ncol(tabletext)] <- apply(tabletext[,2:ncol(tabletext)],2,
                                               function(x){format(as.numeric(x),digits=4)})
        tabletext <- rbind(c("Variable","HR","lower 95%CI","upper 95%CI","pvalue"),tabletext)
        pdf('Fig3A.indep.unicox.pdf',width = 10,height =6)
        forestplot(labeltext=tabletext,
                   mean=c(NA,outTab2$HR),#HR值
                   lower=c(NA,outTab2$HR.95L), #95%置信区间下限
                   upper=c(NA,outTab2$HR.95H),#95%置信区间上限
                   graph.pos=6,#图在表中的列位置
                   graphwidth = unit(.25,"npc"),#图在表中的宽度比
                   fn.ci_norm="fpDrawDiamondCI",#box类型选择钻石
                   col=fpColors(box="#00A896", lines="#02C39A", zero = "black"),#box颜色
                   boxsize=0.2,#box大小固定
                   lwd.ci=1.6,
                   ci.vertices.height = 0.2,ci.vertices=T,#不显示区间
                   zero=1,#zero线横坐标
                   lwd.zero=2,#zero线宽
                   xticks = c(0.1,1,12),#横坐标刻度根据需要可随意设置
                   lwd.xaxis=2,
                   xlab="HR",
                   hrzl_lines=list("1" = gpar(lwd=2, col="black"),#第一行顶部加黑实线
                                   "2" = gpar(lwd=1, col="grey50", lty=2),#第二行顶部加灰色虚线
                                   "9" = gpar(lwd=2, col="black")),#最后一行底部加黑线，""中数字为nrow(tabletext) + 1
                   txt_gp=fpTxtGp(label=gpar(cex=1.2),#各种字体大小设置
                                  ticks=gpar(cex=0.85),
                                  xlab=gpar(cex=1.2),
                                  title=gpar(cex=1.2)),
                   lineheight = unit(1.7,"cm"),#固定行高
                   colgap = unit(0.5,"cm"),
                   mar=unit(rep(1.5, times = 4), "cm"),
                   new_page = F
        )
        dev.off()
        
        #多因素COX
        indep.uni <- outTab2[outTab2$pvalue<=0.05,]$id
        rtttt1 <- rtttt[,c('futime','fustate',unique(indep.uni))]
        multiCox1=coxph(Surv(futime, fustate) ~ ., data = rtttt1)
        multiCox1=step(multiCox1,direction = "both")
        multiCoxSum1=summary(multiCox1)
        outTab3=data.frame()
        outTab3=cbind(
          HR=multiCoxSum1$conf.int[,"exp(coef)"],
          HR.95L=multiCoxSum1$conf.int[,"lower .95"],
          HR.95H=multiCoxSum1$conf.int[,"upper .95"],
          pvalue=multiCoxSum1$coefficients[,"Pr(>|z|)"])
        outTab3=cbind(id=row.names(outTab3),outTab3)
        outTab3=data.frame(outTab3)
        outTab3[,c(2:5)] <- apply(outTab3[,c(2:5)],2,function(x){as.numeric(x)})
        #outTab3 <- outTab3[order(outTab3$pvalue,decreasing = F),]
        indep.final <- outTab3[outTab3$pvalue<0.05,]$id
        write.table(outTab3,'indep.mulCox.xls',row.names = F,col.names = T,sep = '\t',quote = F)
        tabletext <- outTab3
        tabletext[,2:ncol(tabletext)] <- apply(tabletext[,2:ncol(tabletext)],2,
                                               function(x){format(as.numeric(x),digits=4)})
        tabletext <- rbind(c("Variable","HR","lower 95%CI","upper 95%CI","pvalue"),tabletext)
        pdf('Fig3B.indep.multiox.pdf',width = 10,height =5)
        forestplot(labeltext=tabletext,
                   mean=c(NA,outTab3$HR),#HR值
                   lower=c(NA,outTab3$HR.95L), #95%置信区间下限
                   upper=c(NA,outTab3$HR.95H),#95%置信区间上限
                   graph.pos=6,#图在表中的列位置
                   graphwidth = unit(.25,"npc"),#图在表中的宽度比
                   fn.ci_norm="fpDrawDiamondCI",#box类型选择钻石
                   col=fpColors(box="#00A896", lines="#02C39A", zero = "black"),#box颜色
                   boxsize=0.2,#box大小固定
                   lwd.ci=1.6,
                   ci.vertices.height = 0.2,ci.vertices=T,#不显示区间
                   zero=1,#zero线横坐标
                   lwd.zero=2,#zero线宽
                   xticks = c(1,3),#横坐标刻度根据需要可随意设置
                   lwd.xaxis=2,
                   xlab="HR",
                   hrzl_lines=list("1" = gpar(lwd=2, col="black"),#第一行顶部加黑实线
                                   "2" = gpar(lwd=1, col="grey50", lty=2),#第二行顶部加灰色虚线
                                   "5" = gpar(lwd=2, col="black")),#最后一行底部加黑线，""中数字为nrow(tabletext) + 1
                   txt_gp=fpTxtGp(label=gpar(cex=1.2),#各种字体大小设置
                                  ticks=gpar(cex=0.85),
                                  xlab=gpar(cex=1.2),
                                  title=gpar(cex=1.2)),
                   lineheight = unit(2,"cm"),#固定行高
                   colgap = unit(0.5,"cm"),
                   mar=unit(rep(1.5, times = 4), "cm"),
                   new_page = F
        )
        dev.off()
        
        #输出一堆结果
        pvalue.filt <- paste(ifelse(outTab1$pvalue > 0.05,"no","yes"),sep = "",collapse=",")
        HR.model.filt <- data.frame(HR.mul.model.f = paste(ifelse(outTab1$HR > 0.9 & outTab1$HR<1.1,'No','Yeah'),sep = "",collapse=","),
                                    HR.indep.mul.model.f = paste(ifelse(outTab3$HR > 0.9 & outTab3$HR<1.1,'No','Yeah'),sep = "",collapse=","))
        HR.model <- ifelse('riskScore' %in% outTab3$id,outTab3[outTab3$id=='riskScore',]$HR,'NULL')
        #p.indep <- paste(outTab3$pvalue,sep = "",collapse=",")
        model.info <- rbind(model.info,cbind(data.frame(Seed = seed,
                                                        MG=paste(modelgene,sep = "",collapse=","),
                                                        train.p.sur=pValue1,
                                                        value.p.sur=pValue2,
                                                        train.auc=paste(AUC1,sep = "",collapse=","),
                                                        vali.auc=paste(AUC2,sep = "",collapse=","),
                                                        Indep.factor=paste(indep.final,sep = "",collapse=","),
                                                        p.filt = pvalue.filt,
                                                        HR.riskScore = HR.model)))
        
      }
    }
  }
}
write.table(model.info,'model.info.txt',row.names = F,sep = '\t',quote = F)

#训练集里模型snoRNA的生存曲线####
survival <- coxinput[,c('Patient','fustate','futime','U81','SNORA7B')]
survival$exp1 <- ifelse(survival[,4]>median(survival[,4]),
                       'High_Expression','Low_Expression')
survival$exp2 <- ifelse(survival[,5]>median(survival[,5]),
                         'High_Expression','Low_Expression')
fit <- survfit(Surv(futime, fustate) ~ exp1, data = survival)
pdf('Fig1K.U81.survival.pdf',width = 6,height = 6)
ggsurvplot(fit, data = survival ,
                ggtheme = theme_bw(), #想要网格就运行这行
                conf.int = F, #不画置信区间，想画置信区间就把F改成T
                #conf.int.style = "step",#置信区间的类型，还可改为ribbon
                risk.table = TRUE, # 绘制累计风险曲线
                #tables.theme = theme_void(),
                tables.height = 0.25,
                censor = T, #不显示观察值所在的位置
                palette = c("DodgerBlue4","Red3"), #线的颜色对应高、低(自定义调色板)
                xlab = "Time (Days)",
                legend.title = "U81",#基因名写在图例题目的位置
                font.legend = 12,#图例的字体大小
                #font.title = 12,font.x = 10,font.y = 10,#设置其他字体大小
                legend=c(0.8,0.88),
                legend.labs = c('High_Expression','Low_Expression'),pval = T)
dev.off()

fit <- survfit(Surv(futime, fustate) ~ exp2, data = survival)
pdf('Fig1L.SNORA7B.survival.pdf',width = 6,height = 6)
ggsurvplot(fit, data = survival ,
           ggtheme = theme_bw(), #想要网格就运行这行
           conf.int = F, #不画置信区间，想画置信区间就把F改成T
           #conf.int.style = "step",#置信区间的类型，还可改为ribbon
           risk.table = TRUE, # 绘制累计风险曲线
           #tables.theme = theme_void(),
           tables.height = 0.25,
           censor = T, #不显示观察值所在的位置
           palette = c("DodgerBlue4","Red3"), #线的颜色对应高、低(自定义调色板)
           xlab = "Time (Days)",
           legend.title = "SNORA7B",#基因名写在图例题目的位置
           font.legend = 12,#图例的字体大小
           #font.title = 12,font.x = 10,font.y = 10,#设置其他字体大小
           legend=c(0.8,0.88),
           legend.labs = c('High_Expression','Low_Expression'),pval = T)
dev.off()

#风险模型与临床性状相关性####
#卡方检验
#train是包含高低风险分组和临床信息的表格
input <- read.table('OSinput.txt',header = T,sep = '\t',check.names = F,fill = NA)
rownames(input) <- input$Patient
cox.train <- read.table('train.coxinput.txt',header = T,sep = '\t',check.names = F) %>% .$Patient
cox.vali <- read.table('test.coxinput.txt',header = T,sep = '\t',check.names = F) %>% .$Patient
train <- input[cox.train,]
vali <- input[cox.vali,] 
risk.train <-read.table('risk.txt',header = T,sep = '\t',check.names = F) %>% .[,c('Patient','risk')]
risk.vali <-read.table('test.risk.txt',header = T,sep = '\t',check.names = F) %>% .[,c('Patient','risk')]
train <- merge(risk.train,train,by='Patient',sort=F)
vali <- merge(risk.vali,vali,by='Patient',sort=F)
##训练集，单因素Cox检验
#rtttt = train[,-c(1,2)]
rtttt = rbind(train[,1:10],vali[,1:10])
riskTime <- rtttt
riskTime$risk <- factor(riskTime$risk ,levels = c("High_risk","Low_risk","P-value"),labels = c("High_risk","Low_risk","P-value"))
#分类变量标签调整
riskTime$Gender<-factor(riskTime$Gender,levels=c("female","male"),labels =c("female","male"))
riskTime$Age<- factor(ifelse(riskTime$Age>=50,">=50","<50"),labels = c(">=50","<50"))
riskTime$ajcc_pathologic_m<- factor(riskTime$ajcc_pathologic_m)
riskTime$ajcc_pathologic_n<- factor(riskTime$ajcc_pathologic_n)
# riskTime$ajcc_pathologic_stage<- factor(ifelse(riskTime$ajcc_pathologic_stage=='Stage I',"StageI~StageII",
#                                                ifelse(riskTime$ajcc_pathologic_stage=='Stage II',"StageI~StageII","StageIII~StageIII"))
#                                         )
riskTime$ajcc_pathologic_stage <- factor(riskTime$ajcc_pathologic_stage)
riskTime$ajcc_pathologic_t<- factor(riskTime$ajcc_pathologic_t)

## 左侧标签名调整
labels <- list(
  variables=list(Gender="Gender",
                 Age="Age (years)",
                 ajcc_pathologic_m="M stage",
                 ajcc_pathologic_n="N stage",
                 ajcc_pathologic_stage='Pathologic stage',
                 ajcc_pathologic_t='T stage'),
  groups=list("", "Expression",""))

## 设置短横线亚组
strata <- c(list(Total=riskTime), split(riskTime, riskTime$risk))

### 共有三处变量需要修改
rndr <- function(x, name, ...) {
  if (length(x) == 0) {
    y <- riskTime[[name]]##修改riskTime
    s <- rep("", length(render.default(x=y, name=name, ...)))
    if (is.numeric(y)) {
      p <- t.test(y ~ riskTime$risk)$p.value##修改riskTime
    } else {
      p <- chisq.test(table(y, droplevels(riskTime$risk)))$p.value###修改riskTime
    }
    s[2] <- sub("<", "&lt;", format.pval(p, digits=3, eps=0.001))
    s
  } else {
    render.default(x=x, name=name, ...)
  }
}

rndr.strat <- function(label, n, ...) {
  ifelse(n==0, label, render.strat.default(label, n, ...))
}

## 绘制三线表

table1(strata, labels, groupspan=c(1,2,1),droplevels = F,
       render=rndr, render.strat=rndr.strat,
       topclass="Rtable1-zebra")


#热图
df <- cbind(rbind(train[,1:10],vali[,1:10]),
            rbind(train[,modelgene],vali[,modelgene]))
rownames(df) <- df$Patient
df <- df[order(df$risk,decreasing = T),]
exp <- data.frame(t(df[,modelgene]),check.names = F)
annotation_col = data.frame(Age=df$Age,
                            Gender = df$Gender,
                            pStage = df$ajcc_pathologic_stage,
                            pM = df$ajcc_pathologic_m,
                            pN = df$ajcc_pathologic_n,
                            pT = df$ajcc_pathologic_t,
                            risk = factor(df$risk))   ##分组，legend名称：Type
rownames(annotation_col) = rownames(df)
ann_colors = list(Age=c(`<50`='CadetBlue1',`>=50`='DodgerBlue'),
                  Gender = c(female='Gold',male='Red'),
                  pStage= c(`Stage I`='LightSteelBlue', `Stage II`='CornflowerBlue',
                            `Stage III`='SteelBlue',`Stage IV`='NavyBlue'),
                  pN = c(N0='DarkOliveGreen2',N1='GreenYellow',
                         N2='LimeGreen',N3='ForestGreen'),
                  pT = c(T1='RosyBrown1',T2='LightCoral',
                         T3='MediumPurple',T4='DarkOrchid'),
                  pM = c(M0='Orange1',M1='DarkOrange1'),
                  risk = c(High_risk = "Brown1",Low_risk = "GreenYellow")) #给分组设置颜色.注意格式
color.key <- c("#3300CC","#3399FF","white","#FF3333","#CC0000")
pdf(file = 'Fig2B.risk.clinical.heatmap.pdf',width = 10,height = 13)  ##实际作图时，先筛差异基因，再用差异基因画图；或者用变化大的Top几千个基因画图
pheatmap(exp,cellwidth = 0.5, cellheight = 80, fontsize = 14,
         method="spearman", #计算gene或sample之间的相关性的方法，可选"pearson" (default), "kendall", or "spearman"
         scale="row", #为基因做scale
         cluster_rows=T,#为基因做聚类##不要给样本做聚类，否则样本名顺序会变
         cluster_cols=F,#不为sample做聚类
         color = colorRampPalette(color.key)(100),###表征表达量高低的颜色域
         show_colnames=F,show_rownames =T,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,##多个分组时再加上这个参数吧
         #treeheight_row = "0",treeheight_col = "0",#不画树
         border_color = "NA")
dev.off()

#sno表达量和临床性状相关性,结果不好，看看riskscore####
#U81####
colnames(df)[6:9] <- c('pM','pN','pStage','pT')
df$Age<- factor(ifelse(df$Age>=50,">=50","<50"),labels = c(">=50","<50"))
#性别
p1 <- ggboxplot(df,x='Gender',y='U81',color ='Gender',
               ylab = 'Relative Expression',
               xlab ='Gender',palette = "jco",
               add ='jitter',
               size =1)+
  theme(axis.text.x=element_text(size=13),
        axis.text.y=element_text(size=13),
        axis.title = element_text(size=15),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15),
        axis.line = element_line(size=1))+
  stat_compare_means(size = 6,aes(group=Gender),
                     label = 'p.signif',label.x = 1.5)

#年龄
p2 <- ggboxplot(df,x='Age',y='U81',color ='Age',
                ylab = 'Relative Expression',
                xlab ='Age',palette = "jco",
                add ='jitter',
                size =1)+
  theme(axis.text.x=element_text(size=13),
        axis.text.y=element_text(size=13),
        axis.title = element_text(size=15),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15),
        axis.line = element_line(size=1))+
  stat_compare_means(size = 6,aes(group=Age),
                     label = 'p.signif',label.x = 1.5)
#M分期
p3 <- ggboxplot(df,x='pM',y='U81',color ='pM',
                ylab = 'Relative Expression',
                xlab ='pM',palette = "jco",
                add ='jitter',
                size =1)+
  theme(axis.text.x=element_text(size=13),
        axis.text.y=element_text(size=13),
        axis.title = element_text(size=15),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15),
        axis.line = element_line(size=1))+
  stat_compare_means(size = 6,aes(group=pM),
                     label = 'p.signif',label.x = 1.5)

#N分期
df$pN <- factor(df$pN,levels = c('N0','N1','N2','N3'))
group <- list(c('N0','N1'),
              c('N0','N2'),
              c('N0','N3'),
              c('N1','N2'),
              c('N1','N3'),
              c('N2','N3'))
p4 <- ggboxplot(df,x='pN',y='U81',color ='pN',
                ylab = 'Relative Expression',palette = "jco",
                add ='jitter',size =1,axis.line =1)+
  stat_compare_means(size = 4,label = 'p.signif',
                     comparisons = group)+
  stat_compare_means(size=4,method = 'anova',label.x = 2,label.y = 21)+
  theme(legend.text = element_text(size =12),
        legend.title = element_text(size=12))

#T分期
df$pT <- factor(df$pT,levels = c('T1','T2','T3','T4'))
group <- list(c('T1','T2'),
              c('T1','T3'),
              c('T1','T4'),
              c('T2','T3'),
              c('T2','T4'),
              c('T3','T4'))
p5 <- ggboxplot(df,x='pT',y='U81',color ='pT',
                ylab = 'Relative Expression',palette = "jco",
                add ='jitter',size =1,axis.line =1)+
  stat_compare_means(size = 4,label = 'p.signif',
                     comparisons = group)+
  stat_compare_means(size=4,method = 'anova',label.x = 2,label.y = 21)+
  theme(legend.text = element_text(size =12),
        legend.title = element_text(size=12))

#STAGE
df$pStage <- factor(df$pStage,levels = c('Stage I','Stage II','Stage III','Stage IV'))
group <- list(c('Stage I','Stage II'),
              c('Stage I','Stage III'),
              c('Stage I','Stage IV'),
              c('Stage II','Stage III'),
              c('Stage II','Stage IV'),
              c('Stage III','Stage IV'))
p6 <- ggboxplot(df,x='pStage',y='U81',color ='pStage',
                ylab = 'Relative Expression',palette = "jco",
                add ='jitter',size =1,axis.line =1)+
  stat_compare_means(size = 4,label = 'p.signif',
                     comparisons = group)+
  stat_compare_means(size=4,method = 'anova',label.x = 2,label.y = 21)+
  theme(legend.text = element_text(size =12),
        legend.title = element_text(size=12))
pdf('Fig2D.U81.clinical.pdf',width = 15,height = 11)
p1+p2+p3+p4+p5+p6
dev.off()


#SNORA7B####
colnames(df)[6:9] <- c('pM','pN','pStage','pT')
df$Age<- factor(ifelse(df$Age>=50,">=50","<50"),labels = c(">=50","<50"))
#性别
p1 <- ggboxplot(df,x='Gender',y='SNORA7B',color ='Gender',
                ylab = 'Relative Expression',
                xlab ='Gender',palette = "jco",
                add ='jitter',
                size =1)+
  theme(axis.text.x=element_text(size=13),
        axis.text.y=element_text(size=13),
        axis.title = element_text(size=15),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15),
        axis.line = element_line(size=1))+
  stat_compare_means(size = 6,aes(group=Gender),
                     label = 'p.signif',label.x = 1.5)

#年龄
p2 <- ggboxplot(df,x='Age',y='SNORA7B',color ='Age',
                ylab = 'Relative Expression',
                xlab ='Age',palette = "jco",
                add ='jitter',
                size =1)+
  theme(axis.text.x=element_text(size=13),
        axis.text.y=element_text(size=13),
        axis.title = element_text(size=15),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15),
        axis.line = element_line(size=1))+
  stat_compare_means(size = 6,aes(group=Age),
                     label = 'p.signif',label.x = 1.5)
#M分期
p3 <- ggboxplot(df,x='pM',y='SNORA7B',color ='pM',
                ylab = 'Relative Expression',
                xlab ='pM',palette = "jco",
                add ='jitter',
                size =1)+
  theme(axis.text.x=element_text(size=13),
        axis.text.y=element_text(size=13),
        axis.title = element_text(size=15),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15),
        axis.line = element_line(size=1))+
  stat_compare_means(size = 6,aes(group=pM),
                     label = 'p.signif',label.x = 1.5)

#N分期
df$pN <- factor(df$pN,levels = c('N0','N1','N2','N3'))
group <- list(c('N0','N1'),
              c('N0','N2'),
              c('N0','N3'),
              c('N1','N2'),
              c('N1','N3'),
              c('N2','N3'))
p4 <- ggboxplot(df,x='pN',y='SNORA7B',color ='pN',
                ylab = 'Relative Expression',palette = "jco",
                add ='jitter',size =1,axis.line =1)+
  stat_compare_means(size = 4,label = 'p.signif',
                     comparisons = group)+
  stat_compare_means(size=4,method = 'anova',label.x = 2,label.y = 21)+
  theme(legend.text = element_text(size =12),
        legend.title = element_text(size=12))

#T分期
df$pT <- factor(df$pT,levels = c('T1','T2','T3','T4'))
group <- list(c('T1','T2'),
              c('T1','T3'),
              c('T1','T4'),
              c('T2','T3'),
              c('T2','T4'),
              c('T3','T4'))
p5 <- ggboxplot(df,x='pT',y='SNORA7B',color ='pT',
                ylab = 'Relative Expression',palette = "jco",
                add ='jitter',size =1,axis.line =1)+
  stat_compare_means(size = 4,label = 'p.signif',
                     comparisons = group)+
  stat_compare_means(size=4,method = 'anova',label.x = 2,label.y = 21)+
  theme(legend.text = element_text(size =12),
        legend.title = element_text(size=12))

#STAGE
df$pStage <- factor(df$pStage,levels = c('Stage I','Stage II','Stage III','Stage IV'))
group <- list(c('Stage I','Stage II'),
              c('Stage I','Stage III'),
              c('Stage I','Stage IV'),
              c('Stage II','Stage III'),
              c('Stage II','Stage IV'),
              c('Stage III','Stage IV'))
p6 <- ggboxplot(df,x='pStage',y='SNORA7B',color ='pStage',
                ylab = 'Relative Expression',palette = "jco",
                add ='jitter',size =1,axis.line =1)+
  stat_compare_means(size = 4,label = 'p.signif',
                     comparisons = group)+
  stat_compare_means(size=4,method = 'anova',label.x = 2,label.y = 21)+
  theme(legend.text = element_text(size =12),
        legend.title = element_text(size=12))
pdf('Fig2C.SNORA7B.clinical.pdf',width = 15,height = 11)
p1+p2+p3+p4+p5+p6
dev.off()

#riskscore####
input <- read.table('OSinput.txt',header = T,sep = '\t',check.names = F,fill = NA)
rownames(input) <- input$Patient
cox.train <- read.table('train.coxinput.txt',header = T,sep = '\t',check.names = F) %>% .$Patient
cox.vali <- read.table('test.coxinput.txt',header = T,sep = '\t',check.names = F) %>% .$Patient
train <- input[cox.train,]
vali <- input[cox.vali,] 
risk.train <-read.table('risk.txt',header = T,sep = '\t',check.names = F) %>% .[,c('Patient','riskScore')]
risk.vali <-read.table('test.risk.txt',header = T,sep = '\t',check.names = F) %>% .[,c('Patient','riskScore')]
train <- merge(risk.train,train,by='Patient',sort=F)
vali <- merge(risk.vali,vali,by='Patient',sort=F)
df <- cbind(rbind(train[,1:10],vali[,1:10]),
            rbind(train[,modelgene],vali[,modelgene]))
rownames(df) <- df$Patient
df <- df[order(df$risk,decreasing = T),]
colnames(df)[6:9] <- c('pM','pN','pStage','pT')
df$Age<- factor(ifelse(df$Age>=50,">=50","<50"),labels = c(">=50","<50"))
#性别
p1 <- ggboxplot(df,x='Gender',y='riskScore',color ='Gender',
                ylab = 'riskScore',
                xlab ='Gender',palette = "jco",
                add ='jitter',
                size =1)+
  theme(axis.text.x=element_text(size=13),
        axis.text.y=element_text(size=13),
        axis.title = element_text(size=15),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15),
        axis.line = element_line(size=1))+
  stat_compare_means(size = 6,aes(group=Gender),
                     label = 'p.signif',label.x = 1.5)

#年龄
p2 <- ggboxplot(df,x='Age',y='riskScore',color ='Age',
                ylab = 'riskScore',
                xlab ='Age',palette = "jco",
                add ='jitter',
                size =1)+
  theme(axis.text.x=element_text(size=13),
        axis.text.y=element_text(size=13),
        axis.title = element_text(size=15),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15),
        axis.line = element_line(size=1))+
  stat_compare_means(size = 6,aes(group=Age),
                     label = 'p.signif',label.x = 1.5)
#M分期
p3 <- ggboxplot(df,x='pM',y='riskScore',color ='pM',
                ylab = 'riskScore',
                xlab ='pM',palette = "jco",
                add ='jitter',
                size =1)+
  theme(axis.text.x=element_text(size=13),
        axis.text.y=element_text(size=13),
        axis.title = element_text(size=15),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15),
        axis.line = element_line(size=1))+
  stat_compare_means(size = 6,aes(group=pM),
                     label = 'p.signif',label.x = 1.5)

#N分期
df$pN <- factor(df$pN,levels = c('N0','N1','N2','N3'))
group <- list(c('N0','N1'),
              c('N0','N2'),
              c('N0','N3'),
              c('N1','N2'),
              c('N1','N3'),
              c('N2','N3'))
p4 <- ggboxplot(df,x='pN',y='riskScore',color ='pN',
                ylab = 'riskScore',palette = "jco",
                add ='jitter',size =1,axis.line =1)+
  stat_compare_means(size = 4,label = 'p.signif',
                     comparisons = group)+
  stat_compare_means(size=4,method = 'anova',label.x = 2,label.y = 9.6)+
  theme(legend.text = element_text(size =15),
        legend.title = element_text(size=15))

#T分期
df$pT <- factor(df$pT,levels = c('T1','T2','T3','T4'))
group <- list(c('T1','T2'),
              c('T1','T3'),
              c('T1','T4'),
              c('T2','T3'),
              c('T2','T4'),
              c('T3','T4'))
p5 <- ggboxplot(df,x='pT',y='riskScore',color ='pT',
                ylab = 'riskScore',palette = "jco",
                add ='jitter',size =1,axis.line =1)+
  stat_compare_means(size = 4,label = 'p.signif',
                     comparisons = group)+
  stat_compare_means(size=4,method = 'anova',label.x = 2,label.y = 10)+
  theme(legend.text = element_text(size =15),
        legend.title = element_text(size=15))

#STAGE
df$pStage <- factor(df$pStage,levels = c('Stage I','Stage II','Stage III','Stage IV'))
group <- list(c('Stage I','Stage II'),
              c('Stage I','Stage III'),
              c('Stage I','Stage IV'),
              c('Stage II','Stage III'),
              c('Stage II','Stage IV'),
              c('Stage III','Stage IV'))
p6 <- ggboxplot(df,x='pStage',y='riskScore',color ='pStage',
                ylab = 'riskScore',palette = "jco",
                add ='jitter',size =1,axis.line =1)+
  stat_compare_means(size = 4,label = 'p.signif',
                     comparisons = group)+
  stat_compare_means(size=4,method = 'anova',label.x = 2,label.y = 9)+
  theme(legend.text = element_text(size =15),
        legend.title = element_text(size=15))
pdf('Fig2C.riskScore.clinical.pdf',width = 14,height = 11)
p1+p2+p3+p4+p5+p6
dev.off()




#列线图####
rtttt = rbind(train,vali)
riskmerge = rtttt[,-c(5,6,8,10)]
#riskmerge <- riskmerge[order(riskmerge$risk,decreasing = T),]
nomod <- riskmerge
dd<-datadist(nomod)
options(datadist="dd")
options(na.action="na.delete")
summary(nomod$futime)#这里的futime要和nomod一致

coxpbc<-cph(formula = Surv(futime,fustate) ~  riskScore+ajcc_pathologic_stage+Age,
            data = nomod,x=T,y=T,surv = T,na.action=na.delete)  #,time.inc =2920
print(coxpbc)
surv<-Survival(coxpbc)
surv1 <- function(x) surv(365,x)
surv2 <- function(x) surv(365*2,x)
surv3<-function(x) surv(1095,x)
surv4<-function(x) surv(365*4,x)
surv5<-function(x) surv(1825,x)
#surv10<-function(x) surv(3650,x)
x<-nomogram(coxpbc,fun = list(surv1,surv2,surv3,surv4,surv5),lp=T,
            funlabel = c('1-year survival Probability',
                         '2-year survival Probability',
                         '3-year survival Probability',
                         '4-year survival Probability',
                         '5-year survival Probability'),
            maxscale = 100,fun.at = c(0.95,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1))
pdf("Fig3C.indep_nomogram_classical.pdf",width = 15,height = 9)
# png("train.indep_nomogram_classical.png",width = 1200,height = 600)
plot(x, lplabel="Linear Predictor",
     xfrac=.35,varname.label=TRUE, varname.label.sep="=", ia.space=.1,
     tck=NA, tcl=-0.20, lmgp=0.3,
     points.label='Points', total.points.label='Total Points',
     total.sep.page=FALSE,
     cap.labels=FALSE,cex.var = 1.6,cex.axis = 1.2,
     lwd=5,label.every = 1,col.grid = gray(c(0.8, 0.95)))
dev.off()

###绘制calibration curve进行验证
#1年数据
#set.seed(12345)
set.seed(888666)
#set.seed(0622)
f1 <- cph(formula = Surv(futime,fustate) ~ riskScore+ajcc_pathologic_stage+Age,
          data=nomod,x=T,y=T,surv = T,na.action=na.delete,time.inc = 365)
#参数m=50表示每组50个样本进行重复计算
cal1 <- calibrate(f1, cmethod="KM", method="boot",u=365,m=125,B=1000)
#2年数据
set.seed(12345)
f2 <- cph(formula = Surv(futime,fustate) ~ riskScore+ajcc_pathologic_stage+Age,
          data=nomod,x=T,y=T,surv = T,na.action=na.delete,time.inc = 365*2)
cal2 <- calibrate(f2, cmethod="KM", method="boot",u=365*2,m=125,B=1000)
#3年数据
f3 <- cph(formula = Surv(futime,fustate) ~  riskScore+ajcc_pathologic_stage+Age,
          data=nomod,x=T,y=T,surv = T,na.action=na.delete,time.inc = 1095)
cal3 <- calibrate(f3, cmethod="KM", method="boot",u=1095,m=125,B=1000)
#4年数据
f4 <- cph(formula = Surv(futime,fustate) ~  riskScore+ajcc_pathologic_stage+Age,
          data=nomod,x=T,y=T,surv = T,na.action=na.delete,time.inc = 365*4)
cal4 <- calibrate(f4, cmethod="KM", method="boot",u=365*4,m=125,B=1000)
#5年数据
f5 <- cph(formula = Surv(futime,fustate) ~  riskScore+ajcc_pathologic_stage+Age,
          data=nomod,x=T,y=T,surv = T,na.action=na.delete,time.inc = 1825)
cal5 <- calibrate(f5, cmethod="KM", method="boot",u=1825,m=125,B=1000)

mycol <- c('OrangeRed','Gold','GreenYellow','LightSlateBlue','DarkTurquoise')
pdf("Fig3D.train.indep_calibration_compare.pdf",width = 7,height = 7)
plot(cal1,lwd = 2,lty = 1,errbar.col = c(mycol[1]),#lty=1画上bar
     bty = "l", #只画左边和下边框
     xlim = c(0,1),ylim= c(0,1),#xlim可以调整，怎么好看怎么调
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = c(mycol[1]),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal1[,c('mean.predicted',"KM")],
      type = 'b', lwd = 2, col = c(mycol[1]), pch =16)
mtext("")
plot(cal2,lwd = 2,lty = 1,errbar.col = c(mycol[2]),
     xlim = c(0,1),ylim= c(0,1),col = c(mycol[2]),add = T)
lines(cal2[,c('mean.predicted',"KM")],
      type = 'b', lwd = 4, col = c(mycol[2]), pch = 16)
plot(cal3,lwd = 2,lty = 1,errbar.col = c(mycol[3]),
     xlim = c(0,1),ylim= c(0,1),col = c(mycol[3]),add = T)
lines(cal3[,c('mean.predicted',"KM")],
      type = 'b', lwd = 2, col = c(mycol[3]), pch = 16)
plot(cal4,lwd = 2,lty = 1,errbar.col = c(mycol[4]),
     xlim = c(0,1),ylim= c(0,1),col = c(mycol[4]),add = T)
lines(cal4[,c('mean.predicted',"KM")],
      type = 'b', lwd = 4, col = c(mycol[4]), pch = 16)
plot(cal5,lwd = 2,lty = 1,errbar.col = c(mycol[5]),
     xlim = c(0,1),ylim= c(0,1),col = c(mycol[5]),add = T)
lines(cal5[,c('mean.predicted',"KM")],
      type = 'b', lwd = 2, col = c(mycol[5]), pch = 16)

abline(0,1, lwd = 3, lty = 3, col = c("#224444"))

legend("topleft", #图例的位置
       legend = c('1-year',"2-year",'3-year',"3-year","5-year"), #图例文字
       col =mycol, #图例线的颜色，与文字对应
       lwd = 2,#图例中线的粗细
       cex = 1.2,#图例字体大小
       bty = "n")#不显示图例边框
dev.off()
#c指数
nomod <- riskmerge
f <- coxph(Surv(futime,fustate)~ riskScore+ajcc_pathologic_stage+Age,data=nomod)
sum.surv<-summary(f)
c_index<-sum.surv$concordance
c_index
#C      se(C) 
#0.75972520 0.02314671 



#2个snoRNA####
gtf <- read.table('F:/workspace/gencode.v22.annotation.gtf',header = F,
                  fill=NA,sep = '\t',check.names = F)
write.table(gtf[grepl('ENSG00000207088',gtf$V9),],'SNORA7B.txt',sep = '\t',quote = F,row.names = F)
#SNORA7B:ENSG00000207088,位置chr3:129397210-129397348

#U81:U81_chr1_173833283_173833360

#CNV数据处理####
#EXCEL

#CNV相关性分析####
#U81####
U81 <- read.table('chr1.U81.txt',header = T,
                  sep = '\t',check.names = F)
#U81exp <- TCGA_BRCA.sno.unique['U81_chr1_173833283_173833360',]
#U81exp <- data.frame(t(U81exp),check.names = F)
U81exp <- TCGA_BRCA.sno.unique.log['U81_chr1_173833283_173833360',]
U81exp <- data.frame(U81=U81exp)
U81exp$sample <- rownames(U81exp)
U81exp <- merge(U81exp,U81,by='sample')
library(gmodels)
rt <- data.frame(snoRNA=U81exp$U81,
                   CNV=U81exp$value)
p1 <- ggplot(data=rt, aes(x=snoRNA, y=CNV))+
    geom_point(color="blue")+ggtitle('U81')+
    stat_smooth(method="lm",se=FALSE,color='red')+
    stat_cor(data=rt, method = "pearson")+theme_bw()
  
#SNORA7B####
SNORA7B <- read.table('chr3.SNORA7B.txt',header = T,
                  sep = '\t',check.names = F)
#SNORA7Bexp <- TCGA_BRCA.sno.unique['ACA7B;ENSG00000207088;SNORA7B_chr3_129116052_129116191',]
#SNORA7Bexp <- data.frame(t(SNORA7Bexp),check.names = F)
SNORA7Bexp <- TCGA_BRCA.sno.unique.log['ACA7B;ENSG00000207088;SNORA7B_chr3_129116052_129116191',]
SNORA7Bexp <- data.frame(SNORA7B=SNORA7Bexp)
SNORA7Bexp$sample <- rownames(SNORA7Bexp)
SNORA7Bexp <- merge(SNORA7Bexp,SNORA7B,by='sample')
library(gmodels)
rt <- data.frame(snoRNA=SNORA7Bexp$SNORA7B,
                 CNV=SNORA7Bexp$value)
p2 <- ggplot(data=rt, aes(x=snoRNA, y=CNV))+
  geom_point(color="blue")+ggtitle('SNORA7B')+
  stat_smooth(method="lm",se=FALSE,color='red')+
  stat_cor(data=rt, method = "pearson")+theme_bw()
pdf('sno-CNV.corr.pdf',width = 7,height = 4)
p1+p2
dev.off()

#snoRNA与mRNA相关性####
exp <- read.table('mRNAmatrix.symbol.txt',header = T,sep = '\t',check.names = F)
#只要01A和11A
exp.filt <- cbind(exp$id,exp[,substr(colnames(exp),14,16)=='11A'],
                  exp[,substr(colnames(exp),14,16)=='01A'])
colnames(exp.filt)[1] <- 'id'
#去重复病人的样本样本
tcgaReplicateFilter = function(tsb, analyte_target=c("DNA","RNA"), decreasing=TRUE, analyte_position=20, plate=c(22,25), portion=c(18,19), filter_FFPE=FALSE, full_barcode=FALSE){
  # basically, user provide tsb and analyte_target is fine. If you
  # want to filter FFPE samples, please set filter_FFPE and full_barcode
  # all to TRUE, and tsb must have nchar of 28
  
  analyte_target = match.arg(analyte_target)
  # Strings in R are largely lexicographic
  # see ??base::Comparison
  
  # filter FFPE samples
  # provide by <http://gdac.broadinstitute.org/runs/sampleReports/latest/FPPP_FFPE_Cases.html> 
  if(full_barcode & filter_FFPE){
    ffpe = c("TCGA-44-2656-01B-06D-A271-08", "TCGA-44-2656-01B-06D-A273-01", 
             "TCGA-44-2656-01B-06D-A276-05", "TCGA-44-2656-01B-06D-A27C-26", 
             "TCGA-44-2656-01B-06R-A277-07", "TCGA-44-2662-01B-02D-A271-08", 
             "TCGA-44-2662-01B-02D-A273-01", "TCGA-44-2662-01B-02R-A277-07", 
             "TCGA-44-2665-01B-06D-A271-08", "TCGA-44-2665-01B-06D-A273-01", 
             "TCGA-44-2665-01B-06D-A276-05", "TCGA-44-2665-01B-06R-A277-07", 
             "TCGA-44-2666-01B-02D-A271-08", "TCGA-44-2666-01B-02D-A273-01", 
             "TCGA-44-2666-01B-02D-A276-05", "TCGA-44-2666-01B-02D-A27C-26", 
             "TCGA-44-2666-01B-02R-A277-07", "TCGA-44-2668-01B-02D-A271-08", 
             "TCGA-44-2668-01B-02D-A273-01", "TCGA-44-2668-01B-02D-A276-05", 
             "TCGA-44-2668-01B-02D-A27C-26", "TCGA-44-2668-01B-02R-A277-07", 
             "TCGA-44-3917-01B-02D-A271-08", "TCGA-44-3917-01B-02D-A273-01", 
             "TCGA-44-3917-01B-02D-A276-05", "TCGA-44-3917-01B-02D-A27C-26", 
             "TCGA-44-3917-01B-02R-A277-07", "TCGA-44-3918-01B-02D-A271-08", 
             "TCGA-44-3918-01B-02D-A273-01", "TCGA-44-3918-01B-02D-A276-05", 
             "TCGA-44-3918-01B-02D-A27C-26", "TCGA-44-3918-01B-02R-A277-07", 
             "TCGA-44-4112-01B-06D-A271-08", "TCGA-44-4112-01B-06D-A273-01", 
             "TCGA-44-4112-01B-06D-A276-05", "TCGA-44-4112-01B-06D-A27C-26", 
             "TCGA-44-4112-01B-06R-A277-07", "TCGA-44-5645-01B-04D-A271-08", 
             "TCGA-44-5645-01B-04D-A273-01", "TCGA-44-5645-01B-04D-A276-05", 
             "TCGA-44-5645-01B-04D-A27C-26", "TCGA-44-5645-01B-04R-A277-07", 
             "TCGA-44-6146-01B-04D-A271-08", "TCGA-44-6146-01B-04D-A273-01", 
             "TCGA-44-6146-01B-04D-A276-05", "TCGA-44-6146-01B-04D-A27C-26", 
             "TCGA-44-6146-01B-04R-A277-07", "TCGA-44-6146-01B-04R-A27D-13", 
             "TCGA-44-6147-01B-06D-A271-08", "TCGA-44-6147-01B-06D-A273-01", 
             "TCGA-44-6147-01B-06D-A276-05", "TCGA-44-6147-01B-06D-A27C-26", 
             "TCGA-44-6147-01B-06R-A277-07", "TCGA-44-6147-01B-06R-A27D-13", 
             "TCGA-44-6775-01C-02D-A271-08", "TCGA-44-6775-01C-02D-A273-01", 
             "TCGA-44-6775-01C-02D-A276-05", "TCGA-44-6775-01C-02D-A27C-26", 
             "TCGA-44-6775-01C-02R-A277-07", "TCGA-44-6775-01C-02R-A27D-13", 
             "TCGA-A6-2674-01B-04D-A270-10", "TCGA-A6-2674-01B-04R-A277-07", 
             "TCGA-A6-2677-01B-02D-A270-10", "TCGA-A6-2677-01B-02D-A274-01", 
             "TCGA-A6-2677-01B-02D-A27A-05", "TCGA-A6-2677-01B-02D-A27E-26", 
             "TCGA-A6-2677-01B-02R-A277-07", "TCGA-A6-2684-01C-08D-A270-10", 
             "TCGA-A6-2684-01C-08D-A274-01", "TCGA-A6-2684-01C-08D-A27A-05", 
             "TCGA-A6-2684-01C-08D-A27E-26", "TCGA-A6-2684-01C-08R-A277-07", 
             "TCGA-A6-3809-01B-04D-A270-10", "TCGA-A6-3809-01B-04D-A274-01", 
             "TCGA-A6-3809-01B-04D-A27A-05", "TCGA-A6-3809-01B-04D-A27E-26", 
             "TCGA-A6-3809-01B-04R-A277-07", "TCGA-A6-3810-01B-04D-A270-10", 
             "TCGA-A6-3810-01B-04D-A274-01", "TCGA-A6-3810-01B-04D-A27A-05", 
             "TCGA-A6-3810-01B-04D-A27E-26", "TCGA-A6-3810-01B-04R-A277-07", 
             "TCGA-A6-5656-01B-02D-A270-10", "TCGA-A6-5656-01B-02D-A274-01", 
             "TCGA-A6-5656-01B-02D-A27A-05", "TCGA-A6-5656-01B-02D-A27E-26", 
             "TCGA-A6-5656-01B-02R-A277-07", "TCGA-A6-5656-01B-02R-A27D-13", 
             "TCGA-A6-5659-01B-04D-A270-10", "TCGA-A6-5659-01B-04D-A274-01", 
             "TCGA-A6-5659-01B-04D-A27A-05", "TCGA-A6-5659-01B-04D-A27E-26", 
             "TCGA-A6-5659-01B-04R-A277-07", "TCGA-A6-6650-01B-02D-A270-10", 
             "TCGA-A6-6650-01B-02D-A274-01", "TCGA-A6-6650-01B-02D-A27A-05", 
             "TCGA-A6-6650-01B-02D-A27E-26", "TCGA-A6-6650-01B-02R-A277-07", 
             "TCGA-A6-6650-01B-02R-A27D-13", "TCGA-A6-6780-01B-04D-A270-10", 
             "TCGA-A6-6780-01B-04D-A274-01", "TCGA-A6-6780-01B-04D-A27A-05", 
             "TCGA-A6-6780-01B-04D-A27E-26", "TCGA-A6-6780-01B-04R-A277-07", 
             "TCGA-A6-6780-01B-04R-A27D-13", "TCGA-A6-6781-01B-06D-A270-10", 
             "TCGA-A6-6781-01B-06D-A274-01", "TCGA-A6-6781-01B-06D-A27A-05", 
             "TCGA-A6-6781-01B-06R-A277-07", "TCGA-A6-6781-01B-06R-A27D-13", 
             "TCGA-A7-A0DB-01C-02D-A272-09", "TCGA-A7-A0DB-01C-02R-A277-07", 
             "TCGA-A7-A0DB-01C-02R-A27D-13", "TCGA-A7-A13D-01B-04D-A272-09", 
             "TCGA-A7-A13D-01B-04R-A277-07", "TCGA-A7-A13D-01B-04R-A27D-13", 
             "TCGA-A7-A13E-01B-06D-A272-09", "TCGA-A7-A13E-01B-06R-A277-07", 
             "TCGA-A7-A13E-01B-06R-A27D-13", "TCGA-A7-A26E-01B-06D-A272-09", 
             "TCGA-A7-A26E-01B-06D-A275-01", "TCGA-A7-A26E-01B-06D-A27B-05", 
             "TCGA-A7-A26E-01B-06R-A277-07", "TCGA-A7-A26E-01B-06R-A27D-13", 
             "TCGA-A7-A26J-01B-02D-A272-09", "TCGA-A7-A26J-01B-02D-A275-01", 
             "TCGA-A7-A26J-01B-02D-A27B-05", "TCGA-A7-A26J-01B-02D-A27F-26", 
             "TCGA-A7-A26J-01B-02R-A277-07", "TCGA-A7-A26J-01B-02R-A27D-13", 
             "TCGA-B2-3923-01B-10D-A270-10", "TCGA-B2-3923-01B-10R-A277-07", 
             "TCGA-B2-3923-01B-10R-A27D-13", "TCGA-B2-3924-01B-03D-A270-10", 
             "TCGA-B2-3924-01B-03D-A274-01", "TCGA-B2-3924-01B-03D-A27A-05", 
             "TCGA-B2-3924-01B-03D-A27E-26", "TCGA-B2-3924-01B-03R-A277-07", 
             "TCGA-B2-3924-01B-03R-A27D-13", "TCGA-B2-5633-01B-04D-A270-10", 
             "TCGA-B2-5633-01B-04D-A274-01", "TCGA-B2-5633-01B-04D-A27A-05", 
             "TCGA-B2-5633-01B-04D-A27E-26", "TCGA-B2-5633-01B-04R-A277-07", 
             "TCGA-B2-5633-01B-04R-A27D-13", "TCGA-B2-5635-01B-04D-A270-10", 
             "TCGA-B2-5635-01B-04D-A274-01", "TCGA-B2-5635-01B-04D-A27A-05", 
             "TCGA-B2-5635-01B-04D-A27E-26", "TCGA-B2-5635-01B-04R-A277-07", 
             "TCGA-B2-5635-01B-04R-A27D-13", "TCGA-BK-A0CA-01B-02D-A272-09", 
             "TCGA-BK-A0CA-01B-02D-A275-01", "TCGA-BK-A0CA-01B-02D-A27B-05", 
             "TCGA-BK-A0CA-01B-02D-A27F-26", "TCGA-BK-A0CA-01B-02R-A277-07", 
             "TCGA-BK-A0CA-01B-02R-A27D-13", "TCGA-BK-A0CC-01B-04D-A272-09", 
             "TCGA-BK-A0CC-01B-04D-A275-01", "TCGA-BK-A0CC-01B-04D-A27B-05", 
             "TCGA-BK-A0CC-01B-04R-A277-07", "TCGA-BK-A0CC-01B-04R-A27D-13", 
             "TCGA-BK-A139-01C-08D-A272-09", "TCGA-BK-A139-01C-08D-A275-01", 
             "TCGA-BK-A139-01C-08D-A27B-05", "TCGA-BK-A139-01C-08D-A27F-26", 
             "TCGA-BK-A139-01C-08R-A277-07", "TCGA-BK-A139-01C-08R-A27D-13", 
             "TCGA-BK-A26L-01C-04D-A272-09", "TCGA-BK-A26L-01C-04D-A275-01", 
             "TCGA-BK-A26L-01C-04D-A27B-05", "TCGA-BK-A26L-01C-04D-A27F-26", 
             "TCGA-BK-A26L-01C-04R-A277-07", "TCGA-BK-A26L-01C-04R-A27D-13", 
             "TCGA-BL-A0C8-01B-04D-A271-08", "TCGA-BL-A0C8-01B-04D-A273-01", 
             "TCGA-BL-A0C8-01B-04D-A276-05", "TCGA-BL-A0C8-01B-04D-A27C-26", 
             "TCGA-BL-A0C8-01B-04R-A277-07", "TCGA-BL-A0C8-01B-04R-A27D-13", 
             "TCGA-BL-A13I-01B-04D-A271-08", "TCGA-BL-A13I-01B-04D-A276-05", 
             "TCGA-BL-A13I-01B-04R-A277-07", "TCGA-BL-A13I-01B-04R-A27D-13", 
             "TCGA-BL-A13J-01B-04D-A271-08", "TCGA-BL-A13J-01B-04D-A273-01", 
             "TCGA-BL-A13J-01B-04D-A276-05", "TCGA-BL-A13J-01B-04D-A27C-26", 
             "TCGA-BL-A13J-01B-04R-A277-07", "TCGA-BL-A13J-01B-04R-A27D-13")
    
    tsb = setdiff(tsb, tsb[which(tsb %in% ffpe)])
  }
  
  # find repeated samples
  sampleID = substr(tsb, start = 1, stop = 15)
  dp_samples = unique(sampleID[duplicated(sampleID)])
  
  if(length(dp_samples)==0){
    message("ooo Not find any duplicated barcodes, return original input..")
    tsb
  }else{
    uniq_tsb = tsb[! sampleID %in% dp_samples]
    dp_tsb = setdiff(tsb, uniq_tsb)
    
    add_tsb = c()
    
    # analyte = substr(dp_tsb, start = analyte_position, stop = analyte_position)
    # if analyte_target = "DNA"
    # analyte:  D > G,W,X
    if(analyte_target == "DNA"){
      for(x in dp_samples){
        mulaliquots = dp_tsb[substr(dp_tsb,1,15) == x]
        analytes = substr(mulaliquots, 
                          start = analyte_position,
                          stop = analyte_position)
        if(any(analytes == "D") & !(all(analytes == "D"))){
          aliquot = mulaliquots[which(analytes == "D")]
          
          if (length(aliquot) != 1) {
            # Still have repeats
            # Remove the samples and add repeated id back to list
            dp_tsb = c(setdiff(dp_tsb, mulaliquots), aliquot)
          } else {
            # If have no repeats
            # Just remove samples from list and
            # add unique id to result list
            add_tsb = c(add_tsb, aliquot)
            dp_tsb = setdiff(dp_tsb, mulaliquots) 
          }
        }
        
      }
    }else{
      # if analyte_target = "RNA"
      # analyte: H > R > T 
      for(x in dp_samples){
        mulaliquots = dp_tsb[substr(dp_tsb,1,15) == x]
        analytes = substr(mulaliquots, 
                          start = analyte_position,
                          stop = analyte_position)
        if(any(analytes == "H") & !(all(analytes == "H"))){
          aliquot = mulaliquots[which(analytes == "H")]
          
          if (length(aliquot) != 1) {
            # Still have repeats
            # Remove the samples and add repeated id back to list
            dp_tsb = c(setdiff(dp_tsb, mulaliquots), aliquot)
          } else {
            # If have no repeats
            # Just remove samples from list and
            # add unique id to result list
            add_tsb = c(add_tsb, aliquot)
            dp_tsb = setdiff(dp_tsb, mulaliquots) 
          }
          
        }else if(any(analytes == "R") & !(all(analytes == "R"))){
          aliquot = mulaliquots[which(analytes == "R")]
          
          if (length(aliquot) != 1) {
            # Still have repeats
            # Remove the samples and add repeated id back to list
            dp_tsb = c(setdiff(dp_tsb, mulaliquots), aliquot)
          } else {
            # If have no repeats
            # Just remove samples from list and
            # add unique id to result list
            add_tsb = c(add_tsb, aliquot)
            dp_tsb = setdiff(dp_tsb, mulaliquots) 
          }
        }else if(any(analytes == "T") & !(all(analytes == "T"))){
          aliquot = mulaliquots[which(analytes == "T")]
          
          if (length(aliquot) != 1) {
            # Still have repeats
            # Remove the samples and add repeated id back to list
            dp_tsb = c(setdiff(dp_tsb, mulaliquots), aliquot)
          } else {
            # If have no repeats
            # Just remove samples from list and
            # add unique id to result list
            add_tsb = c(add_tsb, aliquot)
            dp_tsb = setdiff(dp_tsb, mulaliquots) 
          }
        }
        
      }
    }
    
    
    if(length(dp_tsb) == 0){
      message("ooo Filter barcodes successfully!")
      c(uniq_tsb, add_tsb)
    }else{
      # filter according to portion number
      sampleID_res = substr(dp_tsb, start=1, stop=15)
      dp_samples_res = unique(sampleID_res[duplicated(sampleID_res)])
      
      for(x in dp_samples_res){
        mulaliquots = dp_tsb[substr(dp_tsb,1,15) == x]
        portion_codes = substr(mulaliquots,
                               start = portion[1],
                               stop = portion[2])
        portion_keep = sort(portion_codes, decreasing = decreasing)[1]
        if(!all(portion_codes == portion_keep)){
          if(length(which(portion_codes == portion_keep)) == 1){
            add_tsb = c(add_tsb, mulaliquots[which(portion_codes == portion_keep)])
            dp_tsb = setdiff(dp_tsb, mulaliquots)
          }else{
            dp_tsb = setdiff(dp_tsb, mulaliquots[which(portion_codes != portion_keep)])
          }
          
        }
      }
      
      if(length(dp_tsb)==0){
        message("ooo Filter barcodes successfully!")
        c(uniq_tsb, add_tsb)
      }else{
        # filter according to plate number
        sampleID_res = substr(dp_tsb, start=1, stop=15)
        dp_samples_res = unique(sampleID_res[duplicated(sampleID_res)])
        for(x in dp_samples_res){
          mulaliquots = dp_tsb[substr(dp_tsb,1,15) == x]
          plate_codes = substr(mulaliquots,
                               start = plate[1],
                               stop = plate[2])
          plate_keep = sort(plate_codes, decreasing = decreasing)[1]
          add_tsb = c(add_tsb, mulaliquots[which(plate_codes == plate_keep)])
          dp_tsb = setdiff(dp_tsb, mulaliquots)
        }
        
        if(length(dp_tsb)==0){
          message("ooo Filter barcodes successfully!")
          c(uniq_tsb, add_tsb)
        }else{
          message("ooo barcodes ", dp_tsb, " failed in filter process, other barcodes will be returned.")
          c(uniq_tsb, add_tsb)
        }
      }
    }
  }
}
dup.filt <- tcgaReplicateFilter(colnames(exp.filt),'RNA')
exp.filt <- exp.filt[,dup.filt]
#去重
exp.filt1 <- aggregate(.~id,exp.filt,mean)
write.table(exp.filt1,'mRNA.unique.txt',quote = F,sep = '\t',row.names = F)
rownames(exp.filt1) <- exp.filt1$id
exp.filt1 <- exp.filt1[,-1]
exp.filt1.log <- apply(exp.filt1,2,function(x){log2(x+1)})

#还得做差异分析看一下，否则无法cMAP####
group_list=c(rep('Normal',99),rep('Tumor',1072))
## 强制限定顺序
group_list <- factor(group_list,levels = c("Normal","Tumor"))
pvalue <-0.05
logFoldChange <- 1
dat <- exp.filt1.log
design = model.matrix(~0+group_list, data=group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(dat)
design
# 构建差异比较矩阵
contrast.matrix <- makeContrasts(Tumor-Normal, levels = design)
fit <- lmFit(dat, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
options(digits = 4)#输出结果的小数点后保留4位
allDiff=topTable(fit2,coef=1,number=Inf)
allDiff <- na.omit(allDiff)
write.table(cbind(Symbol=rownames(allDiff),allDiff),file="mRNA_limmaOut.txt",sep="\t",row.names = F,quote = F)

diffSig = allDiff[(allDiff$P.Value < pvalue & (allDiff$logFC>logFoldChange | allDiff$logFC<(-logFoldChange))),]
write.table(cbind(Symbol=rownames(diffSig),diffSig), file="mRNA_diffSig.txt",sep="\t",row.names = F,quote = F)

diffUp = allDiff[(allDiff$P.Value < pvalue & (allDiff$logFC>logFoldChange)),]
write.table(cbind(Symbol=rownames(diffUp),diffUp), file="mRNA_up.txt",sep="\t",row.names = F,quote = F)

diffDown = allDiff[(allDiff$P.Value < pvalue & (allDiff$logFC<(-logFoldChange))),]
write.table(cbind(Symbol=rownames(diffDown),diffDown), file="mRNA_down.txt",sep="\t",row.names = F,quote = F)
#火山图
pvalue <-0.05
logFC <- 1
allDiff <- read.table('mRNA_limmaOut.txt',header = T,check.names = F,row.names = 1,sep = '\t')
allDiff$Significant <- ifelse(allDiff$P.Value<pvalue & abs(allDiff$logFC)>= logFC,
                              ifelse(allDiff$logFC> logFC,'up','down'),'no')
mycol <- c("#3CB371","#3D3D3D","#FF4500")
pdf(file="DEGs.volcano.pdf",width=6,height=5.5)
p <- ggplot(allDiff, aes(logFC, -log10(P.Value), colour= Significant))+
  geom_point(size=1.2)+theme_bw()+
  scale_color_manual(values = mycol,name='Significant')+
  labs(title="DEGs",x="Log2FC",y="-log10 (P.value)")+
  geom_hline(yintercept = -log10(pvalue),linetype=3,lwd = 1)+
  geom_vline(xintercept = c(-logFC, logFC), linetype=3,lwd = 1)+
  theme(axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15),
        legend.text = element_text(size = 15),text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5))
print(p)
dev.off()

#snoRNA与差异mRNA相关性####
diffexp.t <- data.frame(t(exp.filt1[rownames(diffSig),]),check.names = F)
diffexp.t$sample <- substr(rownames(diffexp.t),1,15)
hubsnoexp <- TCGA_BRCA.sno.unique[c('U81_chr1_173833283_173833360',
                                     'ACA7B;ENSG00000207088;SNORA7B_chr3_129116052_129116191'),]
hubsnoexp <- data.frame(t(hubsnoexp),check.names = F)
hubsnoexp$sample <- rownames(hubsnoexp)
hubsnoexp <- merge(hubsnoexp,diffexp.t,by='sample')
colnames(hubsnoexp)[2:3] <- c('U81','SNORA7B')

#计算
hub <- hubsnoexp[,c(2,3)]
expr <- hubsnoexp[,4:ncol(hubsnoexp)]

#批量计算相关性
corr <- function(gene){
  y <- as.numeric(hub[,gene])
  colnames <- colnames(expr)
  do.call(rbind,lapply(colnames, function(x){
    dd  <- corr.test(as.numeric(expr[,x]), y , method="pearson")
    data.frame(snoRNA=gene,mRNA=x,cor=dd$r,p.value=dd$p )
  }))
}

#以LAMA2为例，测试一下函数
corr("U81")

#批量计算genelist跟免疫浸润相关性的结果
data <- do.call(rbind,lapply(c('U81','SNORA7B'),corr))
head(data)
#保存到文件
write.csv(data, "snoRNA-mRNA.pearson.correlation.csv", quote = F, row.names = F)

#富集分析,分别对U81和SNORA7B的top50mRNA####
#metascape吧，SNORA7B富集不出来
boxdata <- read.delim2('U81.metaEnrichGO.txt',header = T,sep='\t',check.names = F)
boxdata$LogP <- sapply(boxdata$LogP,function(x){as.numeric(x)})
boxdata$GeneRatio <- sapply(boxdata$GeneRatio,function(x){as.numeric(x)})
boxdata <- rbind(subset(boxdata,Category=='GO Biological Processes')[1:10,],
                 subset(boxdata,Category=='GO Molecular Functions')[1:10,],
                 subset(boxdata,Category=='GO Cellular Components')[1:10,]) %>% na.omit(x)
boxdata <- boxdata[order(boxdata$GeneRatio,decreasing = F),]
boxdata$Description <- factor(boxdata$Description,levels = boxdata$Description)
pdf('Fig5C.U81enrichGO.pdf',width = 9.5,height =10.5)
ggplot(boxdata,aes(x=GeneRatio,y=Description,fill = -LogP)) +
  geom_vline(xintercept = c(0.05,0.1),lty=2, lwd=1,col='Azure3')+
  geom_point(shape = 21,aes(size=GeneRatio))+#shape=21是空心圆
  scale_fill_continuous(low="blue",high="red")+
  scale_size_continuous(range = c(8,11))+
  theme_bw()+xlab('GeneRatio')+ylab('')+
  theme(axis.title = element_text(size=18),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+
  scale_y_discrete(labels=function(x) str_wrap(x, width=40))+
  facet_grid(Category~., scale='free')
dev.off()

#KEGG
boxdata <- read.delim2('U81.metaEnrichKEGG.txt',header = T,sep='\t',check.names = F)
boxdata$LogP <- sapply(boxdata$LogP,function(x){as.numeric(x)})
boxdata$Count <- sapply(boxdata$Count,function(x){as.numeric(x)})
boxdata <- rbind(subset(boxdata,Category=='KEGG Pathway')[1:10,],
                 subset(boxdata,Category=='Reactome Gene Sets')[1:11,]) %>% na.omit(x)
pdf('Fig5D.u81enrichKEGG.pdf',width = 10,height = 7)
ggplot(boxdata,aes(x=Count,y=Description)) +
  geom_bar(stat='identity',aes(fill=-(LogP)),colour = 'black')+
  scale_fill_gradient(low="blue",high="red")+
  theme_bw()+xlab('Count')+ylab('')+
  theme(axis.title = element_text(size=18),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+
  facet_grid(Category~., scale='free')+
  scale_y_discrete(labels=function(x) str_wrap(x, width=60))
dev.off()

#overlap <- read.table('corrmRNA.txt',header = T) %>% .$Symbol
#U81
# overlap <- read.table('U81top50.txt',header = T) %>% .$Symbol
# gene <- mapIds(org.Hs.eg.db,keys = overlap,column = 'ENTREZID',keytype = 'SYMBOL',multiVals='filter')
# #save(gene,file='gene.Rdata')
# kk <- enrichGO(gene = gene,
#                OrgDb = org.Hs.eg.db, 
#                pvalueCutoff =0.05, 
#                qvalueCutoff = 0.2,
#                ont="all",
#                readable =T
# )
# k <- enrichKEGG(gene = gene, organism = "hsa",
#                 pvalueCutoff =0.05, qvalueCutoff =0.2)  
# write.table(kk,file="U81enrichGO.txt",quote=F,row.names = F,sep = '\t')               
# write.table(k,file="U81enrichKEGG.txt",quote=F,row.names = F,sep = '\t') 
# 
# pdf(file="Fig5C.U81enrichGO.pdf",width = 10.5,height = 5.5)
# barplot(kk, drop = TRUE, showCategory = 20,split="ONTOLOGY") + 
#   scale_y_discrete(labels=function(x) str_wrap(x, width=60))+
#   facet_grid(ONTOLOGY~., scale='free') 
# dev.off()
# #KEGG八卦图
# kkid_file <- read.table('U81enrichKEGG.txt',header = T,sep='\t',check.names = F)
# diffSig1 <- read.table('U81.corrmRNA.top50.diff.txt',header = T,sep = '\t',check.names = F)
# kkid <- data.frame(Category = 'All',
#                    ID = kkid_file$ID,
#                    Term = kkid_file$Description,
#                    Genes = gsub('/',',',kkid_file$geneID),
#                    adj_pval = kkid_file$p.adjust)
# allentrez <- mapIds(org.Hs.eg.db,keys = diffSig1$Symbol,column = 'ENTREZID',keytype = 'SYMBOL',multiVals='first')
# diffSig1 <- data.frame(ID = allentrez, logFC = diffSig1$logFC) 
# cirkegg <- circle_dat(kkid,diffSig1) 
# pdf('Fig5D.u81enrichKEGG.pdf',width = 10,height = 6)
# GOCircle(cirkegg,rad1=2.5,rad2=3.5,label.size=3.5,nsub=5,zsc.col=c('blue', 'white','yellow'))  #nsub是富集出来的通路数目
# dev.off()


#SNORA7B
#LHFP,10186
#PPAP2B,8613
#SDPR,8436
#C2orf40,84417

boxdata <- read.delim2('SNORA7B.metaEnrichGO.txt',header = T,sep='\t',check.names = F)
boxdata$LogP <- sapply(boxdata$LogP,function(x){as.numeric(x)})
boxdata$GeneRatio <- sapply(boxdata$GeneRatio,function(x){as.numeric(x)})
boxdata <- rbind(subset(boxdata,Category=='GO Biological Processes')[1:10,],
                 subset(boxdata,Category=='GO Molecular Functions')[1:10,],
                 subset(boxdata,Category=='GO Cellular Components')[1:10,]) %>% na.omit(x)
boxdata <- boxdata[order(boxdata$GeneRatio,decreasing = F),]
boxdata$Description <- factor(boxdata$Description,levels = boxdata$Description)
pdf('Fig5E.SNORA7BenrichGO.pdf',width = 9.5,height =11.5)
ggplot(boxdata,aes(x=GeneRatio,y=Description,fill = -LogP)) +
  geom_vline(xintercept = c(0.025,0.05),lty=2, lwd=1,col='Azure3')+
  geom_point(shape = 21,aes(size=GeneRatio))+#shape=21是空心圆
  scale_fill_continuous(low="blue",high="red")+
  scale_size_continuous(range = c(8,11))+
  theme_bw()+xlab('GeneRatio')+ylab('')+
  theme(axis.title = element_text(size=18),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+
  scale_y_discrete(labels=function(x) str_wrap(x, width=40))+
  facet_grid(Category~., scale='free')
dev.off()

#KEGG
boxdata <- read.delim2('SNORA7B.metaEnrichKEGG.txt',header = T,sep='\t',check.names = F)
boxdata$LogP <- sapply(boxdata$LogP,function(x){as.numeric(x)})
boxdata$Count <- sapply(boxdata$Count,function(x){as.numeric(x)})
boxdata <- rbind(subset(boxdata,Category=='KEGG Pathway')[1:10,],
                 subset(boxdata,Category=='Reactome Gene Sets')[1:10,]) %>% na.omit(x)
pdf('Fig5F.SNORA7BenrichKEGG.pdf',width = 10,height = 6)
ggplot(boxdata,aes(x=Count,y=Description)) +
  geom_bar(stat='identity',aes(fill=-(LogP)),colour = 'black')+
  scale_fill_gradient(low="blue",high="red")+
  theme_bw()+xlab('Count')+ylab('')+
  theme(axis.title = element_text(size=18),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        legend.text = element_text(size=16),
        legend.title = element_text(size=16))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+
  facet_grid(Category~., scale='free')+
  scale_y_discrete(labels=function(x) str_wrap(x, width=60))
dev.off()

# overlap <- read.table('SNORA7Btop50.txt',header = T) %>% .$Symbol
# #overlap <- read.table('SNORA7Btop30.txt',header = T) %>% .$Symbol
# gene <- mapIds(org.Hs.eg.db,keys = overlap,column = 'ENTREZID',keytype = 'SYMBOL',multiVals='filter')
# gene[11:14] <- c('10186','8613','8436','84417')
# #save(gene,file='gene.Rdata')
# kk <- enrichGO(gene = gene,
#                OrgDb = org.Hs.eg.db, 
#                pvalueCutoff =0.05, 
#                qvalueCutoff = 0.2,
#                ont="all",
#                readable =T
# )
# k <- enrichKEGG(gene = gene, organism = "hsa",
#                 pvalueCutoff =0.05, qvalueCutoff =0.2)  
# write.table(kk,file="SNORA7BenrichGO.txt",quote=F,row.names = F,sep = '\t')               
# write.table(k,file="SNORA7BenrichKEGG.txt",quote=F,row.names = F,sep = '\t') 
# 
# pdf(file="Fig5E.SNORA7BenrichGO.pdf",width = 10.5,height = 5.5)
# barplot(kk, drop = TRUE, showCategory = 20,split="ONTOLOGY") + 
#   scale_y_discrete(labels=function(x) str_wrap(x, width=60))+
#   facet_grid(ONTOLOGY~., scale='free') 
# dev.off()

#KEGG八卦图
kkid_file <- read.table('SNORA7BenrichKEGG.txt',header = T,sep='\t',check.names = F)
diffSig1 <- read.table('SNORA7B.corrmRNA.top50.diff.txt',header = T,sep = '\t',check.names = F)
kkid <- data.frame(Category = 'All',
                   ID = kkid_file$ID,
                   Term = kkid_file$Description,
                   Genes = gsub('/',',',kkid_file$geneID),
                   adj_pval = kkid_file$p.adjust)
allentrez <- mapIds(org.Hs.eg.db,keys = diffSig1$Symbol,column = 'ENTREZID',keytype = 'SYMBOL',multiVals='first')
diffSig1 <- data.frame(ID = allentrez, logFC = diffSig1$logFC) 
cirkegg <- circle_dat(kkid,diffSig1) 
pdf('Fig5F.SNORA7BenrichKEGG.pdf',width = 10,height = 6)
GOCircle(cirkegg,rad1=2.5,rad2=3.5,label.size=3.5,nsub=5,zsc.col=c('blue', 'white','yellow'))  #nsub是富集出来的通路数目
dev.off()


#mRNA和snoRNA关系图高级版####
library(tidyverse)
library(data.table)
library(ggraph)
library(tidygraph)
source(file = "gather_graph_node.R")
source(file = "gather_graph_edge.R")

#very_easy_input.csv，这里以上面的富集分析结果为例，展示通路和通路里的基因变化倍数FC。
#三列依次是通路-基因-倍数，可以自由替换成“应用场景”里其他需要展示的信息。
#gene_special.txt，要突出显示的基因。第一列是基因名，第二列是类型（例如基因家族信息）。
#我这里的gene_special是两个snoRNA都有关的基因

df <- read.table("snoRNAonly.txt",header = T,sep = '\t',check.names = F)
head(df)
geneSpecial <- read.table("gene_special.txt", header = T)
geneCol <- geneSpecial$Type
names(geneCol) <- geneSpecial$mRNA
geneCol
#使用小丫提供的函数整理数据
#term和gene分别是输入数据snoRNAonly.txt的列名
nodes <- gather_graph_node(df, index = c("snoRNA", "mRNA"), value = "logfc", root="all")
edges <- gather_graph_edge(df, index = c("snoRNA", "mRNA"), root = "all")
nodes <- nodes %>% mutate_at(c("node.level","node.branch"),as.character)
head(nodes, 10)
head(edges, 10)
# 把要突出显示的基因类型信息加到nodes里
nodes$color <- "normal"
nodes[nodes$node.short_name %in% geneSpecial$mRNA,]$color <- geneCol[nodes[nodes$node.short_name %in% geneSpecial$mRNA,]$node.short_name]
nodes[nodes$node.short_name %in% geneSpecial$mRNA,]
nodes$color <- factor(nodes$color, levels = unique(nodes$color))
# 有了节点和边的数据，使用 `tbl_graph()` 便可以得到一个图。
graph <- tbl_graph(nodes, edges)
#开始画图
#自定义配色，直接出图
# 用 `ggraph` 出图
gc <- ggraph(graph, layout = 'dendrogram', circular = TRUE) + 
  # 使用 filter 参数去掉 root（前面设置为"all"）节点及与其相连的边
  geom_edge_diagonal(aes(color = node1.node.branch,
                         filter=node1.node.level!="all"), 
                     alpha = 1/3,edge_width=1) + 
  geom_node_point(aes(size = node.size, 
                      color = node.branch,
                      filter=node.level!="all"), alpha = 1/3) + 
  scale_size(range = c(0.5,3)) + #做均一化处理，让点的大小介于range之间
  theme(legend.position = "none") + #不画图例
  
  # 点和边的配色
  # 如果要改变边的配色，需要同时给边和点变色，否则会对应不上
  scale_edge_color_brewer(palette = "Set1") + #用?scale_color_brewer查看更多配色方案
  scale_color_brewer(palette = "Set1") +
  
  # 添加周围注释文字，此处是基因名gene
  geom_node_text(
    aes(
      x = 1.048 * x, #控制字跟点的距离
      y = 1.048 * y, #控制字跟点的距离
      label = node.short_name,
      angle = -((-node_angle(x, y) + 90) %% 180) + 90,
      filter = leaf,
      color = node.branch
    ),
    size = 6, hjust = 'outward') +
  
  # 添加内环文字，此处是通路名term
  geom_node_text(
    aes(label=node.short_name,
        filter = !leaf & (node.level != "all"),
        color = node.branch),
    fontface="bold",
    size=8,
    family="sans"
  ) + 
  theme(panel.background = element_rect(fill = NA)) +
  coord_cartesian(xlim=c(-1.3,1.3),ylim = c(-1.3,1.3)) #扩大坐标系

gc
ggsave("ccgraph_color.pdf", width = 12, height = 12)



#甲基化数据####
require(Biobase)
library("impute")
info=read.table("sample.txt",sep="\t",header=T)
library(data.table)
b=info
rownames(b)=b[,1]
# 如果你的甲基化信号矩阵，自己在Excel表格里面整理好。
# 就走下面的fread流程
a=fread("methyl\\TCGA.BRCA.sampleMap_HumanMethylation450\\HumanMethylation450",data.table = F )
a[1:4,1:4]
rownames(a)=a[,1]
a=a[,-1]
beta=as.matrix(a)
#服务器跑imput这一步
load('beta2.Rdata')
#beta=impute.knn(beta)
betaData=beta$data
betaData=betaData+0.00001
a=betaData
a[1:4,1:4]
identical(colnames(a),rownames(b))
# 一定要保证，甲基化信号值矩阵，和表型信息，是一一对应的

library(ChAMP)
# beta 信号值矩阵里面不能有NA值
myLoad=champ.filter(beta = a,pd = b)
myLoad
save(myLoad,file = 'step1-output.Rdata')

#step2-check-betaM
library("minfi")
if(F){
  myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=5)
  dim(myNorm) 
  pD=myLoad$pd
  save(myNorm,pD,file = 'step2-champ_myNorm.Rdata')
}
load(file = 'step2-champ_myNorm.Rdata')
# 原来的450K经过质控过滤后是400K啦
#beta.m=myNorm
group_list=myLoad$pd$Group
dim(myNorm) 

#查找U81和SNORA7B的甲基化位点
#SNORA7B:ENSG00000207088,位置chr3:129397210-129397348
#U81:U81_chr1_173833283_173833360
anno <- read.table('probeMap_illuminaMethyl450_hg19_GPL16304_TCGAlegacy',sep = '\t',
                   check.names = F,header = F)

myNorm$V1 <- rownames(myNorm)



#高低风险组的免疫####
#ssgsea+ESTIMATE####
clinical <- cbind(rbind(train[,1:10],vali[,1:10]),
                        rbind(train[,modelgene],vali[,modelgene]))
clinical <- clinical[order(clinical$risk,decreasing = T),]
clinical <- clinical[,c(1,2)]
exp.filt1.tumor.t <- data.frame(t(exp.filt1[,100:ncol(exp.filt1)]),check.names = F)
exp.filt1.tumor.t$Patient <- substr(rownames(exp.filt1.tumor.t),1,12)
ES.exp <- merge(clinical,exp.filt1.tumor.t,by='Patient')
ES.exp <- ES.exp[order(ES.exp$risk,decreasing = T),]
rownames(ES.exp) <- ES.exp$Patient
ES.exp <- ES.exp[,-c(1,2)]
#432低风险，430高风险
ES.exp <- data.frame(t(ES.exp),check.names = F)
write.table(cbind(Symbol=rownames(ES.exp),ES.exp),'Estimate.txt',sep = '\t',quote = F,row.names = F)
filterCommonGenes(input.f='Estimate.txt',
                  output.f="commonGenes.gct",
                  id="GeneSymbol")

estimateScore(input.ds = "commonGenes.gct",
              output.ds="estimateScore.gct",
              platform="illumina")
#输出样本打分表
scores=read.table("estimateScore.gct",skip = 2,header = T)
rownames(scores)=scores[,1]
scores=t(scores)
scores=scores[-1,]
rownames(scores)=gsub("\\.","\\-",rownames(scores))
#out=rbind(ID=colnames(scores),scores)
write.table(scores,file="estimate_scores.txt",sep="\t",quote=F,col.names=F)

#作图
data <- scores[-1,]
data <- apply(data,2,function(x){as.numeric(x)}) 
data <-  data.frame(data,check.names = F)
rownames(data) <- rownames(scores[-1,])
data$Risk <- factor(c(rep('Low_risk',432),rep('High_risk',430)),levels = c('Low_risk','High_risk'))
p1 <- ggviolin(data,x='Risk',y='ImmuneScore',fill ='Risk',
               ylab = 'Fraction',xlab ='Immunocyte',
               palette = "jco",trim = F,
               position=position_dodge(0.9),
               add ='boxplot',add.params = list(fill = "white"))+guides(fill=FALSE)+
  theme(axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title = element_text(size=15),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15),
        plot.title = element_text(size=12))+
  stat_compare_means(size = 6,aes(group=Risk),
                     label = 'p.signif',label.y = 4200,label.x = 1.4)

p2 <- ggviolin(data,x='Risk',y='ESTIMATEScore',fill ='Risk',
               ylab = 'Fraction',xlab ='ESTIMATEScore',
               palette = "jco",trim = F,
               position=position_dodge(0.9),
               add ='boxplot',add.params = list(fill = "white"))+guides(fill=FALSE)+
  theme(axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title = element_text(size=15),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15),
        plot.title = element_text(size=12))+
  stat_compare_means(size = 6,aes(group=Risk),
                     label = 'p.signif',label.y = 5500,label.x = 1.4)

p3 <- ggviolin(data,x='Risk',y='StromalScore',fill ='Risk',
               ylab = 'Fraction',xlab ='StromalScore',
               palette = "jco",trim = F,
               position=position_dodge(0.9),
               add ='boxplot',add.params = list(fill = "white"))+guides(fill=FALSE)+
  theme(axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title = element_text(size=15),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15),
        plot.title = element_text(size=12))+
  stat_compare_means(size = 6,aes(group=Risk),
                     label = 'p.signif',label.y =2200,label.x = 1.4)
pdf('ESTIMATE.pdf',width = 11,height = 5)
p1+p2+p3
dev.off()

#ssGSEA####
inputFile <- 'Estimate.txt'
rt=read.table(inputFile,sep="\t",header=T,check.names=F,row.names = 1)
exp <- as.matrix(rt)
dimnames=list(rownames(exp),colnames(exp))
mat=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
mat=mat[rowMeans(mat)>0,]
# gmtFile <- 'immunity_by_CELL.gmt'
# immunity=getGmt(gmtFile)
load('easy_input_immunity.symbol.rdata')
#ssGSEA
ssgseaScore=gsva(mat, immunity, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
#输出
ssgseaOut=normalize(ssgseaScore)
ssgseaOut=rbind(id=colnames(ssgseaOut),ssgseaOut)
write.table(ssgseaOut,file="immune_ssgseaOut.txt",sep="\t",quote=F,col.names=F)

#检测哪些细胞是差异的
data <- ssgseaOut[-1,]
data <- apply(data,2,function(x){as.numeric(x)})
rownames(data) <- rownames(ssgseaOut[-1,])
data <- data.frame(t(data),check.names = F)
data$Type <- factor(c(rep('Low_risk',432),rep('High_risk',430)),levels = c('Low_risk','High_risk'))
sig <- data.frame()
for (i in colnames(data[,1:(ncol(data)-1)])) {
  rt <- data.frame(immune=data[,i],Type=data$Type)
  rtTest<- wilcox.test(immune~Type, data=rt)
  if (rtTest$p.value <0.05){
    sig <- rbind(sig,data.frame(Immunocyte = i,pvalue = rtTest$p.value))
  }
}

#箱线图
#画图
dbox <- ssgseaOut[-1,]
dbox <- apply(dbox,2,function(x){as.numeric(x)})
rownames(dbox) <- rownames(ssgseaOut)[-1]
dbox <- data.frame(t(dbox),check.names = F)
dbox$Risk <- factor(c(rep('Low_risk',432),rep('High_risk',430)),levels = c('Low_risk','High_risk'))
dbox  <- gather(dbox,Immunocyte,Score,1:(ncol(dbox)-1))
pdf('Fig5A.ssGSEA.pdf',width = 8,height = 5)
ggboxplot(dbox,x='Immunocyte',y='Score',fill ='Risk',
               ylab = 'Score',
               xlab ='Immunocyte',palette = "jco",
               #add ='jitter',add.params = list(size=0.4),
               size =0.6,outlier.shape = NA)+
  rotate_x_text(45)+
  theme(axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=13),
        axis.title = element_text(size=15),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15),
        axis.line = element_line(size=0.8))+
  stat_compare_means(size = 4,aes(group=Risk),
                     label = 'p.signif')
dev.off()

#IPS####
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
#自定义函数，用于计算 Immunophenoscore(IPS)
ipsmap<- function (x) {
  if (x<=0) {
    ips<-0
  } else {
    if (x>=3) {
      ips<-10
    } else {
      ips<-round(x*10/3, digits=0)
    }
  }
  return(ips)
}
#分配颜色
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 1000)
mapcolors<-function (x) {
  za<-NULL
  if (x>=3) {
    za=1000
  } else {
    if (x<=-3) {
      za=1
    } else {
      za=round(166.5*x+500.5,digits=0)
    }
  }
  return(my_palette[za])
}
my_palette2 <- colorRampPalette(c("black", "white"))(n = 1000)
mapbw<-function (x) {
  za2<-NULL
  if (x>=2) {
    za2=1000
  } else {
    if (x<=-2) {
      za2=1
    } else {
      za2=round(249.75*x+500.5,digits=0)
    }
  }
  return(my_palette2[za2])
}
#输入文件
#easy_input_expr.csv，表达矩阵。例如log2(TPM+1)，每行一个基因，每列一个样本。
#easy_input_IPS_genes.txt，IPS基因权重文件。IPS相关基因、分类及其权重。
## 读取肿瘤样本的表达谱文件
# expr <- read.table('Estimate.txt',sep="\t",header=T,check.names=F,row.names = 1)

## 修改几个基因名
## IPS基因权重文件中的基因名跟常用的gene symbol不同，我们把表达矩阵里的gene symbol修改成跟IPS基因列表一致的基因名。
if(is.element("CAVIN2",rownames(expr))) {
  rownames(expr) <- gsub("CAVIN2","SDPR",rownames(expr)) # 根据https://www.genecards.org上的相关信息进行同名转换
}
# 另外注意示例数据的里的CCL3L3的注释源自GENCODE27，GeneCards显示CCL3L1的ENSEMBL为ENSG00000276085，该编码在GENCODE27里注释为CCL3L3，因此修改
if(is.element("CCL3L3",rownames(expr))) {
  rownames(expr) <- gsub("CCL3L3","CCL3L1",rownames(expr))
}

## 读取IPS相关基因以及其权重
IPSG <- read.table("easy_input_IPS_genes.txt",header=TRUE, sep="\t",check.names=FALSE,stringsAsFactors = FALSE,row.names = NULL)
head(IPSG)
unique_ips_genes <- as.vector(unique(IPSG$NAME))

## 初始化数据
IPS <- MHC <- CP <- EC <- SC <- AZ <- NULL

# 获取表达谱里的基因名
GVEC <- row.names(expr)

# 获取IPS基因文件里的基因名
VEC <- as.vector(IPSG$GENE)

# 匹配基因并找到缺失基因
ind <- which(is.na(match(VEC,GVEC)))
MISSING_GENES <- VEC[ind]

dat <- IPSG[ind,]
if (length(MISSING_GENES) > 0) { # 若存在缺失基因报出（如果上面两个基因不修改，这里是要报缺失的，会导致最终结果错误）
  message(paste0("--differently named or missing genes: ",paste(MISSING_GENES,collapse = ", ")))
  print(IPSG[ind,])
  message("please check data and make sure all genes matches!")
} else {
  message("--all genes matched!") # 请确保所有基因都匹配!!!
}
#画图
## 循环每个样本
# sample_names <- colnames(expr)
# tumsam  <- sample_names
#输出每个样本的IPS计算的3个免疫检查位点的得分
mydata <- data.frame()
for (i in 1:length(sample_names)) { 
  GE <- expr[[i]]
  mGE <- mean(GE)
  sGE <- sd(GE)
  Z1 <- (expr[as.vector(IPSG$GENE),i]-mGE)/sGE
  W1 <- IPSG$WEIGHT
  WEIGHT <- MIG <- NULL
  k <- 1
  for (gen in unique_ips_genes) {
    MIG[k] <- mean(Z1[which (as.vector(IPSG$NAME)==gen)],na.rm=TRUE)
    WEIGHT[k] <- mean(W1[which (as.vector(IPSG$NAME)==gen)])
    k <- k + 1
  }
  WG <- MIG * WEIGHT
  MHC[i] <- mean(WG[1:10])
  CP[i] <- mean(WG[11:20])
  EC[i] <- mean(WG[21:24])
  SC[i] <- mean(WG[25:26])
  AZ[i] <- sum(MHC[i],CP[i],EC[i],SC[i])
  IPS[i] <- ipsmap(AZ[i])
  
  # 绘制Immunophenogram
  data_a <- data.frame (start = c(0,2.5,5,7.5,10,15,seq(20,39),0,10,20,30), end = c(2.5,5,7.5,10,15,seq(20,40),10,20,30,40), y1=c(rep(2.6,26),rep(0.4,4)),y2=c(rep(5.6,26),rep(2.2,4)),z=c(MIG[c(21:26,11:20,1:10)],EC[i],SC[i],CP[i],MHC[i]),vcol=c(unlist(lapply(MIG[c(21:26,11:20,1:10)],mapcolors)), unlist(lapply(c(EC[i],SC[i],CP[i],MHC[i]),mapbw))), label = c(unique_ips_genes[c(21:26,11:20,1:10)],"EC","SC","CP","MHC"))
  data_a$label <- factor(data_a$label, levels=unique(data_a$label))
  mydata1 <- data_a[data_a$label=='PD-L1',] %>% .[,c('z','label')]
  mydata2 <- data_a[data_a$label=='PD-1',] %>% .[,c('z','label')]
  mydata3 <- data_a[data_a$label=='CTLA-4',] %>% .[,c('z','label')]
  mydata4 <- rbind(mydata1,mydata2,mydata3)
  mydata4$sample <- colnames(expr)[i]
  mydata4 <- spread(mydata4,label,z)
  mydata <- rbind(mydata,mydata4)
  # plot_a1 <- ggplot() + 
  #   geom_rect(data=data_a, mapping=aes(xmin=start, xmax=end, ymin=y1, ymax=y2, fill=label), size=0.5,color="black", alpha=1) +
  #   coord_polar()  + 
  #   scale_y_continuous(limits = c(0, 6)) +
  #   scale_fill_manual(values =as.vector(data_a$vcol),guide=FALSE) +
  #   theme_bw() + 
  #   theme(panel.spacing = unit(0, 'mm'), 
  #         panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank(),
  #         panel.border = element_blank(),
  #         panel.background = element_blank(),
  #         axis.line = element_line(colour = "white"),
  #         axis.text=element_blank(), 
  #         axis.ticks= element_blank()) + 
  #   geom_text(aes(x=5, y=1.3, label="EC"), size=4) +
  #   geom_text(aes(x=15, y=1.3, label="SC"), size=4) + 
  #   geom_text(aes(x=25, y=1.3, label="CP"), size=4) + 
  #   geom_text(aes(x=35, y=1.3, label="MHC"), size=4)
  # plot_a2 <- plot_a1 +
  #   geom_text(aes(x=1.25, y=4.1, label="+ Act CD4"), angle=78.75, size=4) +
  #   geom_text(aes(x=3.75, y=4.1, label="+ Act CD8"),angle=56.25, size=4) +
  #   geom_text(aes(x=6.25, y=4.1, label="+ Tem CD4"), angle=33.75,size=4) +
  #   geom_text(aes(x=8.75, y=4.1, label="+ Tem CD8"), angle=11.25,size=4) +
  #   geom_text(aes(x=17.5, y=4.1, label="- MDSC"), angle=-67.5,size=4) +
  #   geom_text(aes(x=12.5, y=4.1, label="- Treg"), angle=-22.5,size=4)
  # plot_a3 <- plot_a2 +
  #   geom_text(aes(x=20.5, y=4.1, label="PD-1 -"), angle=85.5, size=4) +
  #   geom_text(aes(x=21.5, y=4.1, label="CTLA4 -"), angle=76.5, size=4) +
  #   geom_text(aes(x=22.5, y=4.1, label="LAG3 -"), angle=67.5, size=4) +
  #   geom_text(aes(x=23.5, y=4.1, label="TIGIT -"), angle=58.5, size=4) +
  #   geom_text(aes(x=24.5, y=4.1, label="TIM3 -"), angle=49.5, size=4) +
  #   geom_text(aes(x=25.5, y=4.1, label="PD-L1 -"), angle=40.5, size=4) +
  #   geom_text(aes(x=26.5, y=4.1, label="PD-L2 -"), angle=31.5, size=4) +
  #   geom_text(aes(x=27.5, y=4.1, label="CD27 +"), angle=22.5, size=4) +
  #   geom_text(aes(x=28.5, y=4.1, label="ICOS +"), angle=13.5, size=4) +
  #   geom_text(aes(x=29.5, y=4.1, label="IDO1 -"), angle=4.5, size=4)
  # plot_a4 <- plot_a3 +
  #   geom_text(aes(x=30.5, y=4.1, label="B2M +"), angle=-4.5, size=4) +
  #   geom_text(aes(x=31.5, y=4.1, label="TAP1 +"), angle=-13.5, size=4) +
  #   geom_text(aes(x=32.5, y=4.1, label="TAP2 +"), angle=-22.5, size=4) +
  #   geom_text(aes(x=33.5, y=4.1, label="HLA-A +"), angle=-31.5, size=4) +
  #   geom_text(aes(x=34.5, y=4.1, label="HLA-B +"), angle=-40.5, size=4) +
  #   geom_text(aes(x=35.5, y=4.1, label="HLA-C +"), angle=-49.5, size=4) +
  #   geom_text(aes(x=36.5, y=4.1, label="HLA-DPA1 +"), angle=-58.5, size=4) +
  #   geom_text(aes(x=37.5, y=4.1, label="HLA-DPB1 +"), angle=-67.5, size=4) +
  #   geom_text(aes(x=38.5, y=4.1, label="HLA-E +"), angle=-76.5, size=4) +
  #   geom_text(aes(x=39.5, y=4.1, label="HLA-F +"), angle=-85.5, size=4)
  # plot_a5 <- plot_a4 +
  #   geom_text(aes(x=0, y=6, label=paste("Immunophenoscore: ",IPS[i],sep="")), angle=0,size=6,vjust=-0.5) + 
  #   theme(axis.title=element_blank())
  # plot_a <- plot_a5 + theme(plot.margin=unit(c(0,0,0,0),"mm")) +
  #   geom_text(vjust=1.15,hjust=0,aes(x=25.5, y=6,label="\n\n\n\n   MHC: Antigen Processing                                 EC: Effector Cells\n   CP: Checkpoints | Immunomodulators              SC: Suppressor Cells\n\n", hjust = 0), size=4)
  # 
  # ## 绘制图例 sample-wise (averaged) z-scores
  # data_b <- data.frame (start = rep(0,23), 
  #                       end = rep(0.7,23), 
  #                       y1=seq(0,22,by=1), 
  #                       y2=seq(1,23,by=1),
  #                       z=seq(-3,3,by=6/22),
  #                       vcol=c(unlist(lapply(seq(-3,3,by=6/22),mapcolors))), 
  #                       label = LETTERS[1:23])
  # data_b_ticks <- data.frame(x = rep(1.2, 7), 
  #                            value = seq(-3,3, by=1), 
  #                            y = seq(0,6, by=1)*(22/6) +0.5)
  # legendtheme <- theme(plot.margin = unit(c(2,0,2,0),"inch"), 
  #                      panel.spacing = unit(0,"null"), 
  #                      panel.grid.major = element_blank(),
  #                      panel.grid.minor = element_blank(),
  #                      panel.border = element_blank(),
  #                      panel.background = element_blank(), 
  #                      axis.line = element_line(colour = "white"), 
  #                      axis.text=element_blank(), 
  #                      axis.ticks= element_blank(), 
  #                      axis.title.x=element_blank())
  # plot_b <- ggplot(hjust=0) + 
  #   geom_rect(data=data_b, 
  #             mapping=aes(xmin=start, xmax=end, ymin=y1, ymax=y2, fill=label), size=0.5,color="black", alpha=1) +
  #   scale_x_continuous(limits = c(0, 1.5),expand = c(0,0)) +
  #   scale_fill_manual(values =as.vector(data_b$vcol),guide=FALSE) +
  #   geom_text(data=data_b_ticks, aes(x=x, y=y, label=value),hjust="inward", size=4) +
  #   theme_bw() +
  #   legendtheme + 
  #   ylab("Sample-wise (averaged) z-score")
  # 
  # ## 绘制图例 weighted z-scores
  # data_c <- data.frame(start = rep(0,23),
  #                      end = rep(0.7,23),
  #                      y1=seq(0,22,by=1),
  #                      y2=seq(1,23,by=1),
  #                      z=seq(-2,2,by=4/22),
  #                      vcol=c(unlist(lapply(seq(-2,2,by=4/22),mapbw))), 
  #                      label = LETTERS[1:23])
  # data_c_ticks <- data.frame(x = rep(1.2, 5), 
  #                            value = seq(-2,2, by=1),
  #                            y = seq(0,4, by=1)*(22/4) +0.5)
  # plot_c <- ggplot() + 
  #   geom_rect(data=data_c, 
  #             mapping=aes(xmin=start, xmax=end, ymin=y1, ymax=y2, fill=label),size=0.5,color="black", alpha=1) + 
  #   scale_x_continuous(limits = c(0, 1.5),expand = c(0,0)) + 
  #   scale_fill_manual(values =as.vector(data_c$vcol),guide=FALSE) + 
  #   geom_text(data=data_c_ticks, aes(x=x, y=y, label=value),hjust="inward", size=4) +
  #   theme_bw() + 
  #   legendtheme + 
  #   ylab("Weighted z-score")
  # 
  # ## 每一个样本保存图片
  # file_name<-paste("IPS/IPS_",sample_names[i],".pdf",sep="")
  # pdf(file_name, width=10, height=8)
  # grid.arrange(plot_a,plot_b,plot_c, ncol=3, widths=c(0.8,0.1,0.1))
  # invisible(dev.off())
  # 
  # message(paste0("--analysis of ", tumsam[i]," done..."))
}
write.table(mydata,'IPS_checkpoint_score.txt',sep = '\t',row.names = F,quote = F)
# 构建结果，包括各项得分以及IPS
DF <- data.frame(SAMPLE=sample_names,
                 MHC=MHC,
                 EC=EC,
                 SC=SC,
                 CP=CP,
                 AZ=AZ,
                 IPS=IPS,
                 stringsAsFactors = F)

# 输出IPS到文件
write.table(DF,file = "output_IPS.txt", row.names = FALSE, col.names = TRUE, quote=FALSE, sep="\t") 

#IPS结果可视化####
rownames(DF) <- DF$SAMPLE
rownames(DF) <- as.character(sapply(DF$SAMPLE, function(x) substr(x, 1, 12)))
# DF <- DF[colnames(ES.exp),]
# DF$Risk <- factor(c(rep('Low_risk',432),rep('High_risk',430)),levels = c('Low_risk','High_risk'))
DF <- merge(cox.data.plot, DF, by = 0)
p1 <- ggboxplot(DF,x='risk',y='MHC',color ='risk',
               ylab = 'Antigen Processing',
               palette = "jco",
               add = 'jitter')+
  theme(axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title = element_text(size=15),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15),
        plot.title = element_text(size=12))+
  stat_compare_means(size = 6,aes(group=risk),
                     label = 'p.signif',label.x = 1.4)

p2 <- ggboxplot(DF,x='risk',y='EC',color ='risk',
               ylab = 'Effector Cells',
               palette = "jco",
               add = 'jitter')+
  theme(axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title = element_text(size=15),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15),
        plot.title = element_text(size=12))+
  stat_compare_means(size = 6,aes(group=risk),
                     label = 'p.signif',label.x = 1.4)

p3 <- ggboxplot(DF,x='risk',y='SC',color ='risk',
               ylab = 'Suppressor Cells',
               palette = "jco",
               add = 'jitter')+
  theme(axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title = element_text(size=15),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15),
        plot.title = element_text(size=12))+
  stat_compare_means(size = 6,aes(group=risk),
                     label = 'p.signif',label.x = 1.4)
p4 <- ggboxplot(DF,x='risk',y='CP',color ='risk',
               ylab = 'Checkpoints | Immunomodulators',
               palette = "jco",
               add = 'jitter')+ylim(-0.2,0.2)+
  theme(axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title = element_text(size=15),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15),
        plot.title = element_text(size=12))+
  stat_compare_means(size = 6,aes(group=risk),
                     label = 'p.signif',label.x = 1.4)
p5 <- ggviolin(DF,x='risk',y='IPS',fill ='risk',
               ylab = 'Fraction',xlab ='Immunophenoscore',
               palette = "jco",trim = F,
               position=position_dodge(0.9),
               add ='boxplot',add.params = list(fill = "white"))+guides(fill=FALSE)+
  theme(axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title = element_text(size=15),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15),
        plot.title = element_text(size=12))+
  stat_compare_means(size = 6,aes(group=risk),
                     label = 'p.signif',label.x = 1.4)
pdf('Fig5B.IPS.pdf',width = 11,height = 5)
p1+p2+p3
dev.off()

#三个checkpoint得分差异
rownames(mydata) <- mydata$sample
mydata <- mydata[colnames(ES.exp),]
mydata$Risk <- factor(c(rep('Low_risk',432),rep('High_risk',430)),levels = c('Low_risk','High_risk'))
p1 <- ggboxplot(mydata,x='Risk',y='PD-1',color ='Risk',
                ylab = 'IPS-PD-1',
                palette = "jco",
                add = 'jitter')+
  theme(axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title = element_text(size=15),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15),
        plot.title = element_text(size=12))+
  stat_compare_means(size = 6,aes(group=Risk),
                     label = 'p.signif',label.x = 1.5,label.y = 0)

p2 <- ggboxplot(mydata,x='Risk',y='PD-L1',color ='Risk',
                ylab = 'IPS-PD-L1',
                palette = "jco",
                add = 'jitter')+
  theme(axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title = element_text(size=15),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15),
        plot.title = element_text(size=12))+
  stat_compare_means(size = 6,aes(group=Risk),
                     label = 'p.signif',label.x = 1.5,label.y = 0.05)

p3 <- ggboxplot(mydata,x='Risk',y='CTLA-4',color ='Risk',
                ylab = 'IPS-CTLA-4',
                palette = "jco",
                add = 'jitter')+
  theme(axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title = element_text(size=15),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15),
        plot.title = element_text(size=12))+
  stat_compare_means(size = 6,aes(group=Risk),
                     label = 'p.signif',label.x = 1.5,label.y = 0)
pdf('Fig5B.IPS.checkpoint.pdf',width = 12,height = 4)
p1+p2+p3
dev.off()

#TIDE####
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
#表达量标准化
TIDE <- expr
# 这里为了得到比较好的结果，采用two direction median centered
TIDE <- sweep(TIDE,2, apply(TIDE,2,median,na.rm=T))
TIDE <- sweep(TIDE,1, apply(TIDE,1,median,na.rm=T))
write.table(TIDE,"TIDE_input.self_subtract.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

#------------------------------------#
# 请参照文件夹中的TIDE使用教程.docx完成该部分#
#------------------------------------#

##TIDE结果可视化
TIDE.res <- read.csv("TIDE_output.csv",header = T,row.names = 1,check.names = F,stringsAsFactors = F)
TIDE.res <- TIDE.res[colnames(ES.exp),]
TIDE.res$Risk <- factor(c(rep('Low_risk',432),rep('High_risk',430)),levels = c('Low_risk','High_risk'))
boxdata <- TIDE.res
pdf('Fig5C.TIDE.boxplot.pdf',width = 6,height = 6)
ggplot(boxdata,aes(x=Risk, y=TIDE, fill=Risk))+
  geom_violin(width=1.4,trim = FALSE,color='white')+
  geom_boxplot(width=0.5,position=position_dodge(0.9),fill='white')+
  stat_compare_means(size = 8,label = "p.signif",label.x = 1.5,label.y = 1.8)+
  theme_bw()+ggtitle('TIDE prediction')+
  ylab('')+xlab('')+scale_fill_jco()+
  theme(axis.text.y=element_text(size=13),
        axis.text.x=element_text(size=13),
        axis.title = element_text(size=13),
        legend.text = element_text(size=10),
        legend.title = element_text(size=10),
        plot.title = element_text(size=17,hjust = 0.5))
dev.off()

# 参照文件夹中TIDE使用教程得到输出文件TIDE_output.csv
ann <- data.frame(Risk=TIDE.res$Risk,TIDE=TIDE.res$Responder)
rownames(ann) <- rownames(TIDE.res)
print(table(ann$TIDE,ann$Risk))
chisq <- chisq.test(table(ann$TIDE,ann$Risk))

#response和非response的survival
clinical <- read.table('OSinput.txt',header = T,sep = '\t',check.names = F) %>% .[,c(1,4,9)]
TIDE.res$Patient <- rownames(TIDE.res)
TIDE.res.clinical <- merge(TIDE.res,clinical,by='Patient')
fit <- survfit(Surv(futime, fustate) ~ Responder, data = TIDE.res.clinical)
###summary(fit)
# text.legend1 <- c(paste0("High_risk: n = ",fit$n[1]),paste0("Low_risk: n = ",fit$n[2]))
# text.legend2 <-paste0(" PValue : ",pValue)
p <- ggsurvplot(fit, data = TIDE.res.clinical ,
                ggtheme = theme_bw(), #想要网格就运行这行
                conf.int = F, #不画置信区间，想画置信区间就把F改成T
                #conf.int.style = "step",#置信区间的类型，还可改为ribbon
                risk.table = TRUE, # 绘制累计风险曲线
                #tables.theme = theme_void(),
                tables.height = 0.25,
                censor = T, #不显示观察值所在的位置
                palette = c("DodgerBlue4","Red3"), #线的颜色对应高、低(自定义调色板)
                xlab = "Time (Days)",
                legend.title = "",#基因名写在图例题目的位置
                font.legend = 15,#图例的字体大小
                #font.title = 12,font.x = 10,font.y = 10,#设置其他字体大小
                legend=c(0.8,0.88),
                legend.labs = c('False','True'),pval = T)


#snoRNA与全部mRNA相关性####
exp.t <- data.frame(t(exp.filt1),check.names = F)
exp.t$sample <- substr(rownames(exp.t),1,15)
hubsnoexp <- TCGA_BRCA.sno.unique[c('U81_chr1_173833283_173833360',
                                    'ACA7B;ENSG00000207088;SNORA7B_chr3_129116052_129116191'),]
hubsnoexp <- data.frame(t(hubsnoexp),check.names = F)
hubsnoexp$sample <- rownames(hubsnoexp)
hubsnoexp <- merge(hubsnoexp,exp.t,by='sample')
colnames(hubsnoexp)[2:3] <- c('U81','SNORA7B')

#计算
hub <- hubsnoexp[,c(2,3)]
expr <- hubsnoexp[,4:ncol(hubsnoexp)]

#批量计算相关性
corr <- function(gene){
  y <- as.numeric(hub[,gene])
  colnames <- colnames(expr)
  do.call(rbind,lapply(colnames, function(x){
    dd  <- corr.test(as.numeric(expr[,x]), y , method="pearson")
    data.frame(snoRNA=gene,mRNA=x,cor=dd$r,p.value=dd$p )
  }))
}

#以LAMA2为例，测试一下函数
corr("U81")

#批量计算genelist跟免疫浸润相关性的结果
data <- do.call(rbind,lapply(c('U81','SNORA7B'),corr))
head(data)
#保存到文件
write.csv(data, "snoRNA-allmRNA.correlation.csv", quote = F, row.names = F)
