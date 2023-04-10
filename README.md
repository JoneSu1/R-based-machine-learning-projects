# R-based-machine-learning-projects
Included analyses: Analysis of variance, WGCNA analysis, GO/KEGG enrichment analysis, GSVA analysis, correlation expression analysis, ROC diagnostic curve, Lasso regression analysis, XGbost analysis, random forest, single-factor COX, multi-factor COX, random vector machine SVM analysis, K-M survival analysis,

关于代码的可用性，在Code文件中我已经做了足够的注释。这个Read文件，我将着重于对项目的思路和项目的结果做一个介绍.

这个项目是为了探究和OVOL1基因相关连的基因，并研究和这些基因有关的生物学功能，以及可能存在的免疫治疗靶点.

1 数据集
从GEO数据库（https://www.ncbi.nlm.nih.gov/geo/）下载恶性皮肤黑色素瘤GSE46517数据集，共有195个样本。73例为转移性黑色素瘤，31例为原发性黑色素瘤，15例为正常样本.
目标基因: OVOL1.
2 分析方法
1)	从GEO数据库（https://www.ncbi.nlm.nih.gov/geo/）下载恶性皮肤黑色素瘤GSE46517数据集，共有195个样本。73例为转移性黑色素瘤，31例为原发性黑色素瘤，15例为正常样本.

2)	差异分析：首先以目标基因OVOL1的表达量中位数进行分组（高表达，低表达），然后使用R包“limma”对数据集进行差异分析，其中差异基因的筛选阈值为|logFC|>0.5&P.adjust <0.05。

3)	功能富集：使用R包“clusterProfiler”对上述得到的差异基因进行GO以及KEGG功能富集。 
4)	单基因GSEA功能富集：使用R包“clusterProfiler”对数据集进行OVOL1的单基因GSEA分析。 

5)	免疫浸润分析：使用ssGSEA计算得到TCGA-GBM中每个样本的免疫细胞的含量，分析目标基因高低表达组之间的差异免疫细胞。对差异的免疫细胞，进一步计算和CARD16的相关性。
6)	免疫检查点相关性分析：计算CARD16和免疫检查点的相关性免疫治疗分析：使用线上数据库TIDE（https://cloud.genepattern.org/gp/pages/login.jsf）检验免疫治疗响应性和亚型是否相关，并在R中可视化结果.
7)	免疫治疗分析：使用线上数据库GenePattern(https://cloud.genepattern.org/gp) 预测亚型的免疫治疗响应性.
3 分析结果
3.1差异分析
先从GEO数据库（https://www.ncbi.nlm.nih.gov/geo/）下载恶性皮肤黑色素瘤GSE46517数据集，共有195个样本。73例为转移性黑色素瘤，31例为原发性黑色素瘤，15例为正常样本.
然后以目标基因OVOL1的表达量中位数（5.974117）进行分组（高表达，低表达），然后使用R包“limma”对数据集进行差异分析，其中差异基因的筛选阈值为|logFC|>0.5&P <0.05，
筛选出1082个差异基因，其中有520个基因下调，562个基因上调。 为了从整体上了解差异基因的分布，用火山图来展示其分布如图1-1 所示； 
  
图1. 差异火山图
注：（横坐标差异倍数(High/low，取2的对数)， 纵坐标表示-log10(p)。 图中每个点代表一个基因，蓝色和红色的点代表显著差异表达基因。 
且红色的点表示其基因表达量是上调的（High-expression样本相对于Low-expression样本），蓝色的点表示其基因表达量是下调的（High-expression样本相对于 Low-expression样本），
灰色的点表示这些基因没有显著差异） 

热图是对实验数据分布情况进行分析的直观可视化方法，可以用来进行实验数据的质量控制和直观展示重点研究对象的表达量数据差异变化情况，还可
以对数据和样品进行聚类，观测样品质量。它有多种形式，但基本的元素却通用的。如图 1-2 所示。
 
图 1-2： 样本组间差异基因的热图
（每个小方格表示每个基因，其颜色表示该基因表达量大小，表达量越大颜色越深（红色为高表达， 蓝色为低表达）。第一行表示样本分组，蓝色表示Low-expression样本，
红色表示 High-expression样本。每行表示每个基因不同样本中的表达量情况， 每列表示每个样品中所有差异基因的表达量情况。 左侧树状图表示对来自不同样本的不同基因的聚类分析结果）



3.2.2、功能富集 
使用上述得到的1056个差异基因进行基于 GO 、 KEGG 通路的富集分析，寻找基因集合内大量基因共同的功能及相关通路。使用统计学方法累计超几何分布分析一组基因在某个功能结点上是否过出现（over-presentation），其计算公式如下： 

其中 N 为注释系统中基因的总数， n 为要考察的结点或通路本身所注释的基因数， M 为差异表达基因集大小， x 为基因集与结点或通路的交集数目。京都基因与基因组百科全书（KEGG）是了解高级功能和生物系统（如细胞、生物和生态系统），
从分子水平信息，尤其是大型分子数据集生成的基因组测序和其他高通量实验技术的实用程序数据库资源，由日本京都大学生物信息学中心的Kanehisa 实验室于 1995 年建立，是国际最常用的生物信息数据库之一，
以“理解生物系统的高级功能和实用程序资源库”著称。基因本体论(Gene ontology，GO)系统包括三个部分：生物学过程（biological process， BP）、分子功能（molecularfunctions， MF）、细胞组分（cellular components， CC）。
我们使用 R 包“clusterProfiler”进行 GO 和 KEGG 功能富集分析，寻找差异表达基因集合内大量基因共同的功能及相关通路。 
在富集分析的过程中，我们设定 GO、 KEGG 通路富集分析的阈值为：P.adjust <0.05。其中共显著富集到GO功能条目1263条,其中BP有1042条，其中CC有101条，MF有120条（附件：GO-enrich.csv），
并绘制气泡图（图2-1）,其中TOP5的BP(epidermis development，skin development，keratinocyte differentiation，epidermal cell differentiation，keratinization，),
TOP5的CC(cornified envelope，desmosome，cell-substrate junction，membrane microdomain，membrane raft),TOP5的MF(structural constituent of skin epidermis，serine-type peptidase activity，serine hydrolase activity，serine-type endopeptidase activity，protein kinase regulator activity)；
共显著富集到KEGG通路48条，其中TOP5的是（Hematopoietic cell lineage，Amoebiasis，Viral protein interaction with cytokine and cytokine receptor，Cytokine-cytokine receptor interaction，Staphylococcus aureus infection）（附件：KEGG-enrich.csv），并绘制气泡图（图2-2）
 

图 2-1： GO功能富集
注：图2-1为差异基因 的GO功能富集气泡图，结果以P.adjust排序选取 BP、MF、CC各展示前 5个通路，横坐 标为注释通路的基因的比例，纵坐标为通路，颜色是以 P.adjust的值来决定的，
大小是以注释通路的基因个数来决定的，详情大图见附件
 
图 2-2：KEGG功能富集
注：图2-2为差异基因的KEGG功能富集气泡图，横坐标为注释通路的基因的比例，纵坐标为通路，颜色是以P.adjust的值来决定的，大小是以注释通路的基因个数来决定的；

3.2.3、OVOL1的单基因GSEA 
Gene Set Enrichment Analysis（GSEA，基因集富集分析）用来评估一个预先定义的基因集的基因在与表型相关度排序的基因表中的分布趋势，从而判断其对表型的贡献。与GO 富集分析不同，
GSEA 分析不需要指定阈值（p 值或 FDR）筛选差异基因，可以将差异表达不显著却有着重要生物学意义的基因包含在内，避免遗漏关键信息。GSEA 功能富集其输入数据包含两部分，
一是已知功能的基因集（可以是GO 注释、MsigDB 的注释或其它符合格式的基因集定义），一是表达矩阵，GSEA 会对基因根据其与表型的关联度从大到小排序，然后判断。
基因集内每条注释下的基因是否富集于表型相关度排序后基因表的上部或下部，从而判断此基因集内基因的协同变化对表型变化的影响。 
根据方法的描述使用CARD16在TCGA-GBM数据集中的表达量中位数将数据集划分高低表达组，使用R包“DESeq2”对高低表达组进行差异分析，使用差异分析得到的logFC对数据集进行排序，
从MSigDB数据库（GSEA-msigdb.org/）获取KEGG、GO的背景基因集（c2.cp. kegg .v2022.1.Hs. symbols.gmt、c5.go.v2022.1.Hs.entrez.gmt）并使用R包“clusterProfiler”进行GSEA-GO、
以及GSEA-KEGG功能富集。 
3.2.3.1、GSEA-GO 
对数据集进行GO功能富集分析，设置pvalue<0.05，共显著富集到2173条通路（附件：GSEA_GO_result.csv），部分结果如表1、图2-3所示.
表1：数据集显著富集到的KEGG通路top10
					
GS	SIZE	ES	NES	NOM p-val
GOBP_T_CELL_SELECTION	44	0.66	1.6	0.035
GOBP_T_HELPER_17_CELL_DIFFERENTIATION	23	0.61	1.6	0.032
GOBP_TRACHEA_DEVELOPMENT	17	0.6	1.59	0.028
GOBP_EMBRYONIC_PLACENTA_MORPHOGENESIS	22	0.67	1.63	0.024
GOBP_MAINTENANCE_OF_CELL_POLARITY	17	0.53	1.63	0.023
GOBP_POSITIVE_T_CELL_SELECTION	32	0.67	1.62	0.022
GOBP_PHOSPHATIDIC_ACID_METABOLIC_PROCESS	28	0.51	1.62	0.021
GOBP_EPIDERMAL_CELL_DIFFERENTIATION	155	0.64	1.64	0.019
GOBP_HEMOGLOBIN_METABOLIC_PROCESS	16	0.69	1.69	0.016
GOBP_C21_STEROID_HORMONE_METABOLIC_PROCESS	34	0.59	1.6	0.012
GOBP_T_CELL_SELECTION	44	0.66	1.6	0.035
注：以Pvalue排序top10
 
图 2-3：GSEA-GO功能富集的结果top1(其他图见附件)
注：图2-3中为数据集的GSEA-GO 功能富集，图中展示以pvalue的大小排序，挑选p值最大的通路进行展示。每幅图都分成 3 个部分，最顶部的绿色折线为基因Enrichment Score的折线图。
纵轴为对应的Running ES, 在折线图中有个峰值，该峰值就是这个基因集的Enrichemnt score，峰值之前的基因就是该基因集下的核心基因。横轴代表此基因集下的每个基因，
对应第二部分类似条形码的竖线。 类似条形码的部分，为Hits，每条竖线对应该基因集下的一个基，为所有基因的rank值分布图，纵坐标为ranked list metric，即该基因排序量的值，
可理解为“公式化处理后的foldchange值”.

3.2.3.2、GSEA-KEGG 
对数据集进行GSEA-KEGG功能富集分析，设置pvalue<0.05，共显著富集到9条通路（附件：GSEA_KEGG_result.csv），部分结果如表2、图2-4所示，详情见附件。
表2：数据集显著富集到的KEGG通路top10
GS	SIZE	ES	NES	NOM P-VAL
KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION	240	0.49	1.45	0.046
KEGG_HEMATOPOIETIC_CELL_LINEAGE	82	0.56	1.48	0.044
KEGG_OLFACTORY_TRANSDUCTION	64	0.54	1.45	0.044
KEGG_GNRH_SIGNALING_PATHWAY	87	0.36	1.39	0.043
KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION	230	0.5	1.44	0.04
KEGG_CALCIUM_SIGNALING_PATHWAY	162	0.42	1.43	0.04
KEGG_LINOLEIC_ACID_METABOLISM	22	0.6	1.52	0.031
KEGG_RENIN_ANGIOTENSIN_SYSTEM	16	0.63	1.58	0.021
KEGG_ARACHIDONIC_ACID_METABOLISM	46	0.64	1.74	0
注：以Pvalue排序top9
 
图 2-4：GSEA-KEGG功能富集的TOP1结果
注：图2-4中为数据集的GSEA-KEGG 功能富集，图中展示以pvalue的大小排序，挑选p值最大的通路进行展示。每幅图都分成 3 个部分，最顶部的绿色折线为基因Enrichment Score的折线图。
纵轴为对应的Running ES, 在折线图中有个峰值，该峰值就是这个基因集的Enrichemnt score，峰值之前的基因就是该基因集下的核心基因。横轴代表此基因集下的每个基因，
对应第二部分类似条形码的竖线。 类似条形码的部分，为Hits，每条竖线对应该基因集下的一个基，为所有基因的rank值分布图，纵坐标为ranked list metric，即该基因排序量的值，
可理解为“公式化处理后的foldchange值”.

3.3.1、免疫浸润分析 
肿瘤微环境由肿瘤细胞、各种免疫浸润细胞、成纤维细胞、许多细胞因子和催化因子组成，是一个动态平衡的复杂系统。免疫应答在肿瘤生长、侵袭和转移中起重要作用，
因此肿瘤浸润免疫细胞（TIICs）可以成为化疗放射治疗的靶点。肿瘤免疫细胞浸润是指免疫细胞从血液中移向肿瘤组织，开始发挥它的作用，可以从肿瘤组织中分离出的浸润免疫细胞。
肿瘤中免疫细胞的浸润与临床结果密切相关，肿瘤中浸润的免疫细胞最有可能作为药物靶标来提高患者的生存率。单样本基因集富集分析（single sample gene set enrichment analysis, ssGSEA），
是针对单个样本无法做 GSEA 而设计的。最早是在 2009 年被提出（Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1）。
它可以由 GSVA 这个 R 包来实现。目前 ssGSEA 常被用于评估肿瘤免疫细胞浸润程度。 
为了研究胶质瘤数据集中OVOL1高低表达组（以OVOL1在数据集中的表达量中位数划分高低表达组）中免疫细胞浸润情况，次分析使用ssGSEA来进行免疫浸润分析，
从数据库获取28个免疫细胞的基因（http://cis.hku.hk/TISIDB/download.php），使用R包“GSVA”计算所有样本的28种免疫细胞含量（附件：ssGSEA_result.csv）。
经ssGSEA算法计算得到每种免疫细胞在各个样本中含量并使用R包“ggpubr”绘制上述高低表达组样本间的差异免疫细胞的箱线图，组间使用wilcox.test计算显著性结果如图3-1所示。
其中28种免疫细胞只有8种免疫细胞存在组间差异。
 
图3-1：组间差异的免疫细胞丰度
注：横坐标为免疫细胞，纵坐标为免疫细胞的丰度，米白色为Low样本的免疫细胞丰度，橙色为High样本的免疫细胞丰度，横坐标中标红的为组间存在显著差异的免疫细胞。图中*代表显著程度，
（“*” 表示 P < 0.05 ，“**”表示 P < 0.01，“***”表示 P < 0.001，“****”表示 P < 0.0001），红色标注的为显著的免疫细胞。

接着我们使用R（version 4.0.3）自带的“cor.test”计算，上述8个差异的免疫细胞("Central memory CD4 T cell","Activated dendritic cell","Myeloid derived suppressor cell"                ,"Type 2 T helper cell" ,"Memory B cell" ,"Effector memeory CD4 T cell","Eosinophil",   "Neutrophil")和目标基因OVOL1的相关性（图 2-12，附件：OVOL1_cor.csv）。共有6种免疫细胞与OOVOL1存在高度相关（cor>0.3&&pvalue<0.05）。
 
图2-12：相关性棒棒糖图
图中横坐标为pearson相关性，纵坐标为6个存在差异的免疫细胞基因,和2个没有差异的。右侧数值为p值，颜色也由p值决定.
3.3.2、免疫治疗分析 
TIDE免疫治疗分析
  近期，通过对TCGA（美国肿瘤基因组图谱）和PRECOG（基因组图谱预测临床结局）数据库的分析，研究人员揭示了不同免疫细胞类型的肿瘤浸润水平对患者总体生存期的影响。
  而预测肿瘤对免疫检查点抑制治疗的反应则需要首先了解肿瘤是如何逃避免疫系统的，因此，即使没有确切的免疫检查点抑制治疗数据，
  通过对公开的肿瘤遗传分子图谱的分析仍可能成为寻找免疫检查点抑制治疗预后评估指标的宝贵方法 近期的研究已经揭示了两种不同的肿瘤免疫逃避机制：
  在有些肿瘤中，虽然细胞毒性T细胞浸润程度高，但这些T细胞往往处于功能紊乱状态。而在另一些肿瘤中，免疫抑制因子可以祛除浸润在肿瘤组织中的T细胞。
  因此，Peng Jiang等人设计了一个全新的计算架构——肿瘤免疫功能障碍和排斥（TIDE）评分，用以综合这两种肿瘤免疫逃逸机制，此次研究涉及189个人类肿瘤研究，
  包含了33197份样本，其结果也被认为能够替代单一的生物标记物来有效预测免疫检查点抑制治疗的效果。
  首先将两个输入文件，处理成TIDE所需的格式，在TIDE中输入文件，然后检验免疫治疗响应性和亚型是否相关，p<0.05表示相关。为
  了更清楚的展示TIDE评分在高低表达组之间的差异，我们绘制了堆叠图.
