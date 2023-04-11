# R-based-machine-learning-projects
Included analyses: Analysis of variance, WGCNA analysis, GO/KEGG enrichment analysis, GSVA analysis, correlation expression analysis, ROC diagnostic curve, Lasso regression analysis, XGbost analysis, random forest, single-factor COX, multi-factor COX, random vector machine SVM analysis, K-M survival analysis,

In this Code file, I have provided sufficient annotations for code usability. In this Read file, I will focus on introducing the project's rationale and results.

The objective of this project is to explore genes associated with the OVOL1 gene and to investigate the biological functions related to these genes, as well as potential immunotherapeutic targets.
![workflow](https://user-images.githubusercontent.com/103999272/231105959-f3f26d9d-b41d-4483-8093-eec3c2694e8a.png)



1 Dataset
Download the GSE46517 dataset of malignant melanoma from the GEO database (https://www.ncbi.nlm.nih.gov/geo/), with a total of 195 samples. There are 73 metastatic melanoma samples, 31 primary melanoma samples, and 15 normal samples.
Target gene: OVOL1.

2 Analysis methods

Download the GSE46517 dataset of malignant melanoma from the GEO database (https://www.ncbi.nlm.nih.gov/geo/), with a total of 195 samples. There are 73 metastatic melanoma samples, 31 primary melanoma samples, and 15 normal samples.

Differential analysis: First, group samples based on the median expression level of the target gene OVOL1 (high expression, low expression). Then, use the R package "limma" to perform differential analysis on the dataset, with a differential gene selection threshold of |logFC|>0.5 & P.adjust <0.05.

Functional enrichment: Use the R package "clusterProfiler" to perform GO and KEGG functional enrichment on the obtained differentially expressed genes.

Single-gene GSEA functional enrichment: Use the R package "clusterProfiler" to perform single-gene GSEA analysis on the OVOL1 dataset.

Immune infiltration analysis: Calculate the immune cell content in each TCGA-GBM sample using ssGSEA and analyze the differences between the high and low expression groups of the target gene. For the differential immune cells, further calculate the correlation with CARD16.

Immune checkpoint correlation analysis: Calculate the correlation between CARD16 and immune checkpoints. Immunotherapy analysis: Use the online database TIDE (https://cloud.genepattern.org/gp/pages/login.jsf) to test the correlation between immunotherapy responsiveness and subtypes, and visualize the results in R.

Immunotherapy analysis: Use the online database GenePattern (https://cloud.genepattern.org/gp) to predict the immunotherapy responsiveness of subtypes.
3 Analysis results
3.1 Differential analysis
First, download the GSE46517 dataset of malignant melanoma from the GEO database (https://www.ncbi.nlm.nih.gov/geo/), with a total of 195 samples. There are 73 metastatic melanoma samples, 31 primary melanoma samples, and 15 normal samples.
Then, group samples based on the median expression level of the target gene OVOL1 (5.974117) (high expression, low expression), and use the R package "limma" to perform differential analysis on the dataset, with a differential gene selection threshold of |logFC|>0.5 & P <0.05.
A total of 1082 differentially expressed genes were screened, including 520 downregulated genes and 562 upregulated genes. To gain an overall understanding of the distribution of differentially expressed genes, a volcano plot is used to display their distribution as shown in Figure 1-1.

Figure 1. Differential volcano plot
Note: (Horizontal axis represents the difference multiple (High/low, taking the logarithm base 2), vertical axis represents -log10(p). Each point in the figure represents a gene, and blue and red points represent significantly differentially expressed genes.
Red points indicate upregulated gene expression (High-expression samples compared to Low-expression samples), blue points indicate downregulated gene expression (High-expression samples compared to Low-expression samples), and gray points indicate no significant differences between these genes.)

A heatmap is an intuitive visualization method for analyzing the distribution of experimental data. It can be used for quality control of experimental data and to intuitively display the expression level differences of key research objects. It can also cluster data and samples to observe sample quality. There are various forms of heatmaps, but the basic elements are universal, as shown in Figure 1-2.

Figure 1-2: Heatmap of differentially expressed genes between sample groups
(Each small square represents a gene, and its color represents the gene expression level; the deeper the color, the higher the expression (red for high expression, blue for low expression). The first row represents the sample grouping, with blue representing Low-expression samples and red representing High-expression samples. Each row represents the expression level of each gene in different samples, and each column represents the expression level of all differentially expressed genes in each sample. The dendrogram on the left represents the clustering analysis results of different genes from different samples.)

3.2.2 Functional enrichment
Use the obtained 1056 differentially expressed genes to perform GO and KEGG pathway enrichment analysis, searching for common functions and related pathways among a large number of genes within the gene set. Use the cumulative hypergeometric distribution to analyze whether a group of genes is overrepresented in a particular functional node, with the following formula:

Here, N is the total number of genes in the annotation system, n is the number of annotated genes in the node or pathway being examined, M is the size of the differentially expressed gene set, and x is the number of genes in the intersection between the gene set and the node or pathway.

The Kyoto Encyclopedia of Genes and Genomes (KEGG) is a useful database resource for understanding high-level functions and biological systems (such as cells, organisms, and ecosystems) from molecular-level information, especially large molecular datasets generated by genome sequencing and other high-throughput experimental techniques. The Gene Ontology (GO) system consists of three parts: biological processes (BP), molecular functions (MF), and cellular components (CC).

We used the R package "clusterProfiler" to perform GO and KEGG functional enrichment analysis, searching for common functions and related pathways among differentially expressed genes. During the enrichment analysis process, we set the threshold for GO and KEGG pathway enrichment analysis as P.adjust <0.05. A total of 1263 significant GO function entries were enriched, including 1042 BPs, 101 CCs, and 120 MFs (Attachment: GO-enrich.csv), and bubble plots were drawn (Figure 2-1). The top 5 BPs were epidermis development, skin development, keratinocyte differentiation, epidermal cell differentiation, and keratinization. The top 5 CCs were cornified envelope, desmosome, cell-substrate junction, membrane microdomain, and membrane raft. The top 5 MFs were structural constituent of skin epidermis, serine-type peptidase activity, serine hydrolase activity, serine-type endopeptidase activity, and protein kinase regulator activity.

A total of 48 significant KEGG pathways were enriched, and the top 5 were Hematopoietic cell lineage, Amoebiasis, Viral protein interaction with cytokine and cytokine receptor, Cytokine-cytokine receptor interaction, and Staphylococcus aureus infection (Attachment: KEGG-enrich.csv). Bubble plots were drawn (Figure 2-2).




Figure 2-1: GO Function Enrichment
Note: Figure 2-1 is a bubble chart of the GO function enrichment of differentially expressed genes, with the results sorted by P.adjust and the top 5 pathways for BP, MF, and CC displayed. The horizontal axis represents the proportion of annotated genes in the pathway, the vertical axis represents the pathway, the color is determined by the P.adjust value, and the size is determined by the number of annotated genes in the pathway. See the attachment for a larger image.

Figure 2-2: KEGG Function Enrichment
Note: Figure 2-2 is a bubble chart of the KEGG function enrichment of differentially expressed genes, with the horizontal axis representing the proportion of annotated genes in the pathway, the vertical axis representing the pathway, the color determined by the P.adjust value, and the size determined by the number of annotated genes in the pathway.

3.2.3 OVOL1 Single-Gene GSEA
Gene Set Enrichment Analysis (GSEA) is used to evaluate the distribution trend of a pre-defined gene set in a gene list sorted by phenotype-related relevance, thereby determining its contribution to the phenotype. Unlike GO enrichment analysis, GSEA does not require a specified threshold (p-value or FDR) to filter differentially expressed genes, thus including genes with insignificant differential expression but important biological significance and avoiding the omission of key information. GSEA function enrichment input data consists of two parts: a known function gene set (such as GO annotations, MsigDB annotations, or other gene set definitions that meet the format requirements) and an expression matrix. GSEA sorts genes by their association with the phenotype and determines whether genes under each annotation in the gene set are enriched in the upper or lower part of the phenotype-related sorted gene list, thereby determining the impact of the coordinated changes of genes within the gene set on the phenotypic changes.

According to the method description, CARD16 was used to divide the TCGA-GBM dataset into high and low expression groups based on the median expression level. The R package "DESeq2" was used for differential analysis of high and low expression groups, and the logFC obtained from the differential analysis was used to sort the dataset. The KEGG and GO background gene sets (c2.cp.kegg.v2022.1.Hs.symbols.gmt, c5.go.v2022.1.Hs.entrez.gmt) were obtained from the MSigDB database (GSEA-msigdb.org/) and the R package "clusterProfiler" was used for GSEA-GO and GSEA-KEGG function enrichment.

3.2.3.1 GSEA-GO
GO function enrichment analysis was performed on the dataset, setting p-value < 0.05, and a total of 2173 pathways were significantly enriched (attachment: GSEA_GO_result.csv). Partial results are shown in Table 1 and Figure 2-3.

Table 1: Top 10 KEGG pathways significantly enriched in the dataset
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
(Note: Top 9 sorted by P-value)

Figure 2-4: TOP1 results of GSEA-KEGG functional enrichment
Note: Figure 2-4 shows the GSEA-KEGG functional enrichment of the dataset, displayed in order of the P-value size, selecting the pathway with the highest P-value. Each graph is divided into 3 parts, with the green line at the top representing the gene Enrichment Score curve. The vertical axis represents the Running ES, and there is a peak in the curve, which is the Enrichment Score of this gene set. The genes before the peak are the core genes in this gene set. The horizontal axis represents each gene in this gene set, corresponding to the second part of the graph, which is similar to a barcode. The barcode-like part, called Hits, consists of vertical lines representing one gene in the gene set, and the vertical coordinate is the ranked list metric, i.e., the value of the gene ranking, which can be understood as the "formula-processed fold change value."

3.3.1 Immune infiltration analysis
The tumor microenvironment is a dynamic and complex system composed of tumor cells, various immune infiltrating cells, fibroblasts, and numerous cytokines and catalytic factors. Immune responses play a crucial role in tumor growth, invasion, and metastasis, making tumor-infiltrating immune cells (TIICs) potential targets for chemotherapy and radiotherapy. Tumor immune cell infiltration refers to immune cells migrating from the blood to the tumor tissue and exerting their functions, and infiltrating immune cells can be isolated from tumor tissues. The infiltration of immune cells in tumors is closely related to clinical outcomes, and immune cells infiltrating tumors are most likely to serve as drug targets to improve patient survival rates. Single sample gene set enrichment analysis (ssGSEA) is designed for single samples that cannot perform GSEA. It was first proposed in 2009 (Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1) and can be implemented using the GSVA R package. Currently, ssGSEA is commonly used to assess the degree of tumor immune cell infiltration.

To investigate the immune cell infiltration in glioma datasets with high and low OVOL1 expression (divided into high and low expression groups based on the median expression of OVOL1 in the dataset), this sub-analysis used ssGSEA for immune infiltration analysis, obtained 28 immune cell genes from the database (http://cis.hku.hk/TISIDB/download.php), and calculated the abundance of 28 immune cells in all samples using the R package "GSVA" (Appendix: ssGSEA_result.csv). The ssGSEA algorithm calculated the abundance of each immune cell in each sample and used the R package "ggpubr" to draw box plots of the differential immune cells between the high and low expression group samples, with group differences calculated using the Wilcox.test, as shown in Figure 3-1. Among the 28 immune cells, only 8 immune cells showed group differences.

Figure 3-1: Differential immune cell abundance between groups
Note: The horizontal axis represents immune cells, the vertical axis represents the abundance of immune cells, beige represents the immune cell abundance of Low samples, and orange represents the immune cell abundance of High samples. The red marks on the horizontal axis indicate immune cells with significant group differences. The asterisk in the figure indicates the significance level, with "" indicating P < 0.05, "" indicating P < 0.01, "" indicating P < 0.001, and "****" indicating P < 0.0001. Redannotations indicate significant immune cells.

Next, we used R's built-in "cor.test" (version 4.0.3) to calculate the correlation between the 8 differential immune cells ("Central memory CD4 T cell", "Activated dendritic cell", "Myeloid derived suppressor cell", "Type 2 T helper cell", "Memory B cell", "Effector memory CD4 T cell", "Eosinophil", and "Neutrophil") and the target gene OVOL1 (Figure 2-12, Appendix: OVOL1_cor.csv). A total of 6 immune cells showed a high correlation with OVOL1 (cor > 0.3 && P-value < 0.05).

Figure 2-12: Correlation lollipop plot
In the figure, the horizontal axis represents Pearson correlation, the vertical axis represents the 6 differential immune cell genes and 2 without differences. The values on the right side represent P-values, and the color is determined by the P-value.
3.3.2 Immune therapy analysis
TIDE immune therapy analysis
Recently, by analyzing TCGA (The Cancer Genome Atlas) and PRECOG (Genomic Landscape Predictive of Clinical Outcomes) databases, researchers have revealed the impact of different immune cell types infiltrating tumor levels on patients' overall survival. To predict the tumor's response to immune checkpoint inhibitor therapy, it is necessary first to understand how the tumor evades the immune system. Therefore, even without precise immune checkpoint inhibitor therapy data, analyzing publicly available tumor genetic molecular landscape data may still be a valuable method for finding prognostic indicators for immune checkpoint inhibitor therapy. Recent studies have revealed two different tumor immune evasion mechanisms: in some tumors, although cytotoxic T cell infiltration is high, these T cells are often functionally disordered. In other tumors, immune inhibitory factors can eliminate T cells infiltrating the tumor tissue. Thus, Peng Jiang and colleagues designed a new computational framework called Tumor Immune Dysfunction and Exclusion (TIDE) score to integrate these two tumor immune escape mechanisms. This study involved 189 human tumor studies, covering 33,197 samples, and its results are considered to replace single biomarkers effectively in predicting the effects of immune checkpoint inhibitor therapy.

First, process the two input files into the format required by TIDE, input the files into TIDE, and then test whether immune therapy responsiveness and subtypes are related, with P < 0.05 indicating relevance. To better illustrate the difference in TIDE scores between high and low expression groups, we created a stacked plot.


  
Note: The gene expression value should be normalized toward a control sample which could be either normal tissues related with a cancer type or mixture sample from diverse tumor samples.

The log2(RPKM+1) values from a RNA-seq experiment may not be meaningful unless a good reference control is available to adjust the batch effect and cancer type difference.

In our study, we used the all sample average in each study as the normalization control.
