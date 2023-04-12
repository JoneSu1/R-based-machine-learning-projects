The analysis content includes:
1. Selection of relevant genes: Differential gene screening, disease trait-related gene module screening (WGCNA), and screening for endoplasmic reticulum stress-related
genes in disease traits and disease-related genes (obtaining a list of endoplasmic reticulum stress-related genes from public databases, and then intersecting with 
differential genes and disease-related module genes).
2. Validation of relevant genes: Disease-related endoplasmic reticulum stress genes are subjected to GO and KEGG analyses in R to examine whether the related BC and 
biological pathways meet the requirements. Due to the large number of obtained genes and weak BC correlation, we chose to use ROC diagnostic curves to screen for 
high-quality related genes (only retaining those with AUC values greater than 0.8), resulting in 17 genes remaining. We imported these 17 genes into the STRING database, 
constructed a PPI network, and then imported it into Cytoscape. Using the Maximal Clique Centrality (MCC) algorithm from an external plugin, we retained the 
top 5 core scoring genes as the most relevant endoplasmic reticulum stress genes for PD.
3. We also searched for external datasets to validate the correlation between these five genes and PD (by drawing single-gene expression box plots).

Dataset
Test dataset: Parkinson's BA9 brain tissue dataset GSE68719 downloaded from the GEO database (https://www.ncbi.nlm.nih.gov/) was used as the test set, with 44 normal samples and 29 disease samples.
Validation dataset: Parkinson's BA9 brain tissue dataset GSE20168 downloaded from the GEO database (https://www.ncbi.nlm.nih.gov/) was used as the validation set, with 15 normal samples and 14 disease samples.
Endoplasmic reticulum stress-related genes: A total of 868 genes with a relevance score greater than 7 were selected from the GeneCards database (https://www.genecards.org/) using the search term "endoplasmic reticulum stress-related genes".

Analysis Methods

1. The Parkinson's BA9 brain tissue dataset GSE68719 (RNA-seq data) downloaded from the GEO database (https://www.ncbi.nlm.nih.gov/) was used as the test set, with 44 normal samples and 29 disease samples. The DESeq2 package was employed for differential analysis between disease and control groups, and the ggplot2 package was utilized for result visualization.
2. The Parkinson's BA9 brain tissue dataset GSE20168 was downloaded from the GEO database (https://www.ncbi.nlm.nih.gov/) as the validation set, with 15 normal samples and 14 disease samples. Endoplasmic reticulum stress-related genes were searched for in the GeneCards database (https://www.genecards.org/) and genes with a relevance score greater than 7 were selected as endoplasmic reticulum stress-related genes (a total of 868 genes).
3. WGCNA analysis was conducted on the test set with PD as the trait, and the most PD-related module genes were obtained.
4. The online Venn tool was used to find the intersection of the test set's differentially expressed genes, WGCNA module genes, and downloaded endoplasmic reticulum stress-related genes, obtaining highly PD-related endoplasmic reticulum stress-related genes.
5. To explore the functions and related pathways of the core genes screened in step 4, the clusterprofile package in R was used for functional enrichment analysis (KEGG, GO), identifying potential biological functions and pathways related to the highly PD-related endoplasmic reticulum stress-related genes.
6. The roc() function in the pROC package was used to calculate the AUC values for the highly PD-related endoplasmic reticulum stress-related genes obtained in step 4 (separately for upregulated and downregulated genes). Seventeen genes with AUC values greater than 0.8 were selected as relevant genes.
7. The 17 relevant genes validated by ROC analysis were input into the STRING database for PPI network construction, and the top 5 genes ranked by the CytoHubba plugin in Cytoscape were selected as core genes.
8. To verify the effectiveness of the core genes related to endoplasmic reticulum stress in PD, the Parkinson's BA9 brain tissue dataset GSE20168 was used as the validation set, with 15 normal samples and 14 disease samples. Then, the expression boxplot of the key genes in normal and disease samples with drug intervention was drawn.
Analysis Results
3.1 Screening and functional prediction of highly relevant endoplasmic reticulum stress-related genes for PD
3.1.1 Screening of differentially expressed genes between Disease and Control
Differential expression analysis aims to identify genes with different expression levels between different sample groups and to further investigate differentially expressed genes. Firstly, we used the R language package "DESeq2" to perform differential analysis between disease samples and control samples in the dataset. The screening criteria for differentially expressed genes were P.value < 0.05 and |logFC| > 1. The results show that there are a total of 3,143 differentially expressed genes, including 1,829 upregulated genes and 1,314 downregulated genes. To make the results more intuitive, we plotted the following volcano plot (Figure 1) and heatmap (Figure 2) of differentially expressed genes.

![1](https://user-images.githubusercontent.com/103999272/231325020-5b601c8d-26e4-4d3d-8ede-a0d12dd9bfb0.png)

Figure 1. Disease vs Control differentially expressed gene volcano plot
Note: Each point in the graph represents a differentially expressed gene, with red solid dots representing upregulated genes, blue solid dots representing downregulated genes, and gray solid dots representing non-significant genes.

The heatmap shows that the expression of genes in Control samples and Disease samples is significantly different (Figure 2).


![2](https://user-images.githubusercontent.com/103999272/231325034-e54c037a-3584-434d-9902-be57251246cc.png)

Figure 2. Disease vs Control differentially expressed heatmap
Note: Each small square represents the expression of a gene, with the color representing the level of expression (red for high expression, blue for low expression). The squares at the top are blue for Control samples and red for Disease samples.

3.1.2 Screening of Parkinson's disease-related module genes (WGCNA)
First, a clustering method was used to screen for outlier samples, and no significantly deviant samples were found (Figure 3).

![3](https://user-images.githubusercontent.com/103999272/231325038-16322ed2-bf10-4f0b-9ef7-56387285c2a3.png)


Figure 3. Sample clustering for outlier detection
Note: The dendrogram shows the clustering analysis results of gene expression from different samples.

To construct a weighted co-expression network, we first need to determine a soft thresholding powers value to establish an adjacency matrix and make the gene distribution conform to the scale-free network according to the connectivity. The pickSoftThreshold() function in the WGCNA package can calculate the weighted network based on the given gene expression matrix according to the expression similarity and analyze the scale-free topology fitting index and connectivity information under the specified powers value (building a vector and passing it to powerVector) to help select a suitable powers value for network construction. Considering these indicators comprehensively, a suitable soft threshold powers value is selected, making R2 as large as possible while ensuring that connectivity is not too low. Generally, the power value is chosen when R2 is greater than 0.85, so in this example, we need to choose powers=16. The pickSoftThreshold() function itself also estimates an optimal reference value, which is also powers=16.


![4](https://user-images.githubusercontent.com/103999272/231325057-dc67bdb4-a06b-4a0c-a29c-3a85e2a93b40.png)

Figure 4. Fitting index vs power value scatter plot and average connectivity vs power value scatter plot;
Note: The vertical axis in the left graph shows the scale-free topology fitting index R2; the vertical axis in the right graph shows the average connectivity of the network.

Using the estimated most suitable power value, a scale-free network is constructed. First, the adjacency matrix is obtained based on the similarity of gene expression, and then the topological overlap matrix (TOM) is calculated. The co-expression similarity matrix between genes records the similarity of gene expression between them. Furthermore, based on this, the distance is calculated, and a clustering tree is plotted to see the overall distribution characteristics of similarity. We can also use the co-expression similarity pattern to divide into different clustering clusters. Genes within the same module have higher co-expression similarity, indicating that they have higher synergy, suggesting that they are involved in similar regulatory pathways or function in similar cellular regions. By using the dynamic tree cutting algorithm to segment modules, and setting the parameter minModuleSize to 100, we obtained 19 modules. Then, by setting the parameter MEDissThres to 0.30, we merged the modules and ultimately obtained 11 gene modules, as shown in Figure 5. Different color modules represent different modules, each of which is composed of genes with similar expression patterns. These modules are commonly referred to by "colors," such as "red module," "lightgreen module," etc., used as module names in the literature.


![5](https://user-images.githubusercontent.com/103999272/231325061-47802da2-2dbb-4a1e-877c-309b353ac6e2.png)

Figure 5. Gene expression clustering tree and co-expression topology heatmap
Note: From Figure 5, it can be clearly seen that the genes within each cluster have higher co-expression similarity, while the co-expression similarity between clusters is lower. Based on the co-expression patterns (cluster features), we divided them into 19 more distinct co-expression modules (Dynamic tree cut). For the divided co-expression modules, there are too many modules, and we used module eigengenes to further merge the modules, aggregating similar modules together, finally reducing the modules to 11 (Merged dynamic).

To determine the relationship between the disease and the modules, we performed a correlation analysis between the two, and the results are shown in Figure 6. The modules with a |correlation coefficient| > 0.5, pink and cyan, have a significantly strong correlation with the disease.



![6](https://user-images.githubusercontent.com/103999272/231325081-b3283124-1aa6-4c20-85e1-ab26f04d877f.png)

Figure 6. Association plot of co-expression modules and clinical phenotypes
Note: Each row represents a different gene module, and each column represents a different clinical phenotype (healthy and diseased). The first row of numbers in each cell represents the correlation coefficient, distinguished by red and green colors for positive and negative correlations, respectively, and the value in brackets in the second row is the significance p-value.

3.1.3 Screening of PD-related endoplasmic reticulum stress genes using online Venn tools
To screen for endoplasmic reticulum stress-related genes associated with Parkinson's disease, we took the intersection of the differentially expressed genes in PD, the highly correlated module genes, and the endoplasmic reticulum stress-related genes downloaded from the Genecards database using online Venn tools, and identified 68 common genes. We considered these common genes as PD-related endoplasmic reticulum stress genes.


![7](https://user-images.githubusercontent.com/103999272/231325101-8aa4cbf7-2ce1-43be-9c87-0191ad4322e3.png)


Figure 7. Venn diagram for the selection of PD-related endoplasmic reticulum stress genes
Note: The circles in the graph represent different gene lists, and the bar chart at the bottom indicates the number of genes in the different gene lists.



3.1.4 Exploration of PD-related endoplasmic reticulum stress gene-associated functions and pathways

To investigate the associated functions and pathways of the identified PD-related endoplasmic reticulum stress genes, we used the clusterProfiler package in R for functional enrichment analysis (KEGG, GO) to find potential biological functions and pathways involved in drug interventions. First, we mapped the PD-related endoplasmic reticulum stress genes to the GO database (http://www.geneontology.org/) terms and calculated the number of genes for each term to obtain a list of genes with specific GO functions and gene count statistics. Then, using a hypergeometric test, we identified GO entries significantly enriched in differentially expressed genes compared to the whole-genome background. After calculating p-values and adjusting them through FDR correction, we defined GO terms with a corrected p-value ≤ 0.05 as significantly enriched in the genes. We found 1332 enriched biological processes (BP), with the top 5 being: interleukin-6 production, regulation of interleukin-6 production, adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains, response to temperature stimulus (Figure 8). There were 77 cellular components (CC) where the target genes mainly function, with the top 5 being: endocytic vesicle, endocytic vesicle membrane, blood microparticle, endoplasmic reticulum lumen, membrane raft (Figure 8). The target genes had 88 molecular functions (MF), with the top 5 being: heat shock protein binding, amide binding, peptide binding, chaperone binding, misfolded protein binding (Figure 8).

![8](https://user-images.githubusercontent.com/103999272/231325799-60f28981-c917-4be0-b8f5-32b0f127c335.png)


Figure 8. Bubble plot of GO enrichment analysis of PD-related endoplasmic reticulum stress genes
Note: Figure 8 shows the top 5 functions in the GO enrichment analysis with the smallest P-values. The bubble plot displays various key statistics of the enrichment analysis results. The Y-axis corresponds to the function, the X-axis represents the ratio of differentially expressed genes to all genes in a specific function, the bubble size indicates the number of differentially expressed genes in the function, and the bubble color changes from blue to purple to red, with smaller p-values indicating a higher degree of significance.

We then used the major public database for Pathways, KEGG. We conducted a pathway significance enrichment analysis using KEGG Pathway as the unit, applying a hypergeometric test to identify pathways significantly enriched in differentially expressed genes compared to the whole-genome background. After calculating p-values and adjusting them through FDR correction, we found 92 pathways significantly enriched in the drug intervention genes with a corrected p-value ≤ 0.05, with the top 5 being: Lipid and atherosclerosis, TNF signaling pathway, Fluid shear stress and atherosclerosis, Kaposi sarcoma-associated herpesvirus infection, Cellular senescence.


![9](https://user-images.githubusercontent.com/103999272/231325808-2879219b-8e5b-4654-a6aa-fff83bfa48eb.png)


Figure 9. Bubble plot of KEGG enrichment analysis for downregulated genes
Note: Figure 9 shows the top 10 pathways in the KEGG enrichment analysis with the smallest P-values. The bubble plot displays various key statistics of the enrichment analysis results. The Y-axis corresponds to the pathway, the X-axis represents the ratio of differentially expressed genes to all genes in a specific pathway, the bubble size indicates the number of differentially expressed genes in the pathway, and the bubble color changes from blue to purple to red, with smaller p-values indicating a higher degree of significance.

3.2 Selection of PD-related endoplasmic reticulum stress core genes
3.2.1 ROC diagnostic screening of high-quality genes
When choosing diagnostictests, many researchers face a difficult trade-off between sensitivity and specificity. So, is it possible to consider both sensitivity and specificity characteristics and evaluate the accuracy of diagnostic tests based on a single indicator? In 1971, Lusted proposed the receiver operating characteristic curve (ROC curve) to describe the intrinsic accuracy of diagnostic tests, which has been widely applied. The ROC curve is plotted with the true positive rate (sensitivity) as the vertical axis and the false positive rate (1-specificity) as the horizontal axis. Each point corresponds to a cut-off point in the diagnostic test, and connecting these possible points creates the empirical ROC curve. This method helps researchers analyze the clinical accuracy of diagnostic tests simply and intuitively. The most important aspect of the ROC curve is the area under the curve (AUC). An AUC of 0.5 indicates random classification with zero discriminative ability, while an area closer to 1 indicates stronger discriminative ability, and an area equal to 1 represents perfect discrimination.

First, we annotated the 68 relevant genes identified in the previous step for upregulation and downregulation in the differentially expressed file. Of these, 57 were upregulated (15 of which had AUC > 0.8), and 11 were downregulated (2 of which had AUC > 0.8). We then used the roc() function in the pROC package to calculate the AUC values for each gene in the highly related PD endoplasmic reticulum stress genes, selecting 17 genes with AUC values greater than 0.8 as highly related genes. Table 1 shows the AUC values for the 17 related genes, and Figure 10 displays the ROC curve for the top AUC gene, AKAP6. (See attachments for others.)

Table 1: AUC values for the 17 related genes
Gene	AUC_value（great than 0.8）
AKAP6	0.857367
MBTPS2	0.844828
CLU	0.811912
HSPA1A	0.833072
HSPB1	0.860502
CAPN2	0.808777
BAG3	0.802508
MAPKAPK2	0.819749
SERPINH1	0.82837
TNFRSF1A	0.80094
MAP2K3	0.826803
CD4	0.826019
TNFRSF10B	0.804075
TGFB1	0.823668
CASP7	0.80721
ABCD1	0.80094
AGER	0.818966
 
 
 ![10](https://user-images.githubusercontent.com/103999272/231327086-77d9def4-e318-49b4-a623-a97216615641.png)

Figure 10. ROC curve of AKAP6
Note: The horizontal axis (X-axis) represents 1 – specificity, also known as the false positive rate (error rate), and the closer the X-axis is to zero, the higher the accuracy; the vertical axis (Y-axis) is called sensitivity, also known as the true positive rate (sensitivity), and the larger the Y-axis, the better the accuracy. The curve divides the entire graph into two parts, with the area below the curve referred to as the AUC (Area Under Curve), which represents the predictive accuracy. The higher the AUC value, the larger the area under the curve, indicating a higher prediction accuracy. The curve is more accurate the closer it is to the upper left corner (smaller X, larger Y). The black point on the curve represents the most comprehensive data point for disease impact when the AUC value is greater than 0.5, with the best threshold of 0.614 and a specificity of 0.966 for PD traits, providing the overall optimal performance.

3.2.2 Screening of PD-related endoplasmic reticulum stress core genes
STRING is a database of known and predicted protein-protein interactions. Interactions include direct (physical) and indirect (functional) associations; they are derived from computational predictions, knowledge transfers between organisms, and interactions aggregated from other (primary) databases. The interactions in STRING come from five major sources: genomic predictions, high-throughput experiments, (conserved) co-expression experiments, automated text mining, and database-related knowledge. The STRING database currently covers 24,584,628 proteins from 5,090 species. Construction of PPI protein interaction networks (PPI, Protein-Protein Interaction Networks) allows us to examine the relationships between these genes and identify key genes.

To screen for PD-related endoplasmic reticulum stress core genes, we first imported the 17 related genes selected in the previous ROC analysis into the STRING database (https://string-db.org/) to construct a PPI network (retaining only connected nodes, with 14 drug intervention genes connected) (Figure 9).


![11](https://user-images.githubusercontent.com/103999272/231327100-0da52c25-3738-47aa-91fa-f344974589cc.png)

Figure 11. Protein interaction network of PD-related endoplasmic reticulum stress genes
Note: In the protein interaction network shown in Figure 11, each circle represents a gene (node), each line represents an interaction between a pair of proteins, and the thickness of the line is based on the confidence level of the interaction in the STRING database (only interactions with a confidence level >0.4 are displayed).

We then imported the PPI network interaction file into Cytoscape and used the CytoHubba plugin in Cytoscape to screen for core genes. The CytoHubba plugin assigns a value to each gene (node) and ranks them by hub-node based on different calculation algorithms. We used the Maximal Clique Centrality (MCC) method and retained the top 5 scored genes as core genes (Figure 12).


![12](https://user-images.githubusercontent.com/103999272/231327117-f74eee82-29a6-486f-9cd2-2fc9efc24e68.png)


Figure 12. Core genes predicted by the CytoHubba plugin
Note: In the core gene interaction network shown in Figure 12, each circle represents a gene (node), each line represents an interaction between a pair of proteins, and the darker the color of the node, the higher the MCC score.

3.3 Validation of PD-related endoplasmic reticulum stress core genes in an external dataset
Using the Parkinson's BA9 brain tissue dataset GSE20168 as a validation set, with 15 normal samples and 14 disease samples. We then plotted the expression box plots of key drug intervention genes in normal and disease samples. As shown in Figure 13, the expression of 5 PD-related endoplasmic reticulum stress core genes in dataset GSE20168 is consistent with the expression trends in the test set. These 9 genes are: HSPBAP1, TNFRSF1A, CD4, HSPA1A, and TNFRSF10, which validates the significance of our selected core genes.


![13](https://user-images.githubusercontent.com/103999272/231327127-08b9d18f-1a32-4f74-abcf-1759000d22f7.png)


Figure 13. Expression box plot of PD-related endoplasmic reticulum stress core genes
Note: The horizontal axis represents the tissue classification of gene expression (Disease (orange) and Control (blue)), the vertical axis represents the expression values of the hub-genes, and the numbers at the top represent the p-values of the differences between the gene expression in the tumor and control groups calculated using the Wilcoxon test algorithm.
