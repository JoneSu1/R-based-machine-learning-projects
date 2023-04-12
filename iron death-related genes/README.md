
The images presented in this project are not the real results, please do not use them directly. Please refer to the Code file if you need to reproduce it.


iron death-related genes
Background:

Coronary atherosclerotic heart disease (CHD), also known as coronary heart disease, is a group of diseases caused by atherosclerotic lesions in coronary vessels that narrow or block the lumen of the vessels, resulting in myocardial ischemia, hypoxia, or necrosis. CHD is one of the most lethal diseases worldwide, with China having the second highest number of CHD-related deaths globally [1]. The pathogenesis of coronary artery disease is complex, with no apparent symptoms in the early stage and only abnormal ST-T changes on exercise electrocardiogram (ECG) examination. Therefore, there is an urgent need for a biomarker detectable in peripheral blood to facilitate early detection of coronary artery disease.

Ferroptosis is an iron-dependent, novel form of programmed cell death, distinct from apoptosis, cell necrosis, and cell autophagy. The primary mechanism of ferroptosis involves catalyzing the high expression of unsaturated fatty acids on the cell membrane in the presence of divalent iron or lipoxygenase, inducing cell death through lipid peroxidation; additionally, it is characterized by a decrease in the expression of antioxidant systems (glutathione GSH and glutathione peroxidase 4-GPX4).

Professor Wang Fulong of Zhejiang University conducted a study using multiple cell death inhibitor treatments and cell death pathway-related knockout mice. They found that only Ferrostatin-1 (Fer-1), a specific inhibitor of ferroptosis, significantly reduced cardiotoxicity caused by the anticancer drug adriamycin (DOX) and effectively increased survival in mice, revealing that an essential mechanism of myocardial injury is ferroptosis. Ferroptosis has also been reported to be involved in cancer immunotherapy, with increased ferroptosis enabling immunotherapy to kill cancer cells more effectively [6]. Consequently, this experiment was conducted to analyze and validate the accuracy of ferroptosis-related genes as biomarkers for coronary heart disease through bioinformatics analysis.

Methods:

In this study, we first downloaded CHD RNA sequencing data and ferroptosis-related genes through the GEO database and FerrDb database. We obtained differentially expressed ferroptosis-related genes in CHD patients using the limma package for screening, constructed a PPI network of differential ferroptosis-related genes using the STRING database, and conducted GO and KEGG analysis of differentially expressed ferroptosis genes participating in biological processes. Subsequently, we further screened ferroptosis genes that could be used as biomarkers by LASSO and SVM-RFE, verified the expression of biomarkers in the GSE dataset and clinical samples, calculated the potential index of ferroptosis genes, determined the composition of 22 immune cells in CHD patients using CIBERSORT, and analyzed the correlation between biomarkers and immune cells. Finally, we performed biomarker analysis using GSEA and GSVA.

Conclusion:

Bioinformatics-based experiments were conducted to analyze and validate that ferroptosis-related genes can be employed as biomarkers for coronary heart disease.


Protocol design:

Flowchart
 
 
 ![workflow](https://user-images.githubusercontent.com/103999272/231342694-989b33dc-cda0-4561-be07-685865dadc01.png)

II, Operation steps (Note: pictures are for reference only)
(i) Data acquisition
1, GEO database to download coronary heart disease RNA sequencing data and clinical information (GEO database can refer to GSE12288, GSE59421);
2, FerrDb database to download 259 iron death-related genes, including 108 driver genes, 69 suppressor genes, and 111 gene markers.
(II) Screening of iron death differential genes
1, Using the "limma" package from R software to calculate the differentially expressed iron death-related genes in normal and coronary heart disease patients, the results were displayed in volcano and heat maps (screened according to the criteria of differential genes);
 ![1](https://user-images.githubusercontent.com/103999272/231343009-321fce99-1c54-462c-8245-b5ab95c6616a.png)

2, Correlation analysis between differentially expressed iron death-related genes in coronary artery disease;
 ![2](https://user-images.githubusercontent.com/103999272/231343014-2315f3bf-25f1-4491-9d4f-e3969d9bbd35.png)

3, Construction of PPI reciprocal network for differential iron death-related genes using STRING database (color differentiation of up and down calls of genes during composition).
 ![3](https://user-images.githubusercontent.com/103999272/231343026-db8d6307-478e-4577-813b-a70af9af03fc.png)

(iii) Functional enrichment
GO and KEGG analysis of differential iron death-associated genes (mainly to see if they are enriched to immune-related pathways), with the results displayed in bubble or bar graphs.
 ![4](https://user-images.githubusercontent.com/103999272/231343080-96d4e642-b730-4129-94b7-20137de8c099.png)

(iv) Screening of iron death biomarkers
The results of LASSO and SVM-RFE are intersected to obtain the key biomarkers (if there is no intersection, one of the methods is chosen, and there should be less than 5 markers).
  
  ![5](https://user-images.githubusercontent.com/103999272/231343128-24a8737d-5be6-43aa-ae81-c67802771ae6.png)

  
(v) ROC curve
ROC curves were constructed for coronary and normal samples in the GEO database and the area under the curve (AUC) was analyzed to assess the validity of iron death-related genes as biomarkers.
 
 ![6](https://user-images.githubusercontent.com/103999272/231343147-db5fba64-0332-42fa-88b8-8a1f876dc0d2.png)

(vi) Iron death potential index
The iron death potential index (FPI) was calculated based on the expression of key biomarkers in the GEO database using the R package "GSVA" to characterize the iron death activity in the samples (PMID: 32629423).
 
 ![7](https://user-images.githubusercontent.com/103999272/231343289-d0f47def-d060-44b3-8f34-1e58649e3deb.png)

(VII) Validation of iron death biomarkers
1, Expression of biomarkers in the GEO database in normal and coronary artery disease;
 ![8](https://user-images.githubusercontent.com/103999272/231343311-c72f70bc-f01c-4553-ab04-c198e9d0116f.png)

2, selection of other GSE datasets to validate the expression of biomarkers (refer to GSE20680);
 ![9](https://user-images.githubusercontent.com/103999272/231343317-08a96652-6f43-46ba-83d8-2436bc3d8b20.png)

3, analyze the expression of biomarkers in normal and coronary heart disease patients using qRT-PCR (Note: samples were sampled once for each patient, the number of samples per group was not less than 10 and repeated 3 times; samples were blood).
 
 
 ![10](https://user-images.githubusercontent.com/103999272/231343361-4b7b5d39-f470-49e3-91fa-590cef52b1e2.png)

(H) Biomarker enrichment analysis
1, GSEA analysis of biomarkers in the GEO database using the R package "clusterprofiler";
![WeChat Screenshot_20230412044342](https://user-images.githubusercontent.com/103999272/231343692-4ce33e05-e424-4157-a12c-d8154e443fb2.png)


2, GSVA analysis of biomarkers in the GEO database using the R package "GSVA".
 ![11](https://user-images.githubusercontent.com/103999272/231343382-8851016a-fcd8-4be5-855a-52b9845199ab.png)

(ix) Immunological infiltration
1, analysis of the proportion of 22 immune cells in the GEO database in patients with coronary artery disease by CIBERSORT algorithm;
  ![12](https://user-images.githubusercontent.com/103999272/231343392-b7af6eb4-2eab-47d3-98bf-eb749e50d3d7.png)


2, Exploration of the differences in immune cell infiltration seen in coronary artery disease and normal samples using PCA clustering plots;
  ![13](https://user-images.githubusercontent.com/103999272/231343416-fef09b24-7207-41bb-a926-0c77f27b1296.png)

3, drawing correlation heat maps and violin plots of 22 immune cell types to study the immune cell types (mainly looking at macrophages) that differ between normal and coronary heart disease patients;
 
![14](https://user-images.githubusercontent.com/103999272/231343425-3bd4b88d-4cf9-4b01-b2fd-ffb01b4b1c1e.png)
4, Correlation analysis of each key biomarker with immune cells.
  ![15](https://user-images.githubusercontent.com/103999272/231343445-f660bee4-281b-4c70-826a-d8e87e78c050.png)


