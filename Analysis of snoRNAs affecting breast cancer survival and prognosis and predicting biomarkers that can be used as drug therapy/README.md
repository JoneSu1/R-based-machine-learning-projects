To prevent misuse of the project results, I will present the project with an example diagram. If you need to reproduce it, see the code for details.

Breast cancer (BRCA) is one of the most common cancers in women, and its incidence has been steadily increasing worldwide in recent years, 
ranking second in female malignant tumors in terms of mortality. BRCA is a complex disease involving a series of epigenetic changes. 
Despite the continuous improvement of diagnostic and therapeutic measures for breast cancer over the past decade, the 5-year survival
rate of breast cancer patients remains low. Therefore, exploring the pathogenesis and progression mechanism of breast cancer and identifying 
new biomarkers for breast cancer treatment is of great significance.

Small nucleolar RNAs (snoRNAs) are small non-coding RNAs widely present in the nucleoli of eukaryotic cells. They are 0-300 nt in length 
and have conserved structural elements. SnoRNAs mainly consist of two families: C/D box snoRNA and H/ACA box snoRNA. The main function of 
snoRNA is to guide ribosomal RNA (rRNA) methylation and pseudouridylation through base pairing. There is also evidence indicating that snoRNA 
participates in the post-transcriptional modification of snRNA, tRNA, and mRNA. The complex formed between snoRNAs and ribonucleoproteins (RNPs) 
is called snoRNPs. Studies have shown that the synergistic effect of snoRNPs and snoRNA can participate in the occurrence and development of cancer
by affecting ribosomes and protein translation. With the continuous development of next-generation sequencing technology, more and more studies have 
shown that snoRNAs are abnormally regulated and play a role in cancer. Currently, the mechanism and clinical value of snoRNA in breast cancer have not
been fully elucidated. Therefore, this study analyzed the biological role of snoRNA in breast cancer by analyzing the RNA-seq sequencing data of breast
cancer from the TCGA database.

In this study, RNA-seq sequencing data and clinical information of BRCA patients were downloaded from the TCGA database. Differential expression 
of snoRNA between the normal control group and the breast cancer patient group was obtained by analyzing the snoRNA extracted from the RNA-seq sequencing data. 
Survival analysis was performed using univariate and multivariate Cox regression analysis, and a survival risk scoring model was constructed. The risk
score of each BRCA patient was calculated, and the patients were divided into high-risk and low-risk groups based on the median value. Kaplan-Meier (KM) 
survival curve analysis was used to analyze the survival difference between the high-risk and low-risk groups, and ROC curve was used to further evaluate 
the accuracy of the prognostic model. A bar line chart was used to predict the survival rate of BRCA patients at 1, 3, and 5 years, and the relationship between
the prognostic model and clinical characteristics was analyzed. SNORi database was used to analyze the correlation between key snoRNAs and copy number 
variations and methylation sites, as well as the functional enrichment analysis of snoRNPs. Finally, CMap was used to analyze potential therapeutic drugs 
for breast cancer.

In conclusion, this study analyzed and validated the role of snoRNA in the survival and prognosis of breast cancer.

Program Design
I. Flowchart
![workflow](https://user-images.githubusercontent.com/103999272/231333135-b39fa6b1-eb80-4e32-ba23-13efa1a3928b.png)

Second, the operation steps (Note: the picture is for reference only)

 Data Download

Download RNA-seq data and clinical information of breast cancer patients from TCGA database.
Divide the breast cancer patient data in the TCGA database into a training set and a validation set in a 7:3 ratio.
 Differential Expression Screening of snoRNA
Figure 1A: Extract snoRNA data information from the RNA-seq data of the breast cancer training set, use the R package "limma" to analyze 
the differential expression of snoRNA between normal controls and breast cancer patients, and display the results in a volcano plot.

Figure 1B: Perform univariate Cox regression analysis on the differential expression of snoRNA in the training set, and display the results
with a forest plot for P<0.05.

Figure 1C: Perform multivariate Cox regression analysis on the differential expression of snoRNA with P<0.05 in the univariate analysis in the training set, 
display the results with a forest plot, and construct a survival risk scoring model (Note: snoRNA used in constructing the model as a biomarker).
![1](https://user-images.githubusercontent.com/103999272/231333523-a179bb32-7622-497d-b316-5300fbf55a19.png)

Figure 1D: Determine the cut-off point for risk scoring by the median, and divide breast cancer patients into high-risk and low-risk groups
(Note: this is done in both the training and validation sets simultaneously).

![2](https://user-images.githubusercontent.com/103999272/231334488-12d7c635-e85f-460f-aef8-2ce1ee4b1223.png)

 Figure 1E: Comparison of survival differences between high and low risk groups by log-rank test (Note: performed in both training and validation sets)
 
 ![3](https://user-images.githubusercontent.com/103999272/231334493-6911bc13-cbcb-4e94-8051-b62723d23933.png)

Figure 1F: Use ROC curve to predict the survival risk scoring model for 1, 3, and 5 years (Note: this is done in both the training and validation sets simultaneously).

 Correlation Between Survival Risk Scoring Model and Clinical Characteristics (Training Set)

![4](https://user-images.githubusercontent.com/103999272/231334567-34ab6e69-254f-428b-b2be-307514b8e793.png)

Figure 2A: Heatmap showing the correlation between survival risk scoring and breast cancer clinicopathological features (selected according to the clinical information in the TCGA database);

![5](https://user-images.githubusercontent.com/103999272/231334593-9684b5a6-babc-4cb1-b127-5250dd9cdd8d.png)

Figure 2B: Box plot showing the correlation between biomarkers and breast cancer clinicopathological features (selected according to the clinical information in the TCGA database).

 Independent Prognostic Analysis of Risk Scoring Model (Training Set)
 
 ![6](https://user-images.githubusercontent.com/103999272/231334621-e6f98e7d-7b32-407d-afbf-e5a06134d8c7.png)

Figure 3A: Univariate and multivariate Cox regression analysis to determine whether the risk score is an independent prognostic factor;

![7](https://user-images.githubusercontent.com/103999272/231334639-640e3999-292b-4d3f-88d8-91ec572e0be9.png)

Figure 3B: Construct a bar line chart based on the risk score and clinical pathological parameters in the TCGA database to analyze the possible 
1, 3, and 5-year survival rates of breast cancer patients;

![8](https://user-images.githubusercontent.com/103999272/231334670-745df3e7-3fc8-45d6-9ebf-b1b8b9881f1d.png)

Figure 3C: Construct a calibration curve for the bar line chart to test the consistency between the predicted survival after bias correction and the
observed survival.

 Analysis of the Correlation between snoRNA and CNV and Methylation
Purpose: CNV and methylation have been reported to be involved in the expression of snoRNAs in multiple cancers. Therefore, analyze the correlation 
between snoRNA and CNV and DNA methylation.

![9](https://user-images.githubusercontent.com/103999272/231334717-5be40452-eef8-4b0c-bd95-ef10e029702e.png)

Figure 4A: Divide snoRNA into H/ACA box and C/D box based on their structure and main functions, and display the distribution of biomarker snoRNAs 
in different types with a graph;

![10](https://user-images.githubusercontent.com/103999272/231334750-e82609d9-5f0f-4eeb-a748-a126dd399a91.png)

Figure 4B: Analyze the correlation between biomarker snoRNA and copy number variation using the SNORi database;

![11](https://user-images.githubusercontent.com/103999272/231334774-72918f5a-b221-43e8-83d4-e5471220b73c.png)

Figure 4C: Display the interaction between biomarker snoRNA and methylation sites with a Sankey diagram;
![12](https://user-images.githubusercontent.com/103999272/231334817-c2b0e212-737e-4549-9b35-2df2dea8d76f.png)

snoRNPs Functional Analysis
![5A](https://user-images.githubusercontent.com/103999272/231335491-3ddb9622-7437-4b1f-a41d-bce08a713f0b.png)

Figure 5A: Determine the RNPs associated with each biomarker snoRNA using Pearson correlation analysis, and display the results as shown in 
the legend (obtained from the SNORic database);
![5B](https://user-images.githubusercontent.com/103999272/231335504-59214930-b2ee-4d92-b57a-c8c913a222dd.png)


Figure 5B: Construct a PPI interaction network of biomarker snoRNAs and their corresponding snpRNP genes using the String software;

![5C](https://user-images.githubusercontent.com/103999272/231335514-e4e637ad-0462-43a9-b4bf-3d4eed4abaa7.png)


Figure 5C: Perform functional enrichment analysis of all biomarker-associated snoRNPs using GO and KEGG;
Identification of Potential Therapeutic Drugs based on CMap Analysis

![6-1](https://user-images.githubusercontent.com/103999272/231335534-6f7947cf-b4c3-4ac2-be4e-ad8386cd8e4b.png)

Table 1: Screen potential compounds targeting breast cancer using CMap. Upload all biomarker-associated snoRNPs to the CMap database, and select
compounds with a certain threshold (reference values: enrichment score absolute value ≥ 0.5, P < 0.05) as potential therapeutic drugs for breast cancer. 
Display the detailed information of potential compounds in a list;


![6-2](https://user-images.githubusercontent.com/103999272/231335550-711726a7-8aff-483c-b977-832d4f85103f.png)

Figure 6: Search for the structure of potential therapeutic drugs using NCBI's PubChem.
Experimental Validation


![7A](https://user-images.githubusercontent.com/103999272/231335578-a1385b5c-f25a-4e8b-aa7b-0176b1043abe.png)

Figure 7A: Analyze the expression levels of biomarker snoRNAs in the serum of normal controls and breast cancer patients using qRT-PCR 
(if the sample size permits, this step can be expanded to include different subtypes of normal and breast cancer samples. Note: samples 
are collected from each patient once, with at least 10 samples per group and repeated 3 times);

![7B](https://user-images.githubusercontent.com/103999272/231335604-d81efc01-3810-4287-ae8c-8f6d6a791bdf.png)

Figure 7B: Analyze the expression levels of biomarker snoRNAs in the tissues of normal controls and breast cancer patients using qRT-PCR 
(Note: samples are collected from each patient once, with at least 10 samples per group and repeated 3 times; this step can be performed 
if clinical tissue samples are available).

