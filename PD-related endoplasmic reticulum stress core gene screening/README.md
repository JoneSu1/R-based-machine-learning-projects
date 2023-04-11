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
