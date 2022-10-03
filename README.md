# MSc_Thesis_Transcriptomics

Repository for all R scripts and results associated with my research dissertation for the MSc Bioinformatics and Computational Biology degree in University College Cork.

#### Project title:
"Transcriptomic analysis of the substantia nigra reveals alterations in key biological systems involved in the pathogenesis of Parkinson’s Disease"

#### Main objectives:

1.	Merge SN gene expression data generated from independent microarray studies for a statistically high-powered transcriptomic analysis of the SN in PD cases and healthy controls. 
2.	Use the merged microarray dataset to perform differential gene expression and functional enrichment analyses of SN tissue in PD cases and healthy controls. 
3.	Perform differential gene expression analysis of individual laser-captured dopaminergic neurons of the substantia nigra pars compacta and compare to the results obtained from the whole SN.
4.	Use the merged microarray dataset to perform weighted gene co-expression network analysis of SN tissue in PD cases and healthy controls, for the detection of biologically relevant differentially co-expressed gene sets. 
5.	Achieve a comprehensive tissue-wide insight into the molecular alterations and mechanisms contributing to PD at the transcriptome-level of the SN. 


#### Datasets: 
Gene expression data used in this study was downloaded from the Gene Expression Omnibus (GEO).

Merging - GSE8397, GSE20292, GSE49036, GSE7621, GSE20186.

Laser Captured - GSE20141


#### Scripts folder:
DE_GO_GSEA_workflow.R - Workflow used for data pre-processing, merging, differential expression analysis, GO over-representation analysis and GSEA.

process_meta.R - Manual assembly and processing of the merged sample metadata.

Laser_dissected.R - Differential expression nalysis of Laser Captured Microdissection (LCM) Dataset - GSE20141

WGCNA_workflow.R - Workflow used to perform WGCNA and module profiling.

Helper_functions.R - All functions utlized in this study.

#### Results folder:


