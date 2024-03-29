# Code supporting Tommasini, D, Fox, R, Ngo, KJ, Hinman, JD, Fogel, BL (2023). Alterations in oligodendrocyte transcriptional networks reveal region-specific vulnerabilities to neurological disease. iScience, 26, 4:106358. 
URL: https://doi.org/10.1016/j.isci.2023.106358

Contents: 

1. RNAseq_parallel_paired_reverse.sh	(Read alignment and quantification of raw sequence files from Gene Expression Omnibus)
2. correlation.R	(Oligodendrocyte gene expression correlation across regions)
3. DESeq2.R (Differential gene expression using negative binomial distribution)
4. hierarchical_clustering.R (Hierarchical clustering analysis)
5. PCA_analysis.R (Principal component analysis)
6. run_salmon.sh (Quantification of OL transcripts at isoform resolution for downstream differential splicing analysis in IsoformSwitchAnalyzeR)
7. WGCNA_2trait.Rmd (Example weighted gene co-expression network analysis)
