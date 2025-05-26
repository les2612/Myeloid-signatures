# Diagnostic Myeloid Signatures in Peripheral Blood of Autoimmune Diseases
> **Ekaterina Lebedeva**</br>
> tg: *@katya_lebe* </br>
> email: *EkaterinaLebedeva2612@gmail.com*

This repository contains the results of the bioinformatics analysis for my scientific thesis project (Bioinformatics Institute, 2024–2025).

**Supervisor**: Alexander Zaytsev
---

## Motivation
Autoimmune diseases and immunotherapy-induced autoimmunity represent major clinical challenges due to their unpredictable onset and complex immune involvement [1]. Current diagnostics often lack sensitivity and specificity, especially for early detection. Gene expression profiling of peripheral blood offers a promising, non-invasive approach to monitor immune activity and identify at-risk patients before symptoms develop [2, 3].

### Aim
To identify robust, clinically meaningful myeloid gene signatures in peripheral blood that are associated with autoimmune activity across diverse conditions and validate their utility in an immunotherapy context.

### Objectives
1. Perform differential expression analysis across 12 autoimmune diseases.

2. Select genes with strong discriminative power (ROC AUC > 0.55).

3. Derive and refine myeloid gene signatures using correlation-based clustering.

4. Functionally annotate signatures via GO enrichment analysis.

5. Validate signature activity in an independent cancer cohort with/without autoimmune complications.

6. Apply immune deconvolution to control for cellular heterogeneity.


## Raw data
This study leveraged publicly available bulk RNA-seq datasets encompassing over 1,000 whole blood samples from both healthy individuals and patients diagnosed with 12 distinct autoimmune diseases. The datasets were retrieved from the NCBI Gene Expression Omnibus (GEO) and include:
| Diagnosis                                  | Datasets                                                                 | Samples |
|--------------------------------------------|--------------------------------------------------------------------------|---------|
| Healthy                                    | GSE117769, GSE112087, GSE123658, GSE159225, GSE112057, GSE120178,        |
|                                            | GSE90081, GSE176260, GSE183701, GSE72509, GSE147339, GSE110685           | 292     |
| Systemic lupus erythematosus               | GSE72509, GSE110685, GSE112087                                           | 196     |
| Rheumatoid arthritis                       | GSE120178, GSE117769, GSE90081                                           | 193     |
| Juvenile idiopathic arthritis              | GSE112057                                                                | 115     |
| Type 1 diabetes mellitus                   | GSE123658, GSE183701                                                     | 65      |
| Crohn disease                              | GSE112057                                                                | 60      |
| Relapsing-remitting multiple sclerosis     | GSE159225                                                                | 20      |
| Systemic inflammatory response syndrome    | GSE176260                                                                | 17      |
| Ulcerative colitis                         | GSE112057                                                                | 15      |
| Arthropathic psoriasis                     | GSE117769                                                                | 10      |
| Psoriasis                                  | GSE147339                                                                | 10      |
| Secondary progressive multiple sclerosis   | GSE159225                                                                | 10      |
| Ankylosing spondylitis                     | GSE117769                                                                | 6       |
| **Total**                                  |                                                                          | **1009** |


## Content
You can find all the results in the `code` folder:

- [*1_getting_expressions.ipynb*] - This notebook loads expression data directly from Amazon S3 using a custom utility module. It generates a raw TPM expression matrix `data/expression_all.csv` for the selected diagnosis, along with the corresponding annotation file `data/annotations_filtered.csv`.
- [*2_expressions_filtration.ipynb*] - This notebook performs exploratory analysis and filtering of gene expression matrices from multiple datasets. It focuses on assessing data quality and generates data/annotations_whole_blood.csv and data/expression_whole_blood_all_genes.csv, which are ready for downstream analyses (e.g., differential expression).
- [*3_diffexpression.ipynb*] - This notebook performs differential gene expression (DE) analysis using the [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) package from [R](https://www.r-project.org/), integrated into Python via [rpy2](https://rpy2.github.io/). It generates `data/final_dif_gene_list.txt`, a list of genes differentially expressed between healthy individuals and autoimmune patients.. It generates `data/final_dif_gene_list.txt`, a list of genes differentially expressed between healthy individuals and autoimmune patients. 


**Note!** R and DESeq2 Setup required for DE analysis

---
                        1. Install R (if not already):
                        - Ubuntu/Debian: `sudo apt install r-base`
                        - macOS: `brew install r`

                        2. Open R and run:
                        ```r
                        if (!requireNamespace("BiocManager", quietly = TRUE))
                            install.packages("BiocManager")
                        BiocManager::install("DESeq2")
                        ```
----
- [*4_genes_filtration_by_ROC_AUC.ipynb*] - This notebook calculates and generates a clustermap of ROC AUC scores for differentially expressed genes (DEGs) across multiple datasets. It also generates `data/stable_up_genes.txt` — a list of stable genes with consistently elevated expression (ROC AUC > 0.5) in autoimmune patients.
- [*5_getting_primary_signatures.ipynb*] - This notebook performs correlation analysis of a selected `stable_up_genes.txt`list across multiple datasets. It builds individual gene-gene correlation clustermaps, then computes a consensus hierarchical clustermap. Stable correlation clusters are identified and extracted as primary gene signatures.
- [*6_getting_final_signatures.ipynb*] - This notebook performs correlation analysis of a selected `stable_up_genes.txt`list across multiple datasets. It builds individual gene-gene correlation clustermaps, then computes a consensus hierarchical clustermap. Stable correlation clusters are identified and extracted as primary gene signatures.
- [*7_GO_analysis.ipynb*] - - This notebook performs Gene Ontology (GO) analysis of final signatures
- [*8_validation.ipynb*]  - This notebook validates gene signatures by computing ssGSEA scores in the GSE28754 dataset, analyzing their expression dynamics over time, and performing immune cell deconvolution using [CIBERSORT](https://github.com/Moonerss/CIBERSORT.git).
CIBERSORT was run with default parameters using the standard LM22 signature matrix and a prepared input file `data/tpm_matrix_for_cibersort.txt` containing TPM-normalized gene expression data.

## Referenses

1. Pisetsky D. S. (2023). Pathogenesis of autoimmune disease. Nature reviews. Nephrology, 19(8), 509–524. https://doi.org/10.1038/s41581-023-00720-1

2. Nagafuchi, Y., Yanaoka, H., & Fujio, K. (2022). Lessons From Transcriptome Analysis of Autoimmune Diseases. Frontiers in immunology, 13, 857269. https://doi.org/10.3389/fimmu.2022.857269

3. Lozano-Rabella M, et al. Gene expression profiles in peripheral blood predict autoimmune toxicities induced by immune checkpoint inhibitors. Sci Rep. 2020;10:7749. https://doi.org/10.1038/s41598-020-63757-3



