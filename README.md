# Bioinformatic Research Project

## Packages used

## Stakeholder relevance

## project "Gene Ontology, Pathway Enrichment analysis on TB infected mice"
This project falls within the realm of biomedical and bioinformatics research using publicly available data of mice model species. <br>
This project is a project for testing reproducibility in the literature and education in bioinformatics. <br>

## Description
This project uses published microarray molecular data to explore the relative abundance of genes and biological processes in mice infected by *Mycobacterium tuberculosis* with three discrete categories of clinical conditions ([source article](https://pubmed.ncbi.nlm.nih.gov/38899881/)). The goal is finding significant differences in genes that are related to various stages of infections (e.g. asymptomatic and chronic infections). <br>
This project contains the R script for replicating the analysis. <br>
An RNotebook is also provided to view partial results througout the analysis without the need to fully run the script. The RNotebook also contains additional comments througout the analysis for teaching and learning purposes. <br>
For non-R users, a rendered html is provided for viewing the code and result of the RNotebook, but the html file is decently large, and will need to be downloaded to view locally. <br>
The raw data is not provided in this project, but it will be downloaded from GEO when script is run. <br>

## License
This repository is licensed under the MIT License.

This script/notebook/rendered html are all free to use; you can redistribute them and/or modify them under the terms of the General Public License,version 3 (GNU), as published by the Free Software Foundation.
This notebook is distributed in the hope that it will serve educational purposes, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose. 
Use of this code may be subjected to licensing agreement and copyright requirement of the Journal, [Infection and Immunity](https://journals.asm.org/).

For more detail information of the MIT and GNU licenses, please refer to [MITlicense website](https://opensource.org/license/mit) and [GNU general public license](www.gnu.org/licenses/gpl-3.0.en.html).

## Context
In the source article, Diversity Outbred mice were infected with *Mycobacterium tuberculosis* (TB), to identify features of asymptomatic TB infection using experimental and computational approaches.
Methods of analysis used by the source article are summarised in the following table:

| **Source of Evidence**                                        | **Method of Analysis**                                                                                                                 | **Connection to Other Parts of the Article**                                                                                                                                                             |
| ------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| <span style="color:grey;">**Lung tissues**                    | <span style="color:grey;">- Histopathology (H&E staining, carbol fuschin staining) <br> - Immunohistochemistry (CD20 staining for B-cells)  | <span style="color:grey;">- Validates transcriptional findings by localizing CD20+ B-cells in perivascular and peribronchiolar regions. <br> - Links tissue-specific cellular features to imaging biomarkers. |
| **Gene expression microarray**                                | - RNA extracted from lung tissue, hybridised to Mouse Gene (v2.0) ST microarrays <br> - Pathway analysis  using Enrichr                     | - Identifies transcriptional correlates of resistance and susceptibility (e.g., B-cell activation, receptor signaling). <br> - Connects to immunohistochemistry findings.                                     |
| <span style="color:grey;">**Cytokines and chemokines**        | <span style="color:grey;">- ELISA on lung homogenates for markers like CXCL1, TNF, IL-10, and MMP8                                     | <span style="color:grey;">- Provides protein-level validation of inflammatory and immune responses identified in transcriptional data.                                                                   |
| <span style="color:grey;">**Anti-*M. tuberculosis* antibodies** | <span style="color:grey;">- ELISA to detect IgG specific to *M. tuberculosis* components                                                 | <span style="color:grey;">- Adds specificity to classification models and complements cytokine data in differentiating disease states.                                                                   |
| <span style="color:grey;">**Deep-learning model**             | <span style="color:grey;">- Attention-based multiple instance learning applied to lung granuloma images                                | <span style="color:grey;">- Imaging biomarkers (e.g., perivascular and peribronchiolar lymphocytic cuffs) connected to histological and transcriptional findings.                              |
| <span style="color:grey;">**Mouse survival data**             | <span style="color:grey;">- Statistical analysis of survival times and disease progression (log-rank test, ANOVA)                      | <span style="color:grey;">- Provides the phenotypic framework for classifying disease states (progressors, controllers, asymptomatic).                                                                   |
| <span style="color:grey;">**Quantitative imaging**            | <span style="color:grey;">- Density and distribution analysis of CD20+ B-cells using entropy-based quantification and machine learning | <span style="color:grey;">- Correlates spatial localization of B-cells to resistance mechanisms and imaging biomarkers identified by deep learning.                                                      |

## Interesting finding from microarray data:
1. 105 genes are identified as overrepresented/enriched, and 5 pathways are found associated with B-cell pathways. <br>
2. Non-parametric Mann-Whitney U-test is chosen to find the differentially expressed genes (DEGs), and the parametric Welch's t-test is computed for comparison in the Supplementary Information. <br>

## Potential issues
 1. **Comparison methods**: 105 overrepresented genes are not recovered by pairwise comparison. Instead, asymptomatic and the combination of chronic controller, progressor, and non-infected mice are compared. This may cause negative impacts such as:<br>
   (1) Introduce noise in the data (e.g. masking distinct GO terms in chronic controller, progressor, and non-infected mice) <br>
   (2) Data from chronic controller, progressor, and non-infected may not be homogenous enough. <br>
 2. **GEO data**: A subset of microarray data (i.e. GSE179417) is said to be published elsewhere, which likely refers to the dataset under a different article. The complete microarray data is not linked in the article. The article name locates [GSE266564](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266564) as the microarray data for this article. <br>
 3. **GEO annotation platform**:  [GSE266564](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266564) uses platform [GPL32068](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL32068) for annotating the microarray data. [GPL32068](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL32068) is described as identical to [GPL16570](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL16570) but the table dimension/attributes and number of rows are not the same. <br>
 4. **R script**: The R script to perform overrepresentation/enrichment analysis is not provided in the article (i.e. no file in Supplementary Information or GitHub link). <br>
 5. **DEG computation**: The choice of Mann-Whitney U-test over Welch's t-test and Bayesian moderation is seemly not specified in the article, though we may presume the choice is based on the lack of assumption on normality of the underlying data. <br>
 6. **Enrichr**: Enrichr does not preserve all parent-child relationship among their hierarchy during computation of the gene ontology terms, and a cut-off level is implemented implicitly. By comparison, enrichGO does preserve all hierarchy. The two implementations have been found to recovered different number of genes in overrepresentation analysis ([Yu Lab](https://mp.weixin.qq.com/s/6lSsg2WMEK2btwve-9C2rA)). Thus, the use of different computational methods may lead to different results of overrepresentation. <br>

## Objectives
 1. Perform Exploratory Data Analysis on the source CEL files, to test whether the choice of Mann-Whitney U-test over other methods such as Bayesian moderation and Welch's t-test is valid based on two criteria: Normality and outliers. <br>
 2. Use EnrichGO() from the package ClusterProfiler and manual computation to perform an overrepresentation analysis on Biological Process using [GSE266564](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266564) and [GPL16570](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL16570), which is an exploration to recreate the results from the article based on what is considered as complete and equivalent data. <br>

 ## Conclusion
 1. The raw CEL files show a high level of correlations even before normalisation, which may suggest the variability between groups is relatively low. By using rma() in this analysis, source data might have been over normalisaion, causing FoldChange computation issues (i.e. identical in value). A better technique with higher sensitivity may be necessary for this source data. <br>
 ![load microarray data-1](https://github.com/user-attachments/assets/cddd6b7f-80ee-47c4-b8be-78f386309d9b)
 2. The source data is non-normal with rare outlier, thus Mann-Whitney U-test is indeed preferred over Bayesian moderation or Welch's t-test. However, sample sizes are imbalanced among group, which may affect Mann-Whitney U-test. Brunner-Munzel test or something similar may be worth trying in the future. <br>
 3. Some overrepresented genes (e.g., Bank1, Rae1) lack annotations on [GPL16570](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL16570), which might be one of the reasons why the source article uses a custom annotation platform [GPL32068](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL32068). <br>
 4. Of the 11 overrepresented genes reported in the asymptomatic group in the source article, 4 are recovered from this project, namely RAE1, BANK1, NUP107, and CD79B. <br>
 5. Of the 11 overrepresented genes reported in the source article, asymptomatic controller does not show any of the 11 genes when compared to chronic controller. <br>
 6. Many pathways are found to be overrepresented between asymptomatic controller and other experimental groups, showing 249 genes with a wide range of pathways, instead of only B-cell pathways. <br>
