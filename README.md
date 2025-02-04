# R script

## project "GO Enrichment analysis TB mice"
This project is a personal entertaiment and learning project.

## Description
This project uses published molecular data to explore the relative abundance of biological processes in mice infected by $\textit{Mycobacterium tuberculosis}$ with three groups of doagnosible conditions. 
This project contains the R script for replicating the analysis and the RNoteBook for running the analysis by step. 
The raw data is not provided in this project, but it is downloaded from GEO as part of the script. 
A rendered html document is also provided for education purposes.

## License
This repository is licensed under the MIT License.

This script/notebook is free to use; you can redistribute it and/or modify it under the terms of the GNU General Public License,version 3, as published by the Free Software Foundation.
This notebook is distributed in the hope that it will serve educational and entertainment purposes, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose. 
Use of this code may be subjected to licensing agreement and copyright requirement of the Journal, Infection and Immunity.

For more details of the license, please refer to: www.gnu.org/licenses/gpl-3.0.en.html.

## Context
In the source article, Diversity Outbred mice were infected with $\textit{Mycobacterium tuberculosis}$ (TB), in order to identify feature of asymptomatic TB infection using computational and statistical approach.
Methods of analysis range are summarised in the following table:

| **Source of Evidence**                                        | **Method of Analysis**                                                                                                                 | **Connection to Other Parts of the Article**                                                                                                                                                             |
| ------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| <span style="color:grey;">**Lung tissues**                    | <span style="color:grey;">- Histopathology (H&E staining, carbol fuschin staining) - Immunohistochemistry (CD20 staining for B-cells)  | <span style="color:grey;">- Validates transcriptional findings by localizing CD20+ B-cells in perivascular and peribronchiolar regions. - Links tissue-specific cellular features to imaging biomarkers. |
| **Gene expression microarray**                                | - RNA extracted from lung tissue, hybridized to Mouse Gene (v2.0) ST microarrays - Pathway analysis  using Enrichr                     | - Identifies transcriptional correlates of resistance and susceptibility (e.g., B-cell activation, receptor signaling). - Connects to immunohistochemistry findings.                                     |
| <span style="color:grey;">**Cytokines and chemokines**        | <span style="color:grey;">- ELISA on lung homogenates for markers like CXCL1, TNF, IL-10, and MMP8                                     | <span style="color:grey;">- Provides protein-level validation of inflammatory and immune responses identified in transcriptional data.                                                                   |
| <span style="color:grey;">**Anti-M. tuberculosis antibodies** | <span style="color:grey;">- ELISA to detect IgG specific to M. tuberculosis components                                                 | <span style="color:grey;">- Adds specificity to classification models and complements cytokine data in differentiating disease states.                                                                   |
| <span style="color:grey;">**Deep-learning model**             | <span style="color:grey;">- Attention-based multiple instance learning applied to lung granuloma images                                | <span style="color:grey;">- Pinpoints imaging biomarkers (e.g., perivascular and peribronchiolar lymphocytic cuffs) connected to histological and transcriptional findings.                              |
| <span style="color:grey;">**Mouse survival data**             | <span style="color:grey;">- Statistical analysis of survival times and disease progression (log-rank test, ANOVA)                      | <span style="color:grey;">- Provides the phenotypic framework for classifying disease states (progressors, controllers, asymptomatic).                                                                   |
| <span style="color:grey;">**Quantitative imaging**            | <span style="color:grey;">- Density and distribution analysis of CD20+ B-cells using entropy-based quantification and machine learning | <span style="color:grey;">- Correlates spatial localization of B-cells to resistance mechanisms and imaging biomarkers identified by deep learning.                                                      |

## Interesting finding from microarray data:
(1) 105 genes related identified as overrepresented, and 5 pathways are found associated with B-cell pathways.
(2) Non-parametric Mann-Whitney U-test is chosen to find the differentially expressed genes (DEGs), and the parametric Welch's t-test is computed for comparison in the Supplementary Information.

## Potential issues
 (1) **Comparison method**: 105 overrepresented genes are not recovered by pairwise comparison. Instead, asymptomatic and the combination of chronic controller, progressor, and non-infected mice are compared. This may:
   1)Introduce noise in the data (e.g. masking distinct GO terms in chronic controller, progressor, and non-infected mice and distinct GO terms between asymptomatic and any one of the remaining groups)
	 2)Data from chronic controller, progressor, and non-infected may not be homogenous enough.
 (2) **GEO data**: A subset of microarray data (i.e. [GSE179417]()) is said to be published elsewhere, which likely refers to the dataset under a different article. The complete microarray data is not linked in the article. The article name locates [GSE266564](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266564) as the microarray data for this article. 
 (3)**GEO annotation platform**:  [GSE266564](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266564) uses platform [GPL32068](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL32068) for annotating the microarray data. [GPL32068](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL32068) is described as identical to [GPL16570](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL16570) but the table dimension attributes and number of rows are not the same. 
 (4) **R script**: The R script to perform overrepresentation analysis is not provided in the article (i.e. no file in Supplementary Information or GitHub link).
 (5) **DEG computation**: The choice of Mann-Whitney U-test over Welch's t-test and Bayesian moderation is seemly not specified in the article, though we may presume the choice is based on the lack of assumption on normality of the underlying data.
 (6) **Enrichr**: Enrichr does not preserve all parent-child relationship among their hierarchy, and a cut-off level is implemented. By comparison, enrichGO does preserve all hierarchy. The two implementations have been found to recovered different number of genes in overrepresentation analysis ([Yu Lab](https://mp.weixin.qq.com/s/6lSsg2WMEK2btwve-9C2rA)). 

## Task
 (1) Perform Exploratory Data Analysis on the source CEL files, to test whether the choice of Mann-Whitney U-test over Bayesian moderation is valid based on two criteria: Normality and outliers. 
 (2) Perform an overrepresentation analysis on Biological Process using [GSE266564](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266564) and [GPL16570](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL16570), which is an exploration to recreate the result from the article based on what may be considered as complete and equivalent data.

 ## Conclusion
 (1) The source data is non-normal with rare outlier, thus Mann-Whitney U-test is indeed preferred over Bayesian moderation or Welch's t-test. However, sample sizes are imbalanced among group, which may affect Mann-Whitney U-test. Brunner-Munzel test may be tried in the future.
 (2) Of the 11 overrepresented genes reported in the asymptomatic group in the source article, 4 are recovered from this project, namely RAE1, BANK1, NUP107, and CD79B.
 (3) Of the 11 overrepresented genes reported in the source article, asymptomatic controller do not show any of the 11 genes when compared to the chronic controller group.
 (4) Many pathways are found to be overrepresented between asymptomatic controller show 249 genes with a wide range of pathways, instead of strict B-cell pathways.
