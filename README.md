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

- Interesting finding from microarray data:
