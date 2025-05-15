# Differential Gene Expression Tutorial

This repository provides a step-by-step R-based tutorial for analyzing microarray gene expression data. It covers essential preprocessing, exploratory data analysis (EDA), and visualization techniques, culminating in differential expression analysis. The tutorial is designed to be accessible to both beginners and intermediate users in bioinformatics.

---

## Table of Contents

1. [Overview](#overview)
2. [Prerequisites](#prerequisites)
3. [Repository Structure](#repository-structure)
4. [Tutorial Workflow](#tutorial-workflow)
5. [Data Sources](#data-sources)
6. [Usage Instructions](#usage-instructions)
7. [Contributing](#contributing)
8. [License](#license)
9. [Contact](#contact)

---

## Overview

^[This tutorial guides users through the analysis of microarray gene expression data, focusing on the following objectives:]({"attribution":{"attributableIndex":"735-11"}})

- **Data Loading**: ^[Importing raw gene expression data into R.]({"attribution":{"attributableIndex":"1216-1"}})
- **Exploratory Data Analysis (EDA)**: ^[Performing principal component analysis (PCA) and identifying top variable genes.]({"attribution":{"attributableIndex":"1216-3"}})
- **Probe-to-Gene Mapping**: ^[Mapping probes to gene symbols.]({"attribution":{"attributableIndex":"1216-5"}})
- **Visualization**: ^[Generating interactive volcano plots and heatmaps.]({"attribution":{"attributableIndex":"1216-7"}})
- **Differential Expression Analysis**: ^[Identifying differentially expressed genes.]({"attribution":{"attributableIndex":"1216-9"}}) [oai_citation:0‡GitHub](https://github.com/bioinfosourabh/Differential-Gene-expression-Analysis/?utm_source=chatgpt.com)

^[The tutorial emphasizes reproducibility and clarity, making it suitable for educational purposes and practical applications in gene expression analysis.]({"attribution":{"attributableIndex":"1629-0"}})

---

## Prerequisites

To run the scripts in this repository, ensure you have the following R packages installed:


```r
install.packages(c("ggplot2", "pheatmap", "EnhancedVolcano"))
BiocManager::install(c("limma", "affy"))
```

# Tutorial Workflow

1. Data Loading

Begin by downloading and loading the raw gene expression data into R using the provided shell script:

```bash
scripts/01_download_data.sh
```
# 2. Exploratory Data Analysis (EDA)

Perform initial data exploration to understand the structure and quality of the data:
	•	Raw Data EDA: Use scripts/02_EDA_for_raw.R to examine the raw data distribution and detect potential issues.
	•	Preprocessed Data EDA: Use scripts/03_EDA_for_preprocessed.R to assess the data after preprocessing steps.

# 3. Differential Expression Analysis

Conduct differential expression analysis using the limma package:
	•	Run scripts/04_Differential_Expression_Analysis.R to identify genes that are differentially expressed between conditions.
	•	The script generates results including log-fold changes, p-values, and adjusted p-values.

# 4. Visualization

Enhance the interpretability of results through visualization:
	•	Volcano Plot: Visualize the significance and magnitude of differential expression using the EnhancedVolcano package.
	•	Heatmap: Generate heatmaps of top differentially expressed genes using the pheatmap package.

These visualizations aid in the identification of key genes and patterns in the data.
