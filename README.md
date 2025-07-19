# 🧬 lncRNA Biomarker Identification in Hepatocellular Carcinoma (HCC)

This repository contains the R scripts and workflow used for identifying **lncRNA biomarkers in HCC** through **differential expression analysis (DEA)** and **clustering-based annotation**.


---
The project was completed as part of my **Final Year Project (FYP)** in the BS Bioinformatics program.
3 Datasets were used in my project, where slight changes were made accoridng to the accession ids. This repository contains codes modified for GSE72170.


## 📌 Objective

To identify potential **lncRNA biomarkers** in **Hepatocellular Carcinoma (HCC)** by:

- Performing differential gene expression analysis on GEO microarray datasets
- Annotating probe sequences using BLAST and GENCODE reference files
- Filtering for lncRNAs through genomic feature mapping
- Applying unsupervised clustering on lncRNA expression profiles using:
  - 🔹 **K-Means**
  - 🔹 **Gaussian Mixture Models (GMM)**
  - 🔹 **Self-Organizing Maps (SOM)**
- Evaluating clustering models using:
  - 📈 **Silhouette Score**
  - 📊 **Partition Coefficient (PC)**
  - 📉 **Partition Entropy (PE)**
---

## 🧪 Workflow Summary

### 🔹 Step 1: Differential Expression Analysis (DEA)

- Dataset retrieved from **NCBI GEO**
- Tools used:
  - `GEOquery` used for retrieving microarray expression data
  - `limma` used for identifying differentially expressed genes (DEGs)

### 🔹 Step 2: Annotation and lncRNA Identification

- Probe sequences corresponding to DEGs were added to the DEG table
- **BLAST** was used to align probe sequences against the **human reference genome (GRCh38)** from GENCODE
- BLAST results were mapped to a GENCODE **annotation file (GFF3 format)**
- Genes overlapping with annotated lncRNAs were filtered and selected

### 🔹 Step 3: Clustering Analysis

- Expression profiles of filtered lncRNAs were subjected to clustering using **Google Colab**
- Clustering methods applied:
  - K-Means
  - Gaussian Mixture Models (GMM)
  - Self-Organizing Maps (SOM)
- Evaluation metrics:
  - Silhouette Score
  - Partition Coefficient (PC)
  - Partition Entropy (PE)

## 📦 R Packages Used

```
# Differential expression
library(GEOquery)
library(limma)

# Annotation and lncRNA extraction
library(GenomicRanges)
library(rtracklayer)
library(dplyr)

```

## 📦 Python Packages Used

```
# K-means
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score
from google.colab import files

# GMM
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score, davies_bouldin_score, calinski_harabasz_score
from google.colab import files


