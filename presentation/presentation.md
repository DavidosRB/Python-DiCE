---
title: Introduction
description: Marp Presentation Example
author: David Buchmann
header: David Buchmann | Knowledge Graphs and Algorithms Lab
footer: 06.03.2026
lang: en
size: 16:9
paginate: true
theme: default
math: mathjax
marp: true
---

# Knowledge Graphs & Algorithms

#### Presentation of Results by David Buchmann

---

# Data

- Chosen Dataset: [GSE239282](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE239282)
- Topic: Age-related cognitive disorders (Alzheimers, etc.) - Treatment with Music 
- Total of 60 samples, before and after treatment
  - 16 disease (ACD), 14 control samples used (timepoint 2)
- Raw data:
  - 38327 unique gene rows across ACD & Control
- GEO2R:
  - 16951 genes with analyzed expression
  - 3091 genes after NaN filtering & HGNC mapping

---

# Methods

- Download of Raw Data & GEO2R Analysis
- Calculation of (Mutual) Information Gain & above average filtering
- Setup of PPI network based on StringDB based on our proteins
- Mapped pearson's correlation as edge distance
- Calculated Betweenness and Eigenvector Centrality
  - EC once using Distance (1-PC), once using PC
  - Only kept above average EC & BC genes
- Ranked (normalized) EC & BC - combined into ensemble ranking
- Gene Set Enrichment Analysis via GSEAPY on KEGG & GO pathways

---

# Results
### Data Filtering

- 3041 raw data gene rows after NaN filtering & HGNC mapping
- 3091 GEO2R genes after NaN filtering & HGNC mapping

---

# Results
### Top 5 Gene Results

---

# Results
### KEGG Pathway Analysis

---

# Results
### KEGG Interpretation

---

# Results
### GO Pathway analysis

---

# Results
### GO Interpretation

---

# Results
### Plots

---
