# Investigating Non-Invasive Biomarkers and Risk Factors Associated with Liver Cirrhosis

This repository contains the R code and report for an MSc course project on liver cirrhosis, gut microbiome composition, and associated risk factors.

- 📄 [Full report (PDF)](results/liver_cirrhosis_report.pdf)  
- 📝 [Report in Markdown](results/liver_cirrhosis_report.md)

---

## Overview
Whole-metagenome shotgun (WMS) sequencing data were profiled with **MetaPhlAn4**, followed by statistical analyses in **R 4.3.1**.  
The goal: identify non-invasive biomarkers and risk factors for cirrhosis through diversity analyses, differential abundance testing, and covariate associations (BMI, albumin, MELD, HBV, etc.).

---
## Repository structure

```text
.
├── R/
│   ├── 1_install_packages.R
│   └── 2_microbiome_analysis.R
├── data/                     # inputs (not tracked in git; keep .gitkeep)
│   ├── metadata_parsed.csv
│   ├── species_abd.csv
│   ├── genus_abd.csv
│   ├── family_abd.csv
│   ├── order_abd.csv
│   ├── class_abd.csv
│   ├── phylum_abd.csv
│   └── kingdom_abd.csv
├── results/
│   ├── tables/
│   ├── figures/
│   └── logs/
├── docs/
│   ├── Investigating Non-Invasive Biomarkers and Risk Factors Associated with Liver Cirrhosis.pdf
└── README.md
```

## Methods summary
Processes MetaPhlAn4 outputs and metadata to compute:
- Alpha diversity (Shannon, Simpson, Chao1) with Wilcoxon tests  
- Beta diversity (Bray–Curtis → PCoA with 95% ellipses)  
- Differential abundance at species and genus levels (Wilcoxon, mean Δ)  
- Covariate correlations (Spearman) with MELD, BMI, Alb, Crea, TB, HBV, age, sex

## Outputs
	•	Alpha diversity boxplots (Simpson, Shannon, Chao1)
	•	PCoA (Bray–Curtis) with 95% ellipses
	•	Stacked barplots (genus/species mean relative abundance)
	•	Top Δ mean species barplot
	•	Key-species boxplots
	•	Covariate scatter plots (diversity vs factors; MELD correlations)
	•	results/tables/
	•	Diversity metrics per sample
	•	Differential abundance results (p-values, Δ means)

## Citation
Primary dataset: Qin N. et al., Nature 2014.
