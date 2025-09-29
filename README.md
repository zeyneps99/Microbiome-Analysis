# Investigating Non-Invasive Biomarkers and Risk Factors Associated with Liver Cirrhosis

This repository contains the R code and report for an MSc course project on liver cirrhosis, gut microbiome composition, and associated risk factors.

- 📄 [Full report (PDF)](results/Investigating Non-Invasive Biomarkers and Risk Factors Associated with Liver Cirrhosis.pdf)  
- 📝 [Report in Markdown](docs/liver_cirrhosis_report.md)

---

## Overview
Whole-metagenome shotgun (WMS) sequencing data were profiled with **MetaPhlAn4**, followed by statistical analyses in **R 4.3.1**.  
The goal: identify non-invasive biomarkers and risk factors for cirrhosis through diversity analyses, differential abundance testing, and covariate associations (BMI, albumin, MELD, HBV, etc.).

---

## Repository structure 
├── src/
│   ├── 1_install_packages.R        # installs required CRAN packages
│   └── 2_microbiome_analysis.R     # end-to-end analysis
├── data/                           # inputs (not tracked in git)
│   ├── metadata_parsed.csv
│   ├── profiled_metagenome.tab     # optional raw MetaPhlAn table
│   ├── species_abd.csv             # main table (species)
│   ├── genus_abd.csv               # optional for genus plots
│   ├── family_abd.csv              # optional
│   ├── order_abd.csv               # optional
│   ├── class_abd.csv               # optional
│   ├── phylum_abd.csv              # optional
│   └── kingdom_abd.csv             # optional
├── results/                     
│   ├── Investigating Non-Invasive Biomarkers and Risk Factors Associated with Liver Cirrhosis.pdf
│   └── figures/
└── README.md

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
