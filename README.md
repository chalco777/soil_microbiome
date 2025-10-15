
# Soil Microbiome Analysis Repository

This repository gathers multiple exploratory and inferential analyses performed on plasmid, viral, and chromosomal fractions recovered from soil metagenomes. Each folder contains self-contained scripts, intermediate tables, and figures that document the statistical pipelines used to quantify abundance patterns, run differential models, and visualize beta diversity.

## Repository Overview

```
.
├── functional_differential_analysis
│   ├── differential.Rmd
│   ├── results_1.0/
│   └── results_differential_2.0/
├── mosdepth_abundance
│   ├── bar_graphic.R
│   ├── mosdepth_plasmid_genes.sh
│   └── *.png
├── preliminar_analysis_eggnog
│   ├── preliminar_analysis.R
│   ├── *.gff / *.pdf / *.png
└── species_pcoa
    ├── pcoa.R
    ├── relative_abundances_per_sample.txt
    └── *.png
```

The following sections detail the purpose of each directory and script, highlighting the statistical methods that support the generated figures.

## `functional_differential_analysis`

### Purpose
This module performs multi-omics differential abundance modeling on EggNOG-derived annotations (COG, GO, KEGG pathways, BRITE hierarchies, and PFAM domains) using the [MaAsLin2](https://huttenhower.sph.harvard.edu/maaslin2/) linear modeling framework. Results are organized into snapshot folders (`results_1.0` and `results_differential_2.0`) that store coefficient tables, residual diagnostics, and curated plots.

### Key Script: `differential.Rmd`
- **Data wrangling**: Imports a combined GFF table of features, parses annotation columns, and aggregates read counts per functional category and barcode.【F:functional_differential_analysis/differential.Rmd†L33-L106】
- **Metadata handling**: Harmonizes samples into experimental stages (`baseline`, `24 hours`, `10 days`) to drive fixed-effect contrasts.【F:functional_differential_analysis/differential.Rmd†L106-L142】
- **API-backed annotation**: Queries KEGG, GO (including AmiGO fallbacks), and BRITE APIs to attach descriptive labels, which are later wrapped for readability in plots and heatmaps.【F:functional_differential_analysis/differential.Rmd†L144-L268】
- **Modeling pipeline**: Constructs CPM matrices, runs MaAsLin2 linear models (TSS normalization + log transform), and ranks features with Kruskal–Wallis tests complemented by pairwise Wilcoxon post hoc analyses.【F:functional_differential_analysis/differential.Rmd†L270-L360】
- **Visualization**: Generates faceted boxplots with statistical annotations and, when significant features are detected (q < 0.25), produces log2 fold-change and −log10(q) heatmaps for interpretability.【F:functional_differential_analysis/differential.Rmd†L360-L452】

### Outputs
- **MaAsLin2 folders** contain model fits (`fitted.rds`, `residuals.rds`) and logs for each annotation × origin combination.
- **PNG summaries** (e.g., [chromosomal pathway CPMs](functional_differential_analysis/results_differential_2.0/pathways_chromosome_cpm_1-9.png)) visualize the highest-ranking functional categories across timepoints.
- **Heatmaps** such as [significant KEGG pathways](functional_differential_analysis/results_differential_2.0/heatmap_significance_KEGG_pathways.png) translate regression statistics into color-coded effect summaries.

## `mosdepth_abundance`

### Purpose
Integrates read-mapping depth information with gene annotations to highlight the most abundant plasmid and viral products. The directory also stores finished barplots that summarize top functional hits.

### Key Script: `mosdepth_plasmid_genes.sh`
- Automates the end-to-end coverage annotation workflow: sorts and indexes BAM files, computes per-gene depth using `mosdepth`, and appends depth values back into the original GFF3 feature table.【F:mosdepth_abundance/mosdepth_plasmid_genes.sh†L1-L79】
- Logs each barcode’s processing steps and errors, ensuring reproducibility when rerunning across large cohorts.

### Key Script: `bar_graphic.R`
- Reads mosdepth-enriched TSVs for plasmids and viruses, filters out hypothetical proteins, and aggregates depth per functional product.【F:mosdepth_abundance/bar_graphic.R†L1-L40】
- Selects the top 10 products per barcode and visualizes them via faceted horizontal barplots, enabling rapid comparison of dominant functional signatures.【F:mosdepth_abundance/bar_graphic.R†L22-L52】【F:mosdepth_abundance/bar_graphic.R†L56-L75】

### Outputs
- [Top plasmid products](mosdepth_abundance/top10_products_plasmids.png) show barcode-specific gene depth rankings with consistent aesthetics.
- [Top viral products](mosdepth_abundance/top10_products_virus.png) mirror the plasmid plots for viral contigs.

## `preliminar_analysis_eggnog`

### Purpose
Provides exploratory statistics for EggNOG annotations prior to differential testing. It focuses on transforming depth tables into LEfSe-ready matrices and visualizations that inspect category-level variability.

### Key Script: `preliminar_analysis.R`
- Cleans GFF-derived tables, resolves ambiguous names, and aggregates depth per functional label across barcodes.【F:preliminar_analysis_eggnog/preliminar_analysis.R†L15-L47】
- Builds wide matrices with status metadata to feed downstream LEfSe discrimination analyses, ensuring only consistently observed features are retained.【F:preliminar_analysis_eggnog/preliminar_analysis.R†L48-L73】
- Constructs exploratory boxplots comparing depth distributions across timepoints for selected COG classes, with optional outlier removal (e.g., excluding `barcode03`).【F:preliminar_analysis_eggnog/preliminar_analysis.R†L75-L128】
- Identifies the nine most abundant COG categories (based on median depth) and renders faceted boxplots to profile shifts across experimental stages.【F:preliminar_analysis_eggnog/preliminar_analysis.R†L133-L181】

### Outputs
- [Faceted COG depth comparison](preliminar_analysis_eggnog/boxplot_all_categories_faceted.png) captures trends across the top functional classes.
- [Energy production and conversion boxplot](preliminar_analysis_eggnog/boxplot_energy_production_conversion.png) highlights a key metabolic pathway over time.

## `species_pcoa`

### Purpose
Assesses beta diversity among metagenomic samples via Principal Coordinates Analysis (PCoA) applied to species-level relative abundance tables.

### Key Script: `pcoa.R`
- Loads taxonomic abundance matrices, filters samples to the barcodes used in co-assemblies, and imposes prevalence/abundance thresholds to remove rare taxa.【F:species_pcoa/pcoa.R†L7-L56】
- Calculates Bray–Curtis dissimilarities, runs classical multidimensional scaling (`cmdscale`), and reports the variance explained by the first two axes.【F:species_pcoa/pcoa.R†L57-L95】
- Joins ordination scores with experimental metadata, applying `ggrepel` labeling for readability and exporting publication-ready figures.【F:species_pcoa/pcoa.R†L96-L133】
- Repeats the analysis on a broader cohort that includes biofilm barcodes, enabling comparison between standard and biofilm conditions.【F:species_pcoa/pcoa.R†L135-L200】

### Outputs
- [PCoA plot (selected barcodes)](species_pcoa/pcoa_analysis_selected_barcodes.png) visualizes temporal clustering among baseline, 24-hour, and 10-day samples.
- [PCoA plot with biofilm barcodes](species_pcoa/pcoa_analysis_selected_barcodes_filtered_sp.png) extends the ordination to include additional ecological states.

## Getting Started
While absolute paths in the scripts point to the original analyst’s workstation, reproducibility can be achieved by:
1. Cloning the repository and adjusting `setwd()` statements to match your local project root.
2. Ensuring dependencies are installed: R packages (`tidyverse`, `Maaslin2`, `vegan`, `ggpubr`, `ggrepel`, etc.) and command-line tools (`samtools`, `mosdepth`).
3. Re-running scripts within each folder to regenerate tables and figures as needed.

## Citation
If you use these workflows, please cite the relevant software (MaAsLin2, EggNOG-mapper, mosdepth, vegan) and acknowledge the project contributors.