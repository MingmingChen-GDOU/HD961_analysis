# HD961 Salt-Stress Translatomic and Genomic Analysis

This repository contains all scripts and resources needed to reproduce the genome assembly, RNA‑seq, and Ribo‑seq analyses presented in "Genomic and translational insights into eIF2B-mediated salt tolerance in sea rice HD961" (Chen et al.).

---

## Table of Contents

1. [Overview](#overview)
2. [Directory Structure](#directory-structure)
3. [Dependencies](#dependencies)
4. [Usage](#usage)

   * [Raw Data Processing](#raw-data-processing)
   * [Genome Assembly & Comparison](#genome-assembly--comparison)
   * [DESeq2 RNA‑seq Analysis](#deseq2-rna-seq-analysis)
   * [Ribo‑seq Analysis](#ribo-seq-analysis)
   * [Translation Efficiency & Codon Usage](#translation-efficiency--codon-usage)
   * [Enrichment & Plotting](#enrichment--plotting)
5. [Outputs](#outputs)
6. [Contact](#contact)

---

## Overview

Sea rice HD961 exhibits exceptional salt tolerance. This repository includes scripts to:

* Assemble and scaffold the HD961 genome (Nanopore + Hi‑C).
* Compare assembly metrics against other published rice genomes.
* Process and quality‑control RNA‑seq and Ribo‑seq raw reads.
* Perform differential expression (DESeq2) on RNA‑seq data.
* Identify differential translation (edgeR/DESeq2) on Ribo‑seq data and calculate TE.
* Analyze codon usage and ribosome stalling (GCG occupancy).
* Conduct GO and KEGG enrichment (clusterProfiler).
* Generate publication‑quality figures (volcano plots, heatmaps, bar/dot‑plots, metagene profiles).

---

## Directory Structure

```
├── README.md
├── rawdata_processing.txt      # Commands for trimming, filtering, mapping raw reads
├── genome_assembly.txt         # Pipeline for Nanopore assembly & Hi‑C scaffolding
├── genome_comparison.R         # Compare assembly stats across rice genomes
├── Deseq.R                     # RNA‑seq differential expression analysis
├── riboseq.R                   # Ribo‑seq DE & translation efficiency analysis
├── GCC of genes.R              # Codon usage & GCG occupancy plotting
├── clusterProfiler.R           # GO/KEGG enrichment and plotting scripts
├── scatter_highlight.R         # Volcano & scatter highlighting routines
├── data/                       # (Optional) CSV inputs and outputs
├── figures/                    # Exported figure files
└── tables/                     # CSV tables for manuscript
```

---

## Dependencies

* R (>= 4.0)
* R packages:

  * `data.table`, `dplyr`, `ggplot2`, `ggrepel`, `pheatmap`, `edgeR`, `DESeq2`, `limma`, `GenomicRanges`, `RColorBrewer`, `gprofiler2`, `clusterProfiler`, `org.Osativa.eg.db`
* Python (optional for auxiliary scripts)
* Bioinformatics tools:

  * FastQC, fastp, HISAT2, Bowtie2, Samtools, StringTie, and others (see `rawdata_processing.txt`).

## Contact

For questions or issues, please open a GitHub Issue or contact:

Mingming Chen ([mingming.chen@gdou.edu.cn](mailto:mingming.chen@gdou.edu.cn))

---

*Last updated: \$(date +"%Y-%m-%d")*
