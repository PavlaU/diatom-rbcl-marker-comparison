# Evaluation of two overlapping rbcL amplicons for routine diatom DNA biomonitoring

This repository contains R scripts used to process and analyse data associated with the accepted manuscript *Evaluation of two overlapping rbcL amplicons for routine diatom DNA biomonitoring* (Urbánková & Štěpka, in press).

The study compares two partially overlapping diatom **rbcL** metabarcoding markers, 263 bp and 331 bp, using the same benthic diatom samples. Raw sequencing data are available from the European Nucleotide Archive under accession PRJEB106074. The analyses cover ASV-level comparison, taxonomic assignment, community composition and IPS-based ecological assessment.

## Repository structure

```text
.
├── scripts/
│   ├── 01_DADA_pipeline.r
│   ├── 02_technical_replicate_consistency_FigS1.R
│   ├── 03_rarefaction_procrustes_FigS2.R
│   ├── 04_ASV_correlation_Fig1.R
│   ├── 05_prevalence_abundance_FigS3.R
│   ├── 06_primer_mismatch_FigS4.R
│   ├── 07_species_taxonomy_community_FigS5.R
│   └── 08_IPS_analysis_Fig2_Fig3.R
└── input/
    ├── ids263.fasta
    ├── ids331.fasta
    ├── group_info.csv
    └── fromAliView.fasta
 
```

## Scripts

The scripts are intended to be run in numerical order.

| Script                                       | Purpose                                                      |
| -------------------------------------------- | ------------------------------------------------------------ |
| `01_DADA_pipeline.r`                         | DADA2 processing of raw FASTQ files for both rbcL markers.   |
| `02_technical_replicate_consistency_FigS1.R` | Technical replicate consistency analysis.                    |
| `03_rarefaction_procrustes_FigS2.R`          | Evaluation of the effect of rarefaction.                     |
| `04_ASV_correlation_Fig1.R`                  | ASV-level comparison of shared 263-bp core sequences.        |
| `05_prevalence_abundance_FigS3.R`            | Prevalence and abundance of shared and marker-specific ASVs. |
| `06_primer_mismatch_FigS4.R`                 | Primer mismatch analysis.                                    |
| `07_species_taxonomy_community_FigS5.R`      | Species-level taxonomic assignment and community comparison. |
| `08_IPS_analysis_Fig2_Fig3.R`                | IPS calculation and decomposition of IPS differences.        |

## Requirements

The analyses were run in R and require packages used in the scripts, including mainly:

* `dada2`
* `ShortRead`
* `Biostrings`
* `phyloseq`
* `vegan`
* `tidyverse`
* `readxl`
* `ggplot2`

Primer trimming in the DADA2 pipeline requires **Cutadapt**: Martin M. (2011) *Cutadapt removes adapter sequences from high-throughput sequencing reads*. EMBnet.journal, 17(1), 10–12. [https://doi.org/10.14806/ej.17.1.200](https://doi.org/10.14806/ej.17.1.200)

## Notes

The first script generates the DADA2 output files used by the downstream analyses. Scripts 02–08 reproduce the figures and analyses reported in the manuscript.

Two Diat.Barcode v15.2 reference files used by scripts 07 and 08 are not stored in this repository and should be downloaded separately to the input/ folder:

2025-09-04-Diat.barcode_release-version 15.2.xlsx — Diat.Barcode v15.2 release file: https://doi.org/10.15454/TOMBYZ

diat_barcode_v15_2_tax_assign_dada2.fa — DADA2-formatted taxonomic reference file derived from Diat.Barcode: https://doi.org/10.57745/SE6GJH

Methodological details and interpretation of results are provided in the manuscript.

## Citation

If you use this repository, please cite the associated article:

Urbánková P. & Štěpka J. (in press). *Evaluation of two overlapping rbcL amplicons for routine diatom DNA biomonitoring*. *Metabarcoding and Metagenomics*.

This citation will be updated with the DOI after publication.

