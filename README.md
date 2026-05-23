# Psilo_Splicing_GomezLab

This repository contains custom R scripts developed during a Master 2 internship in the Gomez Lab at UC Berkeley. The scripts were used to support candidate selection, annotation, and exploratory prioritization of psilocybin-associated skipped-exon events in mouse medial prefrontal cortex (mPFC) cell-type-specific RNA-seq data.

These scripts are intended to document and reproduce exploratory analysis steps performed during the internship. They are not designed as fully validated general-purpose bioinformatics pipelines.

## Scripts

### 1. `filter_and_classify_splicing_candidates.R`

This script filters skipped-exon alternative splicing events from a pre-generated rMATS output table.

Main steps:

- imports skipped-exon alternative splicing data;
- filters for PV interneuron events in the psilocybin condition;
- retains events with an absolute inclusion-level difference greater than or equal to 0.1;
- removes events with non-informative IncLevel values equal to 0, 1, missing, or non-numeric values;
- retrieves gene annotations using Ensembl BioMart;
- collapses annotations per gene;
- classifies candidate genes into exploratory functional categories using keyword-based scoring.

Output files:

- `results/gene_classification/gene_category_table.csv`
- `results/gene_classification/specific_filtering_candidates.csv`

The input file `allcelltypes_SE.csv` is not included in this repository because it belongs to the Gomez Lab dataset. To run this script, place the rMATS skipped-exon output file locally in the `data/` folder and name it `allcelltypes_SE.csv`.

### 2. `annotate_splicing_event_transcripts_and_motifs.R`

This script annotates a candidate alternative splicing event across all Ensembl mouse transcripts for a given gene.

Main steps:

- converts RefSeq mm39 chromosome accessions to Ensembl-style chromosome names;
- retrieves all Ensembl transcripts associated with a candidate gene;
- determines whether the candidate event overlaps an exon, an intron, or no annotated transcript region;
- retrieves protein motif/domain annotations from BioMart-linked resources, including InterPro, PFAM, SMART, PROSITE, and PRINTS;
- exports transcript-level annotation tables for selected candidate events.

Example output:

- `results/Brd4_transcript_event_protein_motif_annotation.csv`

## Folder structure

Psilo_Splicing_GomezLab  
├── README.md  
├── .gitignore  
├── scripts/  
│   ├── filter_and_classify_splicing_candidates.R  
│   └── annotate_splicing_event_transcripts_and_motifs.R  
├── data/  
│   └── .gitkeep  
└── results/  
    └── .gitkeep  

## Requirements

The scripts were developed in RStudio version 2024.12.0+467.

Required R packages:

- dplyr
- stringr
- purrr
- biomaRt
- GenomicRanges
- IRanges
- S4Vectors

## Important notes

These scripts require internet access because they query Ensembl BioMart. If Ensembl BioMart or one of its mirrors is temporarily unavailable or slow, the scripts may fail or take longer to run.

For reproducible execution, restart R before sourcing the scripts.

The input data file is not included in this repository. To run the filtering script, place the rMATS skipped-exon output file in the `data/` folder and name it `allcelltypes_SE.csv`.

## Project context

These scripts were used as part of an internship project investigating psilocybin-associated alternative splicing events in mouse mPFC parvalbumin interneurons. The goal was to move from cell-type-specific rMATS alternative splicing predictions toward experimentally testable candidate isoform events.
