# dada2 amplicon template

Quarto template for processing Illumina paired-end amplicon data with DADA2: primer trimming (cutadapt), filtering, denoising, merging, chimera removal, taxonomy assignment, an optional phylogenetic tree, and export into a `TreeSummarizedExperiment` object.

## Start a new project

Create a project folder (e.g. an RStudio project), then run in its root:

```bash
curl -fsSL https://raw.githubusercontent.com/msuokas-bco/dada2_template/main/new_project.sh | bash
```

This downloads the latest `dada2_template.qmd` and creates the `reads/` and `results/` directories. An existing template is not overwritten unless `FORCE=1` is given:

```bash
curl -fsSL https://raw.githubusercontent.com/msuokas-bco/dada2_template/main/new_project.sh | FORCE=1 bash
```

## Project layout

```
project/
├── dada2_template.qmd
├── metadata.tsv      # tab-separated, must contain a 'sampleid' column (any letter case)
├── reads/            # demultiplexed *_R1_*.fastq.gz / *_R2_*.fastq.gz
└── results/          # generated: tables, rds objects, logs
```

Sample IDs are taken from the first underscore-delimited token of the FASTQ file names and must match `sampleid` in the metadata.

## Render

```bash
quarto render dada2_template.qmd -P amplicon:16SV4
```

Available amplicon targets: `16SV3V4`, `16SV4`, `16SV4V5`, `ITS1`, `ITS2`, `12STELEO2`, `ADCR2`, `RBCL`.

The example analysis section uses `group_var`, which must be set to a column in your metadata.

## Requirements

- R packages: dada2, mia, scater, vegan, Biostrings, tidyverse, ape, phangorn, kableExtra, patchwork, ggthemes, ggrepel, ggpubr, ggsci, hrbrthemes, bslib, sessioninfo
- cutadapt (activated through conda in the trimming chunk)
- mafft and FastTree in `PATH` (only for targets with tree building)
- Reference databases under `~/reference/` (SILVA 138.2 trainset and species file, UNITE general release)

The template version is shown in the report subtitle.
