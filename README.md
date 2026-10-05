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
quarto render dada2_template.qmd -P amplicon:16SV4 -P group:Sample_type
```

Available amplicon targets: `16SV3V4`, `16SV4`, `16SV4V5`, `ITS1`, `ITS2`, `12STELEO2`, `ADCR2`, `RBCL`.

The optional `group` parameter names the metadata column used to group samples in the example analysis (PCoA and Shannon diversity). The column must contain categories shared by several samples; IDs and continuous measurements are rejected. If `group` is omitted, the example analysis is skipped.

## Test data

`testdata/` holds two small synthetic MiSeq i100 datasets (5 samples, ~200 read pairs each, binned Q-scores) for quick test renders:

| Folder | Amplicon | Reads | Tests |
|---|---|---|---|
| `testdata/16SV4/` | 16S V4, gut-like bacterial community | 2×250 | fixed truncation lengths, chloroplast/mitochondria removal, tree building |
| `testdata/ITS2/` | ITS2, soil/rhizosphere fungal community | 2×300 | `truncLen = c(0, 0)`, variable amplicon lengths (273–405 bp), read-through and overhang trimming at merging |

To run a test set, copy its files into a project the same way as real data:

```bash
cp testdata/ITS2/metadata.tsv <project>/
mkdir -p <project>/reads
cp testdata/ITS2/reads/*.fastq.gz <project>/reads/
cd <project> && quarto render dada2_template.qmd -P amplicon:ITS2 -P group:Sample_type
```

- **Removing test data:** run `rm -rf reads/` in the project. This also removes the `reads/trimmed/` subfolder that cutadapt creates. If the trimmed files were left, they would be processed together with the next dataset.
- **Checking results:** compare `results/asv_table.tsv` and `results/taxonomy.tsv` with the folder's `truth.tsv`, which lists the true per-sample counts and sequences. `example_report.html` in each folder is a local reference render and is not committed.

The data are regenerated with `testdata/simulate_reads.py --region 16SV4|ITS2` (needs the local SILVA 138.2 or UNITE fasta and taxonomy files).

The reads are fully synthetic and contain no data from real samples. The 16S amplicon sequences were extracted in silico from the SILVA SSU Ref NR99 138.2 database (Quast et al. 2013, *Nucleic Acids Res.* 41:D590–D596, [https://www.arb-silva.de](https://www.arb-silva.de)), which is licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). The ITS2 sequences were extracted from the UNITE general FASTA release, dynamic species hypotheses (Abarenkov et al. 2024, *Nucleic Acids Res.* 52:D791–D797, [https://unite.ut.ee](https://unite.ut.ee)), which is licensed under [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/). UNITE sequences end at the ITS2/LSU boundary, so a conserved 5′ LSU segment up to the ITS4 primer site was appended to complete the amplicons.

## Requirements

- R packages: dada2, mia, scater, vegan, Biostrings, tidyverse, ape, phangorn, kableExtra, patchwork, ggthemes, ggrepel, ggpubr, ggsci, hrbrthemes, bslib, sessioninfo
- cutadapt, mafft and FastTree in the `python-311` conda environment (`~/.miniconda3/envs/python-311`); tools elsewhere on `PATH` are used as a fallback for mafft and FastTree, which are needed only for targets with tree building
- Reference databases under `~/reference/` (SILVA 138.2 trainset and species file, UNITE general release)

The template version is shown in the report subtitle.
