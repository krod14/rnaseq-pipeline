# Data Dictionary

## Pipeline Overview
This document describes all input and output files produced by the
RNA-seq pipeline, including file formats, naming conventions, and
column definitions where applicable.

All pipeline parameters are defined in `config/config.yaml`. See inline
comments in that file for parameter descriptions.

---

## Setup

Before running the pipeline, update `config/config.yaml` with your
HPC cluster or cloud environment paths. For full deployment instructions
on AWS EC2, see [docs/aws_deployment.md](aws_deployment.md).

There are three sections to configure:

**1. Reference genome paths**
On most HPC clusters, reference genomes are maintained centrally in a
shared directory accessible to all users (e.g. `/scratch/shared/references/`).
Point the `reference` section of `config.yaml` to your cluster's dm6
FASTA, GTF, and STAR index paths.

If these files are not available on your cluster, download them from
Ensembl release 109:
- [dm6 FASTA](https://ftp.ensembl.org/pub/release-109/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa.gz)
- [dm6 GTF](https://ftp.ensembl.org/pub/release-109/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.32.109.gtf.gz)

**2. STAR index path**
If a pre-built dm6 STAR index exists on your cluster, point
`star_index` to it and the index rule will be skipped automatically.
Otherwise the pipeline builds it on first run (~16GB RAM, ~20-30 min).

**3. Project data paths**
Point `raw_data` to the directory containing your staged FASTQ files.
All other paths (`trimmed`, `aligned`, `counts`, `results`, `logs`)
should point to writable project directories on the cluster.

---

## Inputs

### Raw FASTQ Files
**Location:** `config["paths"]["raw_data"]`
**Format:** `.fastq.gz` (gzip-compressed FASTQ)
**Naming:** `{GSM_accession}_R1.fastq.gz`, `{GSM_accession}_R2.fastq.gz`
**Source:** NCBI SRA — see demo dataset accessions in README
**Note:** Files are expected to be staged in the raw data directory
before running the pipeline

| File | Description |
|------|-------------|
| `{sample}_R1.fastq.gz` | Forward reads (read 1 of paired-end) |
| `{sample}_R2.fastq.gz` | Reverse reads (read 2 of paired-end) |

### Reference Genome
**Location:** `config["reference"]["genome_fasta"]`
**Genome:** Drosophila melanogaster dm6 (BDGP Release 6)
**Source:** Ensembl release 109

| File | Description |
|------|-------------|
| `dm6.fa` | Reference genome FASTA |
| `dm6.gtf` | Gene annotation file (defines gene/exon boundaries) |

### STAR Index
**Location:** `config["reference"]["star_index"]`
**Built by:** `workflow/rules/index.smk` — only if index does not
already exist at the configured path
**Note:** On HPC clusters a pre-built index is often available in a
shared reference directory. If so, the index rule is skipped automatically.
Building from scratch requires ~16GB RAM and approximately 20-30 minutes.

---

## Intermediate Files

### Trimmed FASTQ
**Location:** `config["paths"]["trimmed"]`
**Tool:** Trimmomatic

| File | Description |
|------|-------------|
| `{sample}_R1_paired.fastq.gz` | Trimmed forward reads where both pairs survived |
| `{sample}_R2_paired.fastq.gz` | Trimmed reverse reads where both pairs survived |
| `{sample}_R1_unpaired.fastq.gz` | Forward reads whose reverse pair was dropped |
| `{sample}_R2_unpaired.fastq.gz` | Reverse reads whose forward pair was dropped |

### Aligned BAM Files
**Location:** `config["paths"]["aligned"]/{sample}/`
**Tool:** STAR
**Format:** BAM (Binary Alignment Map), sorted by coordinate

| File | Description |
|------|-------------|
| `Aligned.sortedByCoord.out.bam` | Coordinate-sorted alignment file |
| `Aligned.sortedByCoord.out.bam.bai` | BAM index (produced by samtools) |
| `Log.final.out` | STAR alignment summary statistics |

### Count Matrix
**Location:** `config["paths"]["counts"]`
**Tool:** featureCounts

| File | Description |
|------|-------------|
| `counts.txt` | Raw gene-level read counts; genes × samples matrix |
| `counts.txt.summary` | Per-sample read assignment statistics |

---

## Outputs

### QC Reports
**Location:** `config["paths"]["results"]/qc/`

| File | Description |
|------|-------------|
| `{sample}_R1_fastqc.html` | FastQC report for forward reads |
| `{sample}_R2_fastqc.html` | FastQC report for reverse reads |
| `multiqc_report.html` | Aggregated QC report across all samples and pipeline stages |

### Differential Expression
**Location:** `config["paths"]["results"]/diffexp/`

| File | Description |
|------|-------------|
| `deseq2_results.csv` | Full DESeq2 results table (see columns below) |
| `normalized_counts.csv` | VST-normalized counts for visualization |
| `pca_plot.pdf` | PCA of VST-normalized counts colored by condition |
| `volcano_plot.pdf` | Volcano plot of log2FC vs -log10(padj) |
| `heatmap.pdf` | Heatmap of top 50 most variable genes |

#### `deseq2_results.csv` Column Definitions

| Column | Type | Description |
|--------|------|-------------|
| `gene_id` | string | Gene identifier from GTF annotation |
| `baseMean` | float | Mean normalized count across all samples |
| `log2FoldChange` | float | Log2 fold change (treated / untreated) |
| `lfcSE` | float | Standard error of the log2 fold change |
| `stat` | float | Wald test statistic |
| `pvalue` | float | Raw p-value |
| `padj` | float | Benjamini-Hochberg adjusted p-value (FDR) |

#### `normalized_counts.csv` Column Definitions

| Column | Type | Description |
|--------|------|-------------|
| `gene_id` | string | Ensembl/FlyBase gene identifier (e.g. FBgn0000490) |
| `GSM461177` | float | VST-normalized count for untreated replicate 1 |
| `GSM461178` | float | VST-normalized count for untreated replicate 2 |
| `GSM461180` | float | VST-normalized count for treated replicate 1 |
| `GSM461181` | float | VST-normalized count for treated replicate 2 |

---

### Shiny Dashboard
**Location:** `dashboard/`
**Inputs:** `deseq2_results.csv` and `normalized_counts.csv`
**Launch:** `shiny::runApp("dashboard/")` in R or RStudio

| Tab | Description |
|-----|-------------|
| Overview | Summary statistics and PCA plot |
| Volcano Plot | Interactive, filterable by FDR and fold change |
| Heatmap | Top N most variable genes with adjustable scaling |
| Results Table | Full DESeq2 output with column filtering and sorting |

---

## Logs
**Location:** `config["paths"]["logs"]/{rule}/{sample}.log`

Each pipeline rule writes stderr to a dedicated log file. Logs are
useful for debugging failed jobs and auditing pipeline runs.

| Log | Rule | Description |
|-----|------|-------------|
| `logs/star_index/star_index.log` | `index.smk` | STAR index build output |
| `logs/fastqc/{sample}.log` | `qc.smk` | FastQC per-sample output |
| `logs/multiqc/multiqc.log` | `qc.smk` | MultiQC aggregation output |
| `logs/trimmomatic/{sample}.log` | `trim.smk` | Trimmomatic trimming stats |
| `logs/star/{sample}.log` | `align.smk` | STAR alignment output |
| `logs/featurecounts/featurecounts.log` | `quantify.smk` | featureCounts output |
| `logs/deseq2/deseq2.log` | `diffexp.smk` | DESeq2 R session output |