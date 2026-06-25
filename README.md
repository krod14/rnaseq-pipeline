# RNA-seq Analysis Pipeline

A reproducible, end-to-end bulk RNA-seq pipeline built with Snakemake,
conda, and Docker. Processes raw FASTQ reads through QC, trimming,
alignment, and quantification, culminating in differential expression
analysis and an interactive R/Shiny dashboard for exploration of results.

Designed for execution on HPC clusters and cloud environments (AWS, GCP),
with modular, parameterized configuration that adapts to any compute
environment by updating a single configuration file.

---

## Pipeline Overview

```
Raw FASTQs → STAR Index → FastQC/MultiQC → Trimmomatic → STAR → featureCounts → DESeq2 → Shiny Dashboard
```

| Step | Tool | Purpose |
|------|------|---------|
| Index | STAR | Build genome index from dm6 reference FASTA and GTF |
| QC | FastQC + MultiQC | Per-sample and aggregate quality reports |
| Trimming | Trimmomatic | Adapter removal and low-quality base trimming |
| Alignment | STAR | Splice-aware alignment to reference genome |
| Quantification | featureCounts | Gene-level read counting across all samples |
| Differential Expression | DESeq2 | Statistical modeling of count data |
| Visualization | R/Shiny + plotly | Interactive exploration of DE results |

---

## Repository Structure

```
rnaseq-pipeline/
├── workflow/
│   ├── Snakefile                   # Master workflow definition
│   ├── rules/                      # Modular rule files per pipeline stage
│   │   ├── index.smk               # STAR genome index construction
│   │   ├── qc.smk
│   │   ├── trim.smk
│   │   ├── align.smk
│   │   ├── quantify.smk
│   │   └── diffexp.smk
│   └── scripts/
│       └── deseq2.R                # DESeq2 analysis script
├── config/
│   └── config.yaml                 # Single configuration file for all parameters
├── envs/                           # Per-step conda environment definitions
│   ├── qc.yaml
│   ├── trim.yaml
│   ├── align.yaml
│   ├── quantify.yaml
│   └── diffexp.yaml
├── dashboard/                      # R/Shiny interactive results dashboard
│   ├── app.R
│   └── global.R
├── docs/
│   ├── data_dictionary.md          # Input/output documentation
│   └── aws_deployment.md           # AWS EC2/S3 deployment guide
└── Dockerfile                      # Containerized pipeline environment
```

---

## Quickstart

For deployment on AWS EC2, see the [AWS Deployment Guide](docs/aws_deployment.md).

### Step 1: Clone the repository

```bash
git clone https://github.com/krod14/rnaseq-pipeline.git
cd rnaseq-pipeline
```

### Step 2: Install Snakemake

```bash
conda create -n snakemake -c conda-forge -c bioconda snakemake=8.30.0 python=3.11 -y
conda activate snakemake
```

### Step 3: Configure paths

Update `config/config.yaml` with your cluster or cloud paths before
running the pipeline. There are three sections to update:

**Reference genome paths** — point to your cluster's shared reference
directory, or download dm6 from Ensembl release 109:
- Genome FASTA: [Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa.gz](https://ftp.ensembl.org/pub/release-109/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa.gz)
- GTF annotation: [Drosophila_melanogaster.BDGP6.32.109.gtf.gz](https://ftp.ensembl.org/pub/release-109/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.32.109.gtf.gz)

**STAR index path** — if a pre-built dm6 STAR index exists on your
cluster, point `star_index` to it and the index rule will be skipped
automatically. Otherwise the pipeline will build it on first run
(requires ~16GB RAM, ~20-30 minutes).

**Project data paths** — point `raw_data` to the directory containing
your staged FASTQ files.

### Step 4: Stage your FASTQ files

FASTQ files should be staged in the directory specified by
`config["paths"]["raw_data"]` and follow this naming convention:

```
{GSM_accession}_R1.fastq.gz
{GSM_accession}_R2.fastq.gz
```

For the Pasilla demo dataset, the expected files are:

```
GSM461177_R1.fastq.gz   GSM461177_R2.fastq.gz
GSM461178_R1.fastq.gz   GSM461178_R2.fastq.gz
GSM461180_R1.fastq.gz   GSM461180_R2.fastq.gz
GSM461181_R1.fastq.gz   GSM461181_R2.fastq.gz
```

See the [Demo Dataset](#demo-dataset) section for SRA accession numbers.

### Step 5: Run the pipeline

```bash
snakemake --cores 4 --use-conda
```

On an HPC cluster with SLURM:

```bash
snakemake --cores 4 --use-conda --executor slurm
```

### Step 6: Launch the Shiny dashboard

```r
# In R or RStudio
shiny::runApp("dashboard/")
```

---

### Run with Docker

```bash
docker build -t rnaseq-pipeline .

docker run -v $(pwd)/data:/pipeline/data \
           -v $(pwd)/results:/pipeline/results \
           -v $(pwd)/resources:/pipeline/resources \
           rnaseq-pipeline
```

The `-v` flags mount your local data and results directories into the
container so inputs and outputs persist after the container exits.

---

## Configuration

All pipeline parameters are controlled through `config/config.yaml`.
No hardcoded values exist in the rule files. Adapting the pipeline to
a new environment or dataset requires only updating this file.

```yaml
# Reference paths — update for your environment
reference:
  genome_fasta: "/path/to/shared/references/dm6/dm6.fa"
  gtf: "/path/to/shared/references/dm6/dm6.gtf"
  star_index: "/path/to/shared/references/dm6/star_index/"

# Sample conditions
samples:
  GSM461177: "untreated"
  GSM461178: "untreated"
  GSM461180: "treated"
  GSM461181: "treated"

# DESeq2 thresholds
params:
  deseq2:
    fdr_threshold: 0.05
    lfc_threshold: 1
```

---

## Demo Dataset

This pipeline is demonstrated using the **Pasilla dataset** — a
well-characterized *Drosophila melanogaster* RNA-seq benchmark comparing
RNAi knockdown of the *pasilla* gene against controls. It is the
canonical dataset used in DESeq2 documentation and produces reliable,
interpretable results. Original data available at NCBI GEO under
accession [GSE18508](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18508).

| GSM | SRX | SRR Runs | Condition |
|-----|-----|----------|-----------|
| [GSM461177](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM461177) | [SRX014459](https://www.ncbi.nlm.nih.gov/sra/SRX014459) | SRR031714, SRR031715 | untreated |
| [GSM461178](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM461178) | [SRX014460](https://www.ncbi.nlm.nih.gov/sra/SRX014460) | SRR031716, SRR031717 | untreated |
| [GSM461180](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM461180) | [SRX014462](https://www.ncbi.nlm.nih.gov/sra/SRX014462) | SRR031724, SRR031725 | treated |
| [GSM461181](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM461181) | [SRX014463](https://www.ncbi.nlm.nih.gov/sra/SRX014463) | SRR031726, SRR031727 | treated |

Raw FASTQ files can be downloaded from NCBI SRA using `fasterq-dump`
from the [sra-tools](https://github.com/ncbi/sra-tools) package:

```bash
# Example for one sample — repeat for each SRR accession
fasterq-dump SRR031714 --split-files --outdir data/raw/ --threads 4
fasterq-dump SRR031715 --split-files --outdir data/raw/ --threads 4

# Merge runs and compress
cat data/raw/SRR031714_1.fastq data/raw/SRR031715_1.fastq | gzip > data/raw/GSM461177_R1.fastq.gz
cat data/raw/SRR031714_2.fastq data/raw/SRR031715_2.fastq | gzip > data/raw/GSM461177_R2.fastq.gz
```

The above example shows GSM461177 only. Repeat for GSM461178, GSM461180,
and GSM461181 using their respective SRR accessions from the table above.
For complete step-by-step download instructions for all four samples, see the
[AWS Deployment Guide](docs/aws_deployment.md).

---

## Shiny Dashboard

The interactive dashboard provides four views of the DESeq2 results:

- **Overview** — summary statistics and PCA plot for sample QC
- **Volcano Plot** — interactive, filterable by FDR and fold change cutoffs
- **Heatmap** — top N most variable genes with adjustable scaling
- **Results Table** — full DESeq2 output with column filtering and sorting

---

## Reproducibility

Each pipeline step runs in its own isolated conda environment defined
in `envs/`. The entire pipeline is containerized via Docker for
system-level reproducibility. All parameters are version-pinned.

To reproduce results exactly:
```bash
docker build -t rnaseq-pipeline .
docker run -v $(pwd)/data:/pipeline/data \
           -v $(pwd)/results:/pipeline/results \
           -v $(pwd)/resources:/pipeline/resources \
           rnaseq-pipeline
```

This pipeline was validated end-to-end on AWS EC2 (`r5.xlarge`, 4 vCPUs,
32GB RAM) using the four paired-end samples and dm6 reference genome. 
All 17 pipeline steps completed successfully, producing
FastQC reports, a MultiQC summary, DESeq2 differential expression results,
and visualization outputs. See [docs/aws_deployment.md](docs/aws_deployment.md)
for full deployment instructions.

---

## Requirements

- Snakemake ≥ 8.0
- conda or mamba
- Docker (optional, for containerized execution)
- 16GB RAM recommended for STAR index building and alignment
- HPC cluster or cloud instance recommended for full dataset execution

---

## Author

**Kyle Rodrigues**
M.S. Bioinformatics, Georgetown University
[LinkedIn](https://www.linkedin.com/in/rodrigueskyle) | kyle.r.rodrigues@gmail.com