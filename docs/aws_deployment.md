# AWS Deployment Guide

This document describes how to deploy and run the RNA-seq pipeline
on AWS using EC2 for compute and S3 for storage. This is the
recommended approach for full dataset (4 samples, paired-end) execution, 
as the pipeline requires ~16GB RAM for STAR index building and alignment.

---

## Architecture

- **EC2** — compute instance for running the pipeline
- **S3** — persistent storage for results and reference files
- **IAM** — access management for S3 permissions

---

## Prerequisites

- An AWS account
- A key pair for SSH access
- Basic familiarity with the AWS console

---

## Step 1: Create an S3 Bucket

1. Go to **S3 → Create bucket**
2. Name your bucket (e.g. `rnaseq-pipeline-kr`) — must be globally unique
3. Select your region (e.g. `us-east-2`)
4. Keep all other defaults (block all public access)
5. Click **Create bucket**

Create two folders inside the bucket:
- `results/` — for pipeline outputs
- `references/` — for reference genome and STAR index

---

## Step 2: Create an IAM User

1. Go to **IAM → Users → Create user**
2. Name the user (e.g. `rnaseq-pipeline`)
3. Click **Attach policies directly**
4. Search for and attach `AmazonS3FullAccess`
5. Click **Create user**

Then create an access key:
1. Click on your new user → **Security credentials**
2. Click **Create access key**
3. Select **Command Line Interface (CLI)**
4. Download or copy both the **Access Key ID** and **Secret Access Key**

---

## Step 3: Launch an EC2 Instance

1. Go to **EC2 → Launch Instance**
2. Configure:
   - **Name:** `rnaseq-pipeline`
   - **AMI:** Ubuntu Server 24.04 LTS
   - **Instance type:** `r5.xlarge` (4 vCPUs, 32GB RAM)
   - **Key pair:** create or select existing; download `.pem` file
   - **Storage:** 50GB gp3
   - **Security group:** allow SSH (port 22) from your IP
3. Click **Launch Instance**

---

## Step 4: Connect to the Instance

Move your key file and set permissions:
```bash
mv ~/Downloads/your-key.pem ~/.ssh/
chmod 400 ~/.ssh/your-key.pem
```

SSH in (get the public IP from the EC2 console):
```bash
ssh -i ~/.ssh/your-key.pem ubuntu@<public-ipv4>
```

---

## Step 5: Set Up the Environment

**Update the system:**
```bash
sudo apt-get update && sudo apt-get upgrade -y
```

**Install Miniconda:**
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
```

**Create Snakemake environment:**
```bash
conda create -n snakemake -c conda-forge -c bioconda snakemake=8.30.0 python=3.11 -y
conda activate snakemake
```

**Configure AWS CLI:**
```bash
curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
sudo apt install unzip -y
unzip awscliv2.zip
sudo ./aws/install
aws configure
```

Enter your IAM credentials when prompted:
- AWS Access Key ID
- AWS Secret Access Key
- Default region (e.g. `us-east-2`)
- Default output format (press Enter for default)

---

## Step 6: Clone the Repository and Stage Data

```bash
git clone https://github.com/krod14/rnaseq-pipeline.git
cd rnaseq-pipeline
```

Stage your FASTQ files (4 samples, paired-end) in `data/raw/` following
this naming convention:
```
{GSM_accession}_R1.fastq.gz
{GSM_accession}_R2.fastq.gz
```

For the Pasilla demo dataset (4 samples: 2 treated, 2 untreated),
download from NCBI SRA:
```bash
conda install -c conda-forge -c bioconda sra-tools -y

# Example for GSM461177 (repeat for each sample)
fasterq-dump SRR031714 --split-files --outdir data/raw/ --threads 4
fasterq-dump SRR031715 --split-files --outdir data/raw/ --threads 4

# Merge runs and compress
cat data/raw/SRR031714_1.fastq data/raw/SRR031715_1.fastq | gzip > data/raw/GSM461177_R1.fastq.gz
cat data/raw/SRR031714_2.fastq data/raw/SRR031715_2.fastq | gzip > data/raw/GSM461177_R2.fastq.gz

# Clean up
rm data/raw/SRR031714_*.fastq data/raw/SRR031715_*.fastq
```

Repeat for all four samples using the accessions in the README.

Stage the dm6 reference genome:
```bash
wget https://ftp.ensembl.org/pub/release-109/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa.gz \
    -O resources/genome/dm6.fa.gz
gunzip resources/genome/dm6.fa.gz

wget https://ftp.ensembl.org/pub/release-109/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.32.109.gtf.gz \
    -O resources/genome/dm6.gtf.gz
gunzip resources/genome/dm6.gtf.gz
```

---

## Step 7: Configure the Pipeline

Update `config/config.yaml` with your EC2 paths:

```yaml
reference:
  genome_fasta: "/home/ubuntu/rnaseq-pipeline/resources/genome/dm6.fa"
  gtf: "/home/ubuntu/rnaseq-pipeline/resources/genome/dm6.gtf"
  star_index: "/home/ubuntu/rnaseq-pipeline/resources/star_index/"

paths:
  raw_data: "/home/ubuntu/rnaseq-pipeline/data/raw/"
  trimmed: "/home/ubuntu/rnaseq-pipeline/data/trimmed/"
  aligned: "/home/ubuntu/rnaseq-pipeline/data/aligned/"
  counts: "/home/ubuntu/rnaseq-pipeline/data/counts/"
  results: "/home/ubuntu/rnaseq-pipeline/results/"
  logs: "/home/ubuntu/rnaseq-pipeline/logs/"
```

---

## Step 8: Run the Pipeline

Always use tmux to protect against SSH disconnections:

```bash
tmux new -s pipeline
conda activate rnaseq-pipeline
snakemake --cores 4 --use-conda
```

Detach with `Ctrl+B` then `D`. Reattach anytime with:
```bash
tmux attach -t pipeline
```

The pipeline will execute the following steps automatically:
1. Build the STAR genome index (~20-30 min, ~16GB RAM)
2. Run FastQC on all 4 samples in parallel
3. Trim all 4 samples with Trimmomatic
4. Align all 4 samples with STAR
5. Count reads with featureCounts across all samples
6. Run DESeq2 differential expression analysis

---

## Step 9: Sync Results to S3

After the pipeline completes:

```bash
# Sync results
aws s3 sync results/ s3://your-bucket-name/results/

# Sync reference files for future runs
aws s3 sync resources/ s3://your-bucket-name/references/
```

---

## Step 10: Clean Up

When finished, terminate the EC2 instance from the AWS console to
stop incurring compute charges. Your results are safely stored in S3.

**Important:** Click **Terminate** only when you are done. If you plan
to resume work, click **Stop** instead since this preserves your data but
stops compute charges.

---

## Troubleshooting

**SSH connection refused:**
- Check that port 22 is open in your security group
- Verify you are using the correct public IP (changes on restart)
- Ensure key file permissions are set to 400

**conda: command not found:**
```bash
source ~/.bashrc
```

**Snakemake not found:**
```bash
conda activate snakemake
```

**sra-tools SSL error:**
```bash
conda install -c conda-forge -c bioconda sra-tools -y
```

**STAR index not building:**
- Ensure `resources/star_index/` directory does not already exist
- If it does, remove it: `rm -rf resources/star_index/`

---

## Estimated Costs

| Resource | Cost |
|----------|------|
| r5.xlarge EC2 (~2 hours) | ~$0.50 |
| 50GB EBS storage (~2 hours) | ~$0.01 |
| S3 storage (results + references) | ~$0.30/month |
| Data transfer | ~$0.01 |
| **Total for one run** | **~$0.52** |

*Costs based on us-east-2 on-demand pricing as of 2026.*