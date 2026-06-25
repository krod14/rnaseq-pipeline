# ─────────────────────────────────────────
# Dockerfile for RNA-seq Pipeline
# Base: Mambaforge (conda-forge defaults)
# Includes: Snakemake + all pipeline envs
# ─────────────────────────────────────────

# Start from a minimal conda base image
# Mambaforge uses mamba (faster conda solver)
FROM condaforge/mambaforge:23.3.1-1

# Metadata
LABEL maintainer="kyle.r.rodrigues@gmail.com"
LABEL description="Reproducible RNA-seq pipeline using Snakemake"
LABEL version="1.0"

# Set working directory inside the container
WORKDIR /pipeline

# ── Install Snakemake ──────────────────────
# We install Snakemake first as it's the
# orchestration layer for everything else
RUN mamba install -n base -c conda-forge -c bioconda \
        snakemake=8.30.0 \
        -y && \
    mamba clean --all -y

# ── Copy pipeline files ────────────────────
# Copy the full repo into the container
COPY workflow/  ./workflow/
COPY config/    ./config/
COPY envs/      ./envs/
COPY dashboard/ ./dashboard/
COPY docs/      ./docs/

# ── Pre-build conda environments ──────────
# Build all per-step environments at image
# build time so runtime has no setup delay.
# Snakemake will find these cached envs
# automatically when the pipeline runs.
RUN snakemake --cores 1 \
        --use-conda \
        --conda-create-envs-only \
        --snakefile workflow/Snakefile

# ── Default command ────────────────────────
# Running the container executes the pipeline
# Override with: docker run <image> snakemake --help
CMD ["snakemake", \
     "--cores", "4", \
     "--use-conda", \
     "--snakefile", "workflow/Snakefile"]