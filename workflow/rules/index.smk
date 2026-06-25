# ─────────────────────────────────────────
# Index Rule
# Builds the STAR genome index from the
# dm6 reference FASTA and GTF annotation.
#
# IMPORTANT: On most HPC clusters and cloud
# environments, a pre-built STAR index is
# already available in a shared reference
# directory. If the path specified in
# config["reference"]["star_index"] already
# exists, Snakemake will skip this rule
# automatically.
#
# Only runs if the index directory is missing
# at the configured path. In that case,
# ensure you have:
#   ~16GB RAM
#   ~20-30 minutes
# ─────────────────────────────────────────

rule star_index:
    """
    Build STAR genome index from the dm6
    reference FASTA and GTF annotation.

    Skipped automatically if the index
    already exists at the path specified
    in config["reference"]["star_index"].

    --genomeSAindexNbases 12 is required
    for small genomes like Drosophila dm6.
    The default value is tuned for human
    and will produce warnings and suboptimal
    alignments if used with dm6.
    """
    input:
        fasta = config["reference"]["genome_fasta"],
        gtf   = config["reference"]["gtf"]
    output:
        directory(config["reference"]["star_index"])
    conda:
        "../envs/align.yaml"
    log:
        "logs/star_index/star_index.log"
    threads: 8
    shell:
        """
        mkdir -p {output}

        STAR \
            --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeDir {output} \
            --genomeFastaFiles {input.fasta} \
            --sjdbGTFfile {input.gtf} \
            --genomeSAindexNbases 12 \
            2> {log}
        """