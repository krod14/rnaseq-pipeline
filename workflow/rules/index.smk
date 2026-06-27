# ─────────────────────────────────────────
# Index Rule
# Builds the STAR genome index from the
# dm6 reference FASTA and GTF annotation.
#
# IMPORTANT: On most HPC clusters and cloud
# environments, a pre-built STAR index is
# already available in a shared reference
# directory. If so, point star_index to it
# in config.yaml and this rule will be
# skipped automatically.
#
# The pipeline checks for genomeParameters.txt
# inside the star_index directory to determine
# if the index needs to be built. An empty
# directory will still trigger a build.
#
# If building from scratch, ensure you have:
#   ~16GB RAM
#   ~20-30 minutes
# ─────────────────────────────────────────

rule star_index:
    """
    Build STAR genome index from the dm6
    reference FASTA and GTF annotation.

    Skipped automatically if
    genomeParameters.txt already exists
    inside the configured star_index path.
    An empty directory will still trigger
    a build.

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
        dir    = directory(config["reference"]["star_index"]),
        params = config["reference"]["star_index"] + "genomeParameters.txt"
    conda:
        "../../envs/align.yaml"
    log:
        "logs/star_index/star_index.log"
    threads: 8
    shell:
        """
        mkdir -p {output.dir}

        STAR \
            --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeDir {output.dir} \
            --genomeFastaFiles {input.fasta} \
            --sjdbGTFfile {input.gtf} \
            --genomeSAindexNbases 12 \
            2> {log}
        """