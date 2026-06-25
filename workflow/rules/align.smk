# ─────────────────────────────────────────
# Alignment Rule
# STAR: splice-aware alignment of trimmed
# reads to the reference genome
# ─────────────────────────────────────────

rule star_align:
    """
    Align trimmed paired-end reads to the
    reference genome using STAR. Produces a
    sorted BAM file ready for quantification.
    STAR is splice-aware, meaning it correctly
    handles reads that span exon-exon junctions.
    """
    input:
        r1 = config["paths"]["trimmed"] + "{sample}_R1_paired.fastq.gz",
        r2 = config["paths"]["trimmed"] + "{sample}_R2_paired.fastq.gz",
        index = config["reference"]["star_index"]
    output:
        bam = config["paths"]["aligned"] + "{sample}/Aligned.sortedByCoord.out.bam",
        log_final = config["paths"]["aligned"] + "{sample}/Log.final.out"
    conda:
        "../../envs/align.yaml"
    log:
        "logs/star/{sample}.log"
    params:
        prefix  = config["paths"]["aligned"] + "{sample}/",
        threads = config["params"]["star"]["threads"],
        outSAMtype = config["params"]["star"]["outSAMtype"]
    shell:
        """
        STAR \
            --runThreadN {params.threads} \
            --genomeDir {input.index} \
            --readFilesIn {input.r1} {input.r2} \
            --readFilesCommand zcat \
            --outSAMtype {params.outSAMtype} \
            --outFileNamePrefix {params.prefix} \
            --outSAMattributes NH HI AS NM \
            2> {log}

        samtools index {output.bam}
        """