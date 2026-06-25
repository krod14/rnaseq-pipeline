# ─────────────────────────────────────────
# Trimming Rule
# Trimmomatic: removes adapters and low
# quality bases from raw paired-end reads
# ─────────────────────────────────────────

rule trimmomatic:
    """
    Trim adapters and low-quality bases from
    paired-end reads. Produces paired output
    (both reads survived) and unpaired output
    (one read was dropped). We only use paired
    output in downstream steps.
    """
    input:
        r1 = config["paths"]["raw_data"] + "{sample}_R1.fastq.gz",
        r2 = config["paths"]["raw_data"] + "{sample}_R2.fastq.gz"
    output:
        r1_paired   = config["paths"]["trimmed"] + "{sample}_R1_paired.fastq.gz",
        r2_paired   = config["paths"]["trimmed"] + "{sample}_R2_paired.fastq.gz",
        r1_unpaired = config["paths"]["trimmed"] + "{sample}_R1_unpaired.fastq.gz",
        r2_unpaired = config["paths"]["trimmed"] + "{sample}_R2_unpaired.fastq.gz"
    conda:
        "../../envs/trim.yaml"
    log:
        "logs/trimmomatic/{sample}.log"
    params:
        leading    = config["params"]["trimmomatic"]["leading"],
        trailing   = config["params"]["trimmomatic"]["trailing"],
        slidingwindow = config["params"]["trimmomatic"]["slidingwindow"],
        minlen     = config["params"]["trimmomatic"]["minlen"]
    shell:
        """
        trimmomatic PE \
            {input.r1} {input.r2} \
            {output.r1_paired} {output.r1_unpaired} \
            {output.r2_paired} {output.r2_unpaired} \
            LEADING:{params.leading} \
            TRAILING:{params.trailing} \
            SLIDINGWINDOW:{params.slidingwindow} \
            MINLEN:{params.minlen} \
            2> {log}
        """