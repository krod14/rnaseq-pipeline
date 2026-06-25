# ─────────────────────────────────────────
# QC Rules
# FastQC: per-sample quality report
# MultiQC: aggregates all reports into one
# ─────────────────────────────────────────

rule fastqc:
    """
    Run FastQC on raw FASTQ files for each sample.
    Produces an HTML report and a zip of metrics.
    """
    input:
        r1 = config["paths"]["raw_data"] + "{sample}_R1.fastq.gz",
        r2 = config["paths"]["raw_data"] + "{sample}_R2.fastq.gz"
    output:
        html_r1 = "results/qc/{sample}_R1_fastqc.html",
        html_r2 = "results/qc/{sample}_R2_fastqc.html",
        zip_r1  = "results/qc/{sample}_R1_fastqc.zip",
        zip_r2  = "results/qc/{sample}_R2_fastqc.zip"
    conda:
        "../envs/qc.yaml"
    log:
        "logs/fastqc/{sample}.log"
    shell:
        """
        fastqc {input.r1} {input.r2} \
            --outdir results/qc/ \
            --threads 2 \
            2> {log}
        """

rule multiqc:
    """
    Aggregate all FastQC reports into a single
    interactive HTML summary across all samples.
    """
    input:
        expand("results/qc/{sample}_R1_fastqc.zip", sample=SAMPLES),
        expand("results/qc/{sample}_R2_fastqc.zip", sample=SAMPLES)
    output:
        "results/qc/multiqc_report.html"
    conda:
        "../envs/qc.yaml"
    log:
        "logs/multiqc/multiqc.log"
    shell:
        """
        multiqc results/qc/ \
            --outdir results/qc/ \
            --filename multiqc_report.html \
            2> {log}
        """