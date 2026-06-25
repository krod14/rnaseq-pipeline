# ─────────────────────────────────────────
# Quantification Rule
# featureCounts: counts reads per gene
# across all samples simultaneously,
# producing a count matrix for DESeq2
# ─────────────────────────────────────────

rule featurecounts:
    """
    Count reads mapped to each gene across
    all samples using featureCounts. Runs on
    all BAM files simultaneously to produce
    a single count matrix. The GTF annotation
    defines gene boundaries for assignment.
    """
    input:
        bams = expand(
            config["paths"]["aligned"] + "{sample}/Aligned.sortedByCoord.out.bam",
            sample=SAMPLES
        ),
        gtf = config["reference"]["gtf"]
    output:
        counts    = config["paths"]["counts"] + "counts.txt",
        counts_summary = config["paths"]["counts"] + "counts.txt.summary"
    conda:
        "../../envs/quantify.yaml"
    log:
        "logs/featurecounts/featurecounts.log"
    params:
        threads     = config["params"]["featurecounts"]["threads"],
        strandedness = config["params"]["featurecounts"]["strandedness"]
    shell:
        """
        featureCounts \
            -T {params.threads} \
            -p \
            -s {params.strandedness} \
            -a {input.gtf} \
            -o {output.counts} \
            {input.bams} \
            2> {log}
        """