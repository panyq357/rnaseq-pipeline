
# Collect all samples' counts to a single counts matrix.

rule collect_star_counts:
    input:
        expand("results/star_mapping/{sample}.ReadsPerGene.out.tab", sample=config["samples"])
    output:
        "results/collect_star_counts/counts.csv"
    script:
        "../scripts/collect_star_counts.py"


rule deseq2:
    input:
        counts = "results/collect_star_counts/counts.csv"
    output:
        out_rds = "results/deseq2_pipeline/{name}.deseq2_results.rds",
        out_csv = "results/deseq2_pipeline/{name}.deseq2_results.csv"
    params:
        levels = lambda w: config["levels"][w.name],
        sample_groups = config["sample_groups"]
    script:
        "../scripts/deseq2_contrast.R"


rule volcano_plot:
    input:
        de_res = "results/deseq2_pipeline/{name}.deseq2_results.csv"
    output:
        plot = "results/volcano_plot/{name}.padj-{padj_th}.fc-{fc_th}.volcano.pdf"
    script:
        "../scripts/volcano_plot.R"

