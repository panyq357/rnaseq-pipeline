configfile: "config/sample_check-config.yaml"

rule collect_star_counts:
    input:
        expand("results/star_mapping/{sample}.ReadsPerGene.out.tab", sample=config["samples"])
    output:
        "results/collect_star_counts/counts.csv"
    script:
        "../scripts/collect_star_counts.py"


rule sample_check:
    input:
        counts = "results/collect_star_counts/counts.csv"
    output:
        sample_dist_heatmap = "results/sample_check/sample_dist_heatmap.svg"
    params:
        sample_groups = config["sample_groups"]
    script:
        "../scripts/sample_check.R"

