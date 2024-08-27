# Plot a heatmap for checking sample groups.

rule sample_check:
    input:
        counts = "results/collect_star_counts/counts.csv"
    output:
        sample_dist_heatmap = "results/sample_check/sample_dist_heatmap.svg"
    params:
        sample_groups = config["sample_groups"]
    script:
        "../scripts/sample_check.R"

