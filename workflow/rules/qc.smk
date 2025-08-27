rule cpm_cor_heatmap:
    input:
        "results/counts.csv"
    output:
        "results/cpm_cor_heatmap.svg"
    script:
        "../scripts/cpm_cor_heatmap.R"
