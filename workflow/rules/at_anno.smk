configfile: "config/at_anno-config.yaml"


rule at_add_annotation:
    input:
        de_res = "results/deseq2_pipeline/{name}.deseq2_results.csv",
        tair_10_anno = config["tair-10-anno"],
    output:
        de_res_with_anno = "results/at_anno/{name}/{name}.deseq2_results_anno.csv"
    script:
        "../scripts/at_add_annotation.R"


rule at_clusterProfiler_anno:
    input:
        de_res = "results/deseq2_pipeline/{name}.deseq2_results.csv",
    output:
        xlsx = "results/at_anno/{name}/{name}.clusterProfiler.padj-{padj_th}.fc-{fc_th}.xlsx",
        rds = "results/at_anno/{name}/{name}.clusterProfiler.padj-{padj_th}.fc-{fc_th}.rds",
    script:
        "../scripts/at_clusterProfiler_anno.R"

