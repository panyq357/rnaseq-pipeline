configfile: "config/os_anno-config.yaml"


rule add_annotation:
    input:
        de_res = "results/deseq2_pipeline/{name}.deseq2_results.csv",
        anno = config["rap-anno"]
    output:
        de_res_with_anno = "results/os_anno/{name}.deseq2_results_anno.csv"
    script:
        "../scripts/os_add_annotation.R"


rule oryzabase_anno:
    input:
        de_res = "results/deseq2_pipeline/{name}.deseq2_results.csv",
        oryzabase_xlsx = config["oryzabase-xlsx"],
    output:
        enricher_xlsx = "results/os_anno/{name}/{name}.oryzabase-enricher.padj-{padj_th}.fc-{fc_th}.xlsx",
        enricher_rds = "results/os_anno/{name}/{name}.oryzabase-enricher.padj-{padj_th}.fc-{fc_th}.rds",
        GSEA_xlsx = "results/os_anno/{name}/{name}.oryzabase-GSEA.padj-{padj_th}.fc-{fc_th}.xlsx",
        GSEA_rds = "results/os_anno/{name}/{name}.oryzabase-GSEA.padj-{padj_th}.fc-{fc_th}.rds"
    params:
        oryzabase_sheet = config["oryzabase-sheet"]
    script:
        "../scripts/oryzabase_anno.R"


rule kegg_anno:
    input:
        de_res = "results/deseq2_pipeline/{name}.deseq2_results.csv",
    output:
        enrich_csv = "results/os_anno/{name}/{name}.kegg-enrich.padj-{padj_th}.fc-{fc_th}.csv",
        enrich_rds = "results/os_anno/{name}/{name}.kegg-enrich.padj-{padj_th}.fc-{fc_th}.rds",
        gse_csv = "results/os_anno/{name}/{name}.kegg-gse.padj-{padj_th}.fc-{fc_th}.csv",
        gse_rds = "results/os_anno/{name}/{name}.kegg-gse.padj-{padj_th}.fc-{fc_th}.rds"
    script:
        "../scripts/os_kegg_anno.R"

