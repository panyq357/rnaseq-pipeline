configfile: "config/os_anno-config.yaml"


rule add_annotation:
    input:
        de_res = "results/deseq2_pipeline/{name}.deseq2_results.csv",
        anno = config["rap-anno"]
    output:
        de_res_with_anno = "results/os_anno/{name}/{name}.deseq2_results_anno.csv"
    script:
        "../scripts/os_add_annotation.R"


rule oryzabase_gsea:
    input:
        de_res = "results/deseq2_pipeline/{name}.deseq2_results.csv",
        oryzabase_xlsx = config["oryzabase-xlsx"],
    output:
        out_xlsx = "results/os_anno/{name}/{name}.oryzabase-gsea.xlsx",
        out_rds = "results/os_anno/{name}/{name}.oryzabase-gsea.rds",
        out_pdf = "results/os_anno/{name}/{name}.oryzabase-gsea.pdf"
    params:
        oryzabase_sheet = config["oryzabase-sheet"]
    script:
        "../scripts/oryzabase_gsea.R"


rule oryzabase_enricher:
    input:
        de_res = "results/deseq2_pipeline/{name}.deseq2_results.csv",
        oryzabase_xlsx = config["oryzabase-xlsx"],
    output:
        out_xlsx = "results/os_anno/{name}/{name}.oryzabase-enricher.padj-{padj_th}.fc-{fc_th}.xlsx",
        out_rds = "results/os_anno/{name}/{name}.oryzabase-enricher.padj-{padj_th}.fc-{fc_th}.rds",
        out_pdf = "results/os_anno/{name}/{name}.oryzabase-enricher.padj-{padj_th}.fc-{fc_th}.pdf"
    params:
        oryzabase_sheet = config["oryzabase-sheet"]
    script:
        "../scripts/oryzabase_enricher.R"


rule os_kegg_gse:
    input:
        de_res = "results/deseq2_pipeline/{name}.deseq2_results.csv",
        dosa_id_converter = config["dosa-id-converter"]
    output:
        out_csv = "results/os_anno/{name}/{name}.kegg-gse.csv",
        out_rds = "results/os_anno/{name}/{name}.kegg-gse.rds",
        out_pdf = "results/os_anno/{name}/{name}.kegg-gse.pdf"
    script:
        "../scripts/os_kegg_gse.R"


rule os_kegg_enrich:
    input:
        de_res = "results/deseq2_pipeline/{name}.deseq2_results.csv",
        dosa_id_converter = config["dosa-id-converter"]
    output:
        out_xlsx = "results/os_anno/{name}/{name}.kegg-enrich.padj-{padj_th}.fc-{fc_th}.xlsx",
        out_rds = "results/os_anno/{name}/{name}.kegg-enrich.padj-{padj_th}.fc-{fc_th}.rds",
        out_pdf = "results/os_anno/{name}/{name}.kegg-enrich.padj-{padj_th}.fc-{fc_th}.pdf"
    script:
        "../scripts/os_kegg_enrich.R"

