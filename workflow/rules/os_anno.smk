
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
        oryzabase_xlsx = config["oryzabase"]["xlsx"],
    output:
        out_xlsx = "results/os_anno/{name}/{name}.oryzabase-gsea.xlsx",
        out_rds = "results/os_anno/{name}/{name}.oryzabase-gsea.rds",
        out_pdf = "results/os_anno/{name}/{name}.oryzabase-gsea.pdf"
    params:
        oryzabase_sheet = config["oryzabase"]["sheet"]
    script:
        "../scripts/oryzabase_gsea.R"


rule oryzabase_enricher:
    input:
        de_res = "results/deseq2_pipeline/{name}.deseq2_results.csv",
        oryzabase_xlsx = config["oryzabase"]["xlsx"],
    output:
        out_xlsx = "results/os_anno/{name}/{name}.oryzabase-enricher.padj-{padj_th}.fc-{fc_th}.xlsx",
        out_rds = "results/os_anno/{name}/{name}.oryzabase-enricher.padj-{padj_th}.fc-{fc_th}.rds",
        out_pdf = "results/os_anno/{name}/{name}.oryzabase-enricher.padj-{padj_th}.fc-{fc_th}.pdf"
    params:
        oryzabase_sheet = config["oryzabase"]["sheet"]
    script:
        "../scripts/oryzabase_enricher.R"


rule os_kegg_gsea:
    input:
        de_res = "results/deseq2_pipeline/{name}.deseq2_results.csv",
        pathway_to_gene = config["kegg"]["pathway_to_gene"],
        pathway_to_description = config["kegg"]["pathway_to_description"]
    output:
        out_xlsx = "results/os_anno/{name}/{name}.kegg-gsea.xlsx",
        out_rds = "results/os_anno/{name}/{name}.kegg-gsea.rds",
        out_pdf = "results/os_anno/{name}/{name}.kegg-gsea.pdf"
    script:
        "../scripts/os_kegg_gsea.R"


rule os_kegg_enricher:
    input:
        de_res = "results/deseq2_pipeline/{name}.deseq2_results.csv",
        pathway_to_gene = config["kegg"]["pathway_to_gene"],
        pathway_to_description = config["kegg"]["pathway_to_description"]
    output:
        out_xlsx = "results/os_anno/{name}/{name}.kegg-enricher.padj-{padj_th}.fc-{fc_th}.xlsx",
        out_rds = "results/os_anno/{name}/{name}.kegg-enricher.padj-{padj_th}.fc-{fc_th}.rds",
        out_pdf = "results/os_anno/{name}/{name}.kegg-enricher.padj-{padj_th}.fc-{fc_th}.pdf"
    script:
        "../scripts/os_kegg_enricher.R"

