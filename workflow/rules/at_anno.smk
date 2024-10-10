rule at_add_annotation:
    input:
        de_res = "results/deseq2_pipeline/{name}.deseq2_results.csv",
        anno = config["at_anno"]
    output:
        de_res_with_anno = "results/at_anno/{name}/{name}.deseq2_results_anno.csv"
    script:
        "../scripts/at_add_annotation.R"


rule at_gsea:
    input:
        de_res = "results/deseq2_pipeline/{name}.deseq2_results.csv",
        at_go = config["at_go"],
    output:
        out_xlsx = "results/at_anno/{name}/{name}.at-gsea.xlsx",
        out_rds = "results/at_anno/{name}/{name}.at-gsea.rds",
        out_pdf = "results/at_anno/{name}/{name}.at-gsea.pdf"
    script:
        "../scripts/at_gsea.R"


rule at_enricher:
    input:
        de_res = "results/deseq2_pipeline/{name}.deseq2_results.csv",
        at_go = config["at_go"],
    output:
        out_xlsx = "results/at_anno/{name}/{name}.at-enricher.{p_column}-{p_th}.fc-{fc_th}.xlsx",
        out_rds = "results/at_anno/{name}/{name}.at-enricher.{p_column}-{p_th}.fc-{fc_th}.rds",
        out_pdf = "results/at_anno/{name}/{name}.at-enricher.{p_column}-{p_th}.fc-{fc_th}.pdf"
    script:
        "../scripts/at_enricher.R"


rule at_kegg_gsea:
    input:
        de_res = "results/deseq2_pipeline/{name}.deseq2_results.csv",
        pathway_to_gene = config["at_kegg"]["pathway_to_gene"],
        pathway_to_description = config["at_kegg"]["pathway_to_description"]
    output:
        out_xlsx = "results/at_anno/{name}/{name}.kegg-gsea.xlsx",
        out_rds = "results/at_anno/{name}/{name}.kegg-gsea.rds",
        out_pdf = "results/at_anno/{name}/{name}.kegg-gsea.pdf"
    script:
        "../scripts/kegg_gsea.R"


rule at_kegg_enricher:
    input:
        de_res = "results/deseq2_pipeline/{name}.deseq2_results.csv",
        pathway_to_gene = config["at_kegg"]["pathway_to_gene"],
        pathway_to_description = config["at_kegg"]["pathway_to_description"]
    output:
        out_xlsx = "results/at_anno/{name}/{name}.kegg-enricher.{p_column}-{p_th}.fc-{fc_th}.xlsx",
        out_rds = "results/at_anno/{name}/{name}.kegg-enricher.{p_column}-{p_th}.fc-{fc_th}.rds",
        out_pdf = "results/at_anno/{name}/{name}.kegg-enricher.{p_column}-{p_th}.fc-{fc_th}.pdf"
    script:
        "../scripts/kegg_enricher.R"

