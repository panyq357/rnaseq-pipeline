library(clusterProfiler)


config <- list(
    de_res = snakemake@input$de_res,
    oryzabase_xlsx = snakemake@input$oryzabase_xlsx,
    padj_th = as.numeric(snakemake@wildcards$padj_th),
    abs_l2fc_th = log2(as.numeric(snakemake@wildcards$fc_th)),
    oryzabase_sheet = snakemake@params$oryzabase_sheet,
    out_xlsx = snakemake@output$out_xlsx,
    out_rds = snakemake@output$out_rds,
    out_pdf = snakemake@output$out_pdf
)


main <- function(){

    de_res <- read.csv(config$de_res, row.names=1)

    universe <- row.names(de_res)
    gene <- list(
        All  = de_res |> subset(padj < config$padj_th & abs(log2FoldChange) >  config$abs_l2fc_th) |> row.names(),
        Up   = de_res |> subset(padj < config$padj_th & log2FoldChange      >  config$abs_l2fc_th) |> row.names(),
        Down = de_res |> subset(padj < config$padj_th & log2FoldChange      < -config$abs_l2fc_th) |> row.names()
    )

    oryzabase_enricher_res_list = lapply(gene, oryzabase_enricher, universe, config$oryzabase_xlsx, config$oryzabase_sheet) |> unlist()

    out_df_list <- c(
        lapply(oryzabase_enricher_res_list, as.data.frame),
        setNames(lapply(gene, function(x) data.frame(Gene=x)), sprintf("GeneUsed.%s", names(gene)))
    )
    writexl::write_xlsx(out_df_list, config$out_xlsx)
    saveRDS(oryzabase_enricher_res_list, config$out_rds)

    pdf(config$out_pdf)
        for (name in names(oryzabase_enricher_res_list)) {
            print(
                enrichplot::dotplot(oryzabase_enricher_res_list[[name]]) +
                ggplot2::labs(title=name)
            )
        }
    dev.off()

}


df_enricher <- function(gene, universe, onto_df, gene_col, term_col, name_col) {
    res <- clusterProfiler::enricher(
        gene = gene,
        universe = universe,
        TERM2GENE = onto_df[c(term_col, gene_col)],
        TERM2NAME = onto_df[c(term_col, name_col)],
        pvalueCutoff = 1,
        qvalueCutoff = 1
    )
    return(res)
}


oryzabase_enricher <- function(gene, universe, oryzabase_xlsx, sheet_list) {
    res_list <- list()
    for (sheet in sheet_list) {
        onto_df <- readxl::read_excel(oryzabase_xlsx, sheet=sheet)
        res <- df_enricher(gene, universe, onto_df, "GeneID", "OntoID", "Description")
        res_list[[sheet]] <- res
    }
    return(res_list)
}


main()

