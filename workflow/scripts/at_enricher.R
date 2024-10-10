library(clusterProfiler)


config <- list(
    de_res = snakemake@input$de_res,
    at_go = snakemake@input$at_go,
    p_column = snakemake@wildcards$p_column,
    p_th = as.numeric(snakemake@wildcards$p_th),
    abs_l2fc_th = log2(as.numeric(snakemake@wildcards$fc_th)),
    out_xlsx = snakemake@output$out_xlsx,
    out_rds = snakemake@output$out_rds,
    out_pdf = snakemake@output$out_pdf
)


main <- function(){

    de_res <- read.csv(config$de_res, row.names=1)

    at_go <- readr::read_csv(config$at_go)

    universe <- row.names(de_res)
    gene <- list(
        All  = de_res |> subset(de_res[[config$p_column]] < config$p_th & abs(log2FoldChange) >  config$abs_l2fc_th) |> row.names(),
        Up   = de_res |> subset(de_res[[config$p_column]] < config$p_th & log2FoldChange      >  config$abs_l2fc_th) |> row.names(),
        Down = de_res |> subset(de_res[[config$p_column]] < config$p_th & log2FoldChange      < -config$abs_l2fc_th) |> row.names()
    )

    at_enricher_res_list = lapply(gene, at_enricher, universe, at_go) |> unlist()

    out_df_list <- c(
        lapply(at_enricher_res_list, as.data.frame),
        setNames(lapply(gene, function(x) data.frame(Gene=x)), sprintf("GeneUsed.%s", names(gene)))
    )
    writexl::write_xlsx(out_df_list, config$out_xlsx)
    saveRDS(at_enricher_res_list, config$out_rds)

    pdf(config$out_pdf)
        for (name in names(at_enricher_res_list)) {
            print(
                enrichplot::dotplot(at_enricher_res_list[[name]]) +
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


at_enricher <- function(gene, universe, at_go, category_list=c("BP", "MF", "CC")) {
    res_list <- list()
    for (cat in category_list) {
        onto_df <- subset(at_go, Category == cat)
        res <- df_enricher(gene, universe, onto_df, "GeneID", "OntoID", "Description")
        res_list[[cat]] <- res
    }
    return(res_list)
}


main()

