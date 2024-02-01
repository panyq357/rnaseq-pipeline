library(clusterProfiler)


config <- list(
    de_res = snakemake@input$de_res,
    oryzabase_xlsx = snakemake@input$oryzabase_xlsx,
    oryzabase_sheet = snakemake@params$oryzabase_sheet,
    out_xlsx = snakemake@output$out_xlsx,
    out_rds = snakemake@output$out_rds,
    out_pdf = snakemake@output$out_pdf
)


main <- function(){

    de_res <- read.csv(config$de_res, row.names=1)

    # Rank all genes by log2FoldChange for GSEA.
    gene_rank <- order(de_res$log2FoldChange, decreasing=T)
    ranked_gene_list <- setNames(de_res$log2FoldChange[gene_rank], row.names(de_res)[gene_rank])

    oryzabase_gsea_res_list = oryzabase_gsea(ranked_gene_list, config$oryzabase_xlsx, config$oryzabase_sheet)

    writexl::write_xlsx(lapply(oryzabase_gsea_res_list, as.data.frame), config$out_xlsx)
    saveRDS(oryzabase_gsea_res_list, config$out_rds)
    
    pdf(config$out_pdf, width=10, height=10)
        for (name in names(oryzabase_gsea_res_list)) {
            print(
                enrichplot::ridgeplot(oryzabase_gsea_res_list[[name]]) +
                ggplot2::labs(title=name, x="log2(Fold Change)")
            )
        }
    dev.off()

}


df_gsea <- function(ranked_gene_list, onto_df, gene_col, term_col, name_col) {
    res <- clusterProfiler::GSEA(
        geneList = ranked_gene_list,
        TERM2GENE = onto_df[c(term_col, gene_col)],
        TERM2NAME = onto_df[c(term_col, name_col)],
        pvalueCutoff = 1,
    )
    return(res)
}


oryzabase_gsea <- function(ranked_gene_list, oryzabase_xlsx, sheet_list) {
    res_list <- list()
    for (sheet in sheet_list) {
        onto_df <- readxl::read_excel(oryzabase_xlsx, sheet=sheet)
        res <- df_gsea(ranked_gene_list, onto_df, "GeneID", "OntoID", "Description")
        res_list[[sheet]] <- res
    }
    return(res_list)
}


main()

