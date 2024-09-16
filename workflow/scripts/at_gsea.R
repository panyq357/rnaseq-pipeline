library(clusterProfiler)


config <- list(

    # input
    de_res = snakemake@input$de_res,
    at_go = snakemake@input$at_go,

    # output
    out_xlsx = snakemake@output$out_xlsx,
    out_rds = snakemake@output$out_rds,
    out_pdf = snakemake@output$out_pdf
)


main <- function(){

    de_res <- read.csv(config$de_res, row.names=1)

    at_go <- readr::read_csv(config$at_go)

    # Rank all genes by log2FoldChange for GSEA.
    gene_rank <- order(de_res$log2FoldChange, decreasing=T)
    ranked_gene_list <- setNames(de_res$log2FoldChange[gene_rank], row.names(de_res)[gene_rank])

    at_gsea_res_list = at_gsea(ranked_gene_list, at_go)

    writexl::write_xlsx(lapply(at_gsea_res_list, as.data.frame), config$out_xlsx)
    saveRDS(at_gsea_res_list, config$out_rds)
    
    pdf(config$out_pdf, width=10, height=10)
        for (name in names(at_gsea_res_list)) {
            print(
                enrichplot::ridgeplot(at_gsea_res_list[[name]], showCategory=20) +
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


at_gsea <- function(ranked_gene_list, at_go, category_list=c("BP", "MF", "CC")) {
    res_list <- list()
    for (cat in category_list) {
        onto_df <- subset(at_go, Category == cat)
        res <- df_gsea(ranked_gene_list, onto_df, "GeneID", "OntoID", "Description")
        res_list[[cat]] <- res
    }
    return(res_list)
}


main()

