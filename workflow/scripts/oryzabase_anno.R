library(clusterProfiler)


config <- list(
    de_res = snakemake@input$de_res,
    oryzabase_xlsx = snakemake@input$oryzabase_xlsx,
    padj_th = as.numeric(snakemake@wildcards$padj_th),
    abs_l2fc_th = log2(as.numeric(snakemake@wildcards$fc_th)),
    oryzabase_sheet = snakemake@params$oryzabase_sheet,
    output = list(
        enricher_xlsx = snakemake@output$enricher_xlsx,
        enricher_rds = snakemake@output$enricher_rds,
        GSEA_xlsx = snakemake@output$GSEA_xlsx,
        GSEA_rds = snakemake@output$GSEA_rds
    )
)

main <- function(){

    de_res <- read.csv(config$de_res, row.names=1)

    universe <- row.names(de_res)
    gene <- de_res |> subset(padj < config$padj_th & abs(log2FoldChange) > config$abs_l2fc_th) |> row.names()

    oryzabase_enricher_res_list = oryzabase_enricher(gene, universe, config$oryzabase_xlsx, config$oryzabase_sheet)

    writexl::write_xlsx(lapply(oryzabase_enricher_res_list, as.data.frame), config$output$enricher_xlsx)
    saveRDS(oryzabase_enricher_res_list, config$output$enricher_rds)

    # Rank all genes by log2FoldChange for GSEA.
    gene_rank <- order(de_res$log2FoldChange, decreasing=T)
    ranked_gene_list <- setNames(de_res$log2FoldChange[gene_rank], row.names(de_res)[gene_rank])

    oryzabase_GSEA_res_list = oryzabase_GSEA(ranked_gene_list, config$oryzabase_xlsx, config$oryzabase_sheet)

    writexl::write_xlsx(lapply(oryzabase_GSEA_res_list, as.data.frame), config$output$GSEA_xlsx)
    saveRDS(oryzabase_GSEA_res_list, config$output$GSEA_rds)
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

df_GSEA <- function(ranked_gene_list, onto_df, gene_col, term_col, name_col) {
    res <- clusterProfiler::GSEA(
        geneList = ranked_gene_list,
        TERM2GENE = onto_df[c(term_col, gene_col)],
        TERM2NAME = onto_df[c(term_col, name_col)],
        pvalueCutoff = 1,
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

oryzabase_GSEA <- function(ranked_gene_list, oryzabase_xlsx, sheet_list) {
    res_list <- list()
    for (sheet in sheet_list) {
        onto_df <- readxl::read_excel(oryzabase_xlsx, sheet=sheet)
        res <- df_GSEA(ranked_gene_list, onto_df, "GeneID", "OntoID", "Description")
        res_list[[sheet]] <- res
    }
    return(res_list)
}

main()

