library(clusterProfiler)

config <- list(
    de_res = snakemake@input$de_res,

    padj_th = as.numeric(snakemake@wildcards$padj_th),
    abs_l2fc_th = log2(as.numeric(snakemake@wildcards$fc_th)),

    xlsx = snakemake@output$xlsx,
    rds = snakemake@output$rds
)


main <- function() {

    de_res <- read.csv(config$de_res, row.names=1)

    universe <- row.names(de_res)

    gene_list <- list(
        Up = subset(de_res, padj < config$padj_th & log2FoldChange > config$abs_l2fc_th) |> row.names(),
        Down = subset(de_res, padj < config$padj_th & log2FoldChange < config$abs_l2fc_th) |> row.names()
    )
    gene_list$All = c(gene_list$Up, gene_list$Down)

    # Genes overlapped with diff peaks.
    enrich_go_res_list <- lapply(gene_list, at_enrich_go)
    enrich_go_res_list <- unlist(enrich_go_res_list, recursive=FALSE)
    enrich_kegg_res_list <- lapply(gene_list, at_enrich_kegg)

    l2fc <- setNames(de_res$log2FoldChange, row.names(de_res)) |> sort(decreasing=TRUE)

    gse_go_res_list <- at_gse_go(l2fc)
    gse_kegg_res <- at_gse_kegg(l2fc)

    names(enrich_go_res_list) <- sprintf("EnrichGO.%s", names(enrich_go_res_list))
    names(enrich_kegg_res_list) <- sprintf("EnrichEKGG.%s", names(enrich_kegg_res_list))
    names(gse_go_res_list) <- sprintf("GSEA.GO.%s", names(gse_go_res_list))

    all_res <- c(enrich_go_res_list, enrich_kegg_res_list, gse_go_res_list, list(GSEA.KEGG = gse_kegg_res))

    writexl::write_xlsx(lapply(all_res, as.data.frame), config$xlsx)
    saveRDS(all_res, config$rds)

}

at_enrich_go <- function(gene, universe=NULL) {
    ego <- list()
    for (ont in c("BP", "MF", "CC")) {
        ego[[ont]] <- clusterProfiler::enrichGO(
            gene = gene,
            universe = universe,
            OrgDb = org.At.tair.db::org.At.tair.db,
            keyType = 'TAIR',
            ont = ont,
            pvalueCutoff = 1,
            qvalueCutoff = 1
        )
    }
    return(ego)
}

at_enrich_kegg <- function(gene, universe=NULL) {
    kk <- clusterProfiler::enrichKEGG(
        gene = gene,
        universe = universe,
        organism = 'ath',
        pvalueCutoff = 1,
        qvalueCutoff = 1
    )
    return(kk)
}

at_gse_go <- function(l2fc) {
    gsego <- list()
    for (ont in c("BP", "MF", "CC")) {
        gsego[[ont]] <- clusterProfiler::gseGO(
            geneList = l2fc,
            OrgDb = org.At.tair.db::org.At.tair.db,
            keyType = 'TAIR',
            ont = ont,
            pvalueCutoff = 1
        )
    }
    return(gsego)
}

at_gse_kegg <- function(l2fc) {
    gsekk <- clusterProfiler::gseKEGG(
        geneList = l2fc,
        organism = 'ath',
        pvalueCutoff = 1
    )
    return(gsekk)
}

main()
