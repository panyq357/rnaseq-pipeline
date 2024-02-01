config <- list(
    de_res = snakemake@input$de_res,
    dosa_id_converter = snakemake@input$dosa_id_converter,
    out_csv = snakemake@output$out_csv,
    out_rds = snakemake@output$out_rds,
    out_pdf = snakemake@output$out_pdf
)

main <- function(){

    de_res <- read.csv(config$de_res, row.names=1)

    organism <- "dosa"  # Organism codename list: https://www.genome.jp/kegg/catalog/org_list.html

    id_converter <- read.table(config$dosa_id_converter)
    id_converter <- setNames(id_converter[[1]], id_converter[[2]])

    de_res$tx <- id_converter[row.names(de_res)]

    rank <- order(de_res$log2FoldChange, decreasing=T)
    ranked_gene <- setNames(de_res$log2FoldChange, de_res$tx)[rank]

    gse_kegg_res <- clusterProfiler::gseKEGG(
        geneList = ranked_gene,
        organism = organism,
        pvalueCutoff = 1
    )

    write.csv(gse_kegg_res, config$out_csv, row.names=F)
    saveRDS(gse_kegg_res, config$out_rds)

    pdf(config$out_pdf, width=10, height=30)
    print(
        enrichplot::ridgeplot(gse_kegg_res) +
        ggplot2::labs(x="log2(Fold Change)")
    )
    dev.off()
}


main()

