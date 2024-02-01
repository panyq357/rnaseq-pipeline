config <- list(
    de_res = snakemake@input$de_res,
    dosa_id_converter = snakemake@input$dosa_id_converter,
    padj_th = as.numeric(snakemake@wildcards$padj_th),
    abs_l2fc_th = log2(as.numeric(snakemake@wildcards$fc_th)),
    out_xlsx = snakemake@output$out_xlsx,
    out_rds = snakemake@output$out_rds,
    out_pdf = snakemake@output$out_pdf
)


main <- function(){

    de_res <- read.csv(config$de_res, row.names=1)

    organism="dosa"  # Organism codename list: https://www.genome.jp/kegg/catalog/org_list.html

    id_converter <- read.table(config$dosa_id_converter)
    id_converter <- setNames(id_converter[[1]], id_converter[[2]])

    de_res$tx <- id_converter[row.names(de_res)]
    gene <- list(
        All  = (de_res |> subset(padj < config$padj_th & abs(log2FoldChange) >  config$abs_l2fc_th))$tx |> na.omit(),
        Up   = (de_res |> subset(padj < config$padj_th & log2FoldChange      >  config$abs_l2fc_th))$tx |> na.omit(),
        Down = (de_res |> subset(padj < config$padj_th & log2FoldChange      < -config$abs_l2fc_th))$tx |> na.omit() 
    )

    universe <- na.omit(de_res$tx)

    enrich_kegg_res_list <- lapply(gene, function(gene) {
        clusterProfiler::enrichKEGG(
            gene = gene,
            universe = universe,
            organism = organism,
            pvalueCutoff = 1,
            qvalueCutoff = 1
        )
    })

    out_df_list <- c(
        lapply(enrich_kegg_res_list, as.data.frame),
        setNames(lapply(gene, function(x) data.frame(Gene=x)), sprintf("GeneUsed.%s", names(gene)))
    )

    writexl::write_xlsx(out_df_list, config$out_xlsx)
    saveRDS(enrich_kegg_res_list, config$out_rds)

    pdf(config$out_pdf)
    for (name in names(enrich_kegg_res_list)) {
        print(
            enrichplot::dotplot(enrich_kegg_res_list[[name]]) + ggplot2::labs(title=name)
        )
    }
    dev.off()

}


main()

