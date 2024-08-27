library(clusterProfiler)


config <- list(

    de_res = snakemake@input$de_res,
    pathway_to_gene = snakemake@input$pathway_to_gene,
    pathway_to_description = snakemake@input$pathway_to_description,

    padj_th = as.numeric(snakemake@wildcards$padj_th),
    abs_l2fc_th = log2(as.numeric(snakemake@wildcards$fc_th)),

    out_xlsx = snakemake@output$out_xlsx,
    out_rds = snakemake@output$out_rds,
    out_pdf = snakemake@output$out_pdf

)


main <- function(){

    de_res <- read.csv(config$de_res, row.names=1)
    pathway_to_gene <- read.table(config$pathway_to_gene, col.names=c("Pathway", "Gene"), sep="\t")
    pathway_to_description <- read.table(config$pathway_to_description, col.names=c("Pathway", "Description"), sep="\t")

    universe <- row.names(de_res)
    gene <- list(
        All  = de_res |> subset(padj < config$padj_th & abs(log2FoldChange) >  config$abs_l2fc_th) |> row.names(),
        Up   = de_res |> subset(padj < config$padj_th & log2FoldChange      >  config$abs_l2fc_th) |> row.names(),
        Down = de_res |> subset(padj < config$padj_th & log2FoldChange      < -config$abs_l2fc_th) |> row.names()
    )

    res_list <- lapply(gene, function(gene) {
        clusterProfiler::enricher(
            gene = gene,
            universe = universe,
            TERM2GENE = pathway_to_gene,
            TERM2NAME = pathway_to_description,
            pvalueCutoff = 1,
            qvalueCutoff = 1
        )
    })

    out_df_list <- c(
        lapply(res_list, as.data.frame),
        setNames(lapply(gene, function(x) data.frame(Gene=x)), sprintf("GeneUsed.%s", names(gene)))
    )
    writexl::write_xlsx(out_df_list, config$out_xlsx)
    saveRDS(res_list, config$out_rds)

    pdf(config$out_pdf)
        for (name in names(res_list)) {
            print(
                enrichplot::dotplot(res_list[[name]]) +
                ggplot2::labs(title=name)
            )
        }
    dev.off()

}

main()

