library(clusterProfiler)


config <- list(
    de_res = snakemake@input$de_res,
    pathway_to_gene = snakemake@input$pathway_to_gene,
    pathway_to_description = snakemake@input$pathway_to_description,
    out_xlsx = snakemake@output$out_xlsx,
    out_rds = snakemake@output$out_rds,
    out_pdf = snakemake@output$out_pdf
)


main <- function(){

    de_res <- read.csv(config$de_res, row.names=1)
    pathway_to_gene <- read.table(config$pathway_to_gene, col.names=c("Pathway", "Gene"), sep="\t")
    pathway_to_description <- read.table(config$pathway_to_description, col.names=c("Pathway", "Description"), sep="\t")

    # Rank all genes by log2FoldChange for GSEA.
    gene_rank <- order(de_res$log2FoldChange, decreasing=T)
    ranked_gene_list <- setNames(de_res$log2FoldChange[gene_rank], row.names(de_res)[gene_rank])

    gsea_res <- clusterProfiler::GSEA(
        geneList = ranked_gene_list,
        TERM2GENE = pathway_to_gene,
        TERM2NAME = pathway_to_description,
        pvalueCutoff = 1
    )

    writexl::write_xlsx(list(KEGG_GSEA=as.data.frame(gsea_res)), config$out_xlsx)
    saveRDS(gsea_res, config$out_rds)
    
    pdf(config$out_pdf, width=10, height=10)
        print(
            enrichplot::ridgeplot(gsea_res, showCategory=20) +
            ggplot2::labs(title="KEGG GSEA", x="log2(Fold Change)")
        )
    dev.off()

}

main()

