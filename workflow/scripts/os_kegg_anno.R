#' KEGG enrichment analysis
#'
#' @importFrom clusterProfiler  enrichKEGG
#' @importFrom biomaRt          useEnsemblGenomes getBM
#' 

config <- list(
    de_res = snakemake@input$de_res,
    padj_th = as.numeric(snakemake@wildcards$padj_th),
    l2fc_th = log2(as.numeric(snakemake@wildcards$fc_th)),
    out_list = list(
        enrich_csv = snakemake@output$enrich_csv,
        enrich_rds = snakemake@output$enrich_rds,
        gse_csv = snakemake@output$gse_csv,
        gse_rds = snakemake@output$gse_rds
    )
)


main <- function(){

    de_res <- read.csv(config$de_res, row.names=1)

    organism="dosa"  # Organism codename list: https://www.genome.jp/kegg/catalog/org_list.html

    de_res$tx <- os_gene_to_tx(row.names(de_res))

    universe <- de_res$tx
    gene <- subset(de_res, padj < config$padj_th & abs(log2FoldChange) > config$l2fc_th)$tx

    enrich_kegg_res <- clusterProfiler::enrichKEGG(
        gene = gene,
        universe = universe,
        organism = organism,
        pvalueCutoff = 1,
        qvalueCutoff = 1
    )

    write.csv(enrich_kegg_res, config$out_list$enrich_csv, row.names=F)
    saveRDS(enrich_kegg_res, config$out_list$enrich_rds)

    rank <- order(de_res$log2FoldChange, decreasing=T)
    ranked_gene <- setNames(de_res$log2FoldChange, de_res$tx)[rank]

    gse_kegg_res <- clusterProfiler::gseKEGG(
        geneList = ranked_gene,
        organism = organism,
        pvalueCutoff = 1
    )

    write.csv(gse_kegg_res, config$out_list$gse_csv, row.names=F)
    saveRDS(gse_kegg_res, config$out_list$gse_rds)
}

# KEGG use transcript ID as entry.
# For each locus, there's only one corresponding transcript ID.
# It seems that KEGG use the first, rather then the longest. (e.g. Os05g0591900)
ensembl_osativa <- NULL
os_gene_to_tx <- function(locus_id_vect) {

    if (is.null(ensembl_osativa)) {
        cat("Connecting to EnsemblGenomes osativa ... ")
        assign("ensembl_osativa", biomaRt::useEnsemblGenomes(biomart = "plants_mart", dataset = "osativa_eg_gene"), envir=globalenv())
        cat("Done.\n")
    }

    res_df = biomaRt::getBM(
        attributes=c("ensembl_gene_id", "ensembl_transcript_id"),
        filters="ensembl_gene_id",
        values=locus_id_vect,
        mart=ensembl_osativa
    )

    # Sort by gene ID and then transcript ID.
    res_df = res_df[with(res_df, order(ensembl_gene_id, ensembl_transcript_id)),]

    # Keep the first transcript ID for each gene.
    res_df = res_df[!duplicated(res_df$ensembl_gene_id),]

    # match the order of original locus_id_vect.
    tx_id_vect <- with(res_df, ensembl_transcript_id[match(locus_id_vect, ensembl_gene_id)])

    return(tx_id_vect)
}


main()

