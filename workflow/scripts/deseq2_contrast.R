library(DESeq2)


config <- list(
    cts = snakemake@input$counts,
    levels = snakemake@params$levels,  # c("ref_group", "exp_group")
    sample_groups = snakemake@params$sample_groups,
    name = snakemake@wildcards$name,
    out_rds = snakemake@output$out_rds,
    out_csv = snakemake@output$out_csv
)


main <- function() {

    cts <- read.csv(config$cts, check.names=F, row.names=1)

    # Make coldata.
    sample_groups <- config$sample_groups[config$levels]
    group_vect <- c()
    for (group in names(sample_groups)) {
        group_vect <- c(group_vect, rep(group, length(sample_groups[[group]])))
    }
    coldata <- data.frame(Sample=unlist(sample_groups), Group=group_vect)
    row.names(coldata) <- coldata$Sample

    coldata$Group <- factor(coldata$Group, levels=config$levels)
    coldata <- coldata[order(coldata$Group),]
    cts <- cts[coldata$Sample]

    dds <- DESeqDataSetFromMatrix(cts, coldata, ~ Group)

    # Pre-filtering.
    dds <- dds[rowSums(counts(dds)) >= 10,]

    # DESeq.
    dds <- DESeq(dds)

    # Use lfcShrink to get a better Log2FoldChange for downstream GSEA.
    # colnames(coef(dds))[2] looks like: "Group_mutent_vs_wildtype".
    res <- lfcShrink(dds, colnames(coef(dds))[2], type="apeglm")

    res_df <- data.frame(res, check.names=F)
    
    # Combine normalized counts.
    nc <- counts(dds, normalized=TRUE)
    res_df <- cbind(res_df, nc[match(row.names(res_df), row.names(nc)),])

    # Sort by padj
    res_df <- res_df[order(res_df$padj, decreasing=FALSE),]

    ## Output DESeq results.
    saveRDS(res, config$out_rds)
    write.csv(res_df, config$out_csv)
}

main()
