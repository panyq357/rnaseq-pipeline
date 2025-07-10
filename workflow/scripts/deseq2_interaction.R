log <- file(snakemake@log[[1]], open="wt")
sink(log)  # STDOUT
sink(log, type = "message")  # STDERR


library(DESeq2)


cts <- read.csv(snakemake@input$counts, check.names = F, row.names = 1)
coldata <- read.csv(snakemake@input$coldata, check.names = F)

levels <- snakemake@params$levels

cols <- names(levels)

for (col in cols) {
  coldata <- subset(coldata, coldata[[col]] %in% snakemake@params$levels[[col]])
  coldata[[col]] <- factor(coldata[[col]], levels = snakemake@params$levels[[col]])
}

# Use first two col to do analysis.
coldata <- coldata[order(coldata[[cols[1]]], coldata[[cols[2]]]),]

cts <- cts[coldata$sample]

dds <- DESeqDataSetFromMatrix(cts, coldata, as.formula(sprintf("~ %s + %s + %s:%s", cols[1], cols[2], cols[1], cols[2])))

# Pre-filtering.
dds <- dds[rowSums(counts(dds)) >= 20,]

# DESeq.
dds <- DESeq(dds)

# Use lfcShrink to get a better Log2FoldChange for downstream GSEA.
# colnames(coef(dds))[2] looks like: "Group_mutent_vs_wildtype".
res <- lfcShrink(dds, colnames(coef(dds))[2], type="apeglm")

res_df <- data.frame(res, check.names=F)

# Combine normalized counts.
nc <- counts(dds, normalized=TRUE)

# Add gene annotation
gene_info <- readr::read_tsv(snakemake@input$gene_info, col_names=c("geneID", "symbol", "description"))
gene_info <- gene_info[match(row.names(res_df), gene_info$geneID),]

res_df <- cbind(res_df, nc[match(row.names(res_df), row.names(nc)),], gene_info[-1])

# Sort by padj
res_df <- res_df[order(res_df$pvalue, decreasing=FALSE),]

## Output DESeq results.
tibble::rownames_to_column(res_df, var = "geneID") |>
  readr::write_csv(snakemake@output[[1]])
