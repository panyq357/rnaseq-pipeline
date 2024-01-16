library(ComplexHeatmap)

config <- list(
    counts = snakemake@input$counts,
    sample_groups = snakemake@params$sample_groups,
    sample_dist_heatmap = snakemake@output$sample_dist_heatmap
)

main <- function() {

    coldata <- data.frame(
        Group=rep(names(config$sample_groups), sapply(config$sample_groups, length)),
        Sample=unlist(config$sample_groups)
    )

    counts <- read.csv(config$counts, row.names=1, check.names=F)

    log2_cpm_plus_1 <- log2(apply(counts, 2, function(x) { x / sum(x) * 1E6 }) + 1)

    svg(config$sample_dist_heatmap)
        sample_dist_mat <- as.matrix(dist(t(log2_cpm_plus_1)))
        print(Heatmap(
            sample_dist_mat,
            column_title = "Euclidean Distance of Samples", name = " ",
            col = circlize::colorRamp2(range(sample_dist_mat), c("#1772b4","#d8e3e7")),
        ))
    dev.off()
}

main()
