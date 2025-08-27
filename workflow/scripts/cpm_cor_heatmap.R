counts <- read.csv(snakemake@input[[1]], row.names = 1, check.names = FALSE)

log2_cpm_plus_1 <- log2(apply(counts, 2, function(x) { x / sum(x) * 1E6 }) + 1)

log2_cpm_plus_1 <- log2_cpm_plus_1[order(rowSums(log2_cpm_plus_1), decreasing=T),][1:(nrow(log2_cpm_plus_1)*0.9),]

svg(snakemake@output[[1]])
sample_dist_mat <- as.matrix(dist(t(log2_cpm_plus_1)))
print(
  ComplexHeatmap::Heatmap(
    sample_dist_mat,
    column_title = "Euclidean Distance of Samples", name = " ",
    col = circlize::colorRamp2(range(sample_dist_mat), c("#1772b4", "#d8e3e7")),
  )
)
dev.off()
