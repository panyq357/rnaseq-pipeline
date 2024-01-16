config <- list(
    de_res = snakemake@input$de_res,
    tair_10_anno = snakemake@input$tair_10_anno,
    de_res_with_anno = snakemake@output$de_res_with_anno
)

main <- function() {

    anno <- readr::read_tsv(config$tair_10_anno, col_names=c("ID", "Type", "Brief", "Description", "OtherAnnotation"))
    anno$ID <- sub("(AT[12345CM]G[0-9]{5}).*", "\\1", anno$ID)

    res <- read.csv(config$de_res, row.names=1, check.names=F)
    res <- merge(res, anno, by.x="row.names", by.y="ID")

    write.csv(res, config$de_res_with_anno, row.names=F)
}

main()
