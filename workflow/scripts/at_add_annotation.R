config <- list(

    # input
    de_res = snakemake@input$de_res,
    anno = snakemake@input$anno,  # https://www.arabidopsis.org/api/download-files/download?filePath=Genes/TAIR10_genome_release/TAIR10_functional_descriptions

    # output
    de_res_with_anno = snakemake@output$de_res_with_anno
)

main <- function() {

    anno <- readr::read_tsv(config$anno)

    anno$Model_name <- sub("([^.]+)\\..*", "\\1", anno$Model_name)

    res <- read.csv(config$de_res, row.names=1, check.names=F)

    res[c("Short_description", "Curator_summary")] <- anno[match(row.names(res), anno$Model_name), c("Short_description", "Curator_summary")]

    write.csv(res, config$de_res_with_anno)
}

main()
