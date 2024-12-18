config <- list(

    # input
    de_res = snakemake@input$de_res,
    anno = snakemake@input$anno,  # https://www.arabidopsis.org/api/download-files/download?filePath=Genes/TAIR10_genome_release/TAIR10_functional_descriptions

    # output
    de_res_with_anno = snakemake@output$de_res_with_anno
)

main <- function() {

    anno <- readxl::read_excel(config$anno)

    res <- read.csv(config$de_res, row.names=1, check.names=F)

    res[c("Name", "Description")] <- anno[match(row.names(res), anno$Locus), c("Other Name(Type)", "Description")]

    write.csv(res, config$de_res_with_anno)
}

main()
