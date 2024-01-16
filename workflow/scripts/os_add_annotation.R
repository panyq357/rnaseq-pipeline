
config <- list(
    de_res = snakemake@input$de_res,
    anno = snakemake@input$anno,
    de_res_with_anno = snakemake@output$de_res_with_anno,
    anno_reader = readr::read_tsv,
    anno_idx = "Locus_ID",
    anno_cols = list(
        "Symbol" = "Oryzabase Gene Symbol Synonym(s)",
        "Description" = "Description"
    )
)

main <- function() {

    anno <- config$anno_reader(config$anno)

    res <- read.csv(config$de_res, row.names=1, check.names=F)
    idx <- match(row.names(res), anno[[config$anno_idx]])
    for (col_name in names(config$anno_cols)) {
        res[[col_name]] <- anno[[config$anno_cols[[col_name]]]][idx]
    }
    write.csv(res, config$de_res_with_anno)
}

main()
