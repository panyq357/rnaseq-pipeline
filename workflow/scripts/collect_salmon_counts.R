gtf <- rtracklayer::import(snakemake@input$gtf)

tx2gene <- GenomicRanges::mcols(gtf)[c("transcript_id", "gene_id")] |>
  as.data.frame() |>
  dplyr::distinct() |>
  dplyr::filter(!is.na(transcript_id))

txi <- tximport::tximport(snakemake@input$quant_files, type = "salmon", tx2gene = tx2gene)

colnames(txi$counts) <- snakemake@params$sample_id_list

write.csv(txi$counts, snakemake@output[[1]])
