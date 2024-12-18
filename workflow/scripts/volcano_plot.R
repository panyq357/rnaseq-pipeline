library(ggplot2)

config <- list(
    de_res = snakemake@input$de_res,
    p_column = snakemake@wildcards$p_column,
    p_th = as.numeric(snakemake@wildcards$p_th),
    l2fc_th = log2(as.numeric(snakemake@wildcards$fc_th)),
    name = snakemake@wildcards$name,
    outfile = snakemake@output$plot
)

main <- function() {

    de_res <- read.csv(config$de_res)

    de_res$color <- "none"

    is_up   <- de_res[[config$p_column]] < config$p_th & de_res$log2FoldChange >  config$l2fc_th
    is_down <- de_res[[config$p_column]] < config$p_th & de_res$log2FoldChange < -config$l2fc_th

    de_res$color[is_up] <- "up"
    de_res$color[is_down] <- "down"

    pdf(config$outfile, width=5, height=5)

    plot(
        ggplot(de_res) +
            geom_point(aes(x=log2FoldChange, y=-log10(.data[[config$p_column]]), color=color)) +
            theme_bw() +
            scale_color_manual(
                name = "",
                values = c("none"="grey", "down"="blue", "up"="red"),
                labels = c(
                    "none" = sprintf("none: %s", sum(de_res$color == "none")),
                    "down" = sprintf("down: %s", sum(de_res$color == "down")),
                    "up" = sprintf("up: %s", sum(de_res$color == "up"))
                )
            ) +
            geom_vline(xintercept=c(-config$l2fc_th, config$l2fc_th)) +
            geom_hline(yintercept=-log10(config$p_th)) +
            theme(legend.position="top") +
            labs(title = config$name, y = config$p_column)
    )

    dev.off()
}

main()
