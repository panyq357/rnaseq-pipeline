rule deseq2_contrast:
    input:
        counts = "results/counts.csv",
        coldata = config["coldata"],
        gene_info = config["gene_info"]
    output:
        "results/de_analysis/contrast/{name}.deseq2_contrast.csv"
    log:
        "results/de_analysis/contrast/{name}.deseq2_contrast.log"
    params:
        levels = lambda w: config["de_analysis_jobs"]["contrast"][w.name]
    script:
        "../scripts/deseq2_contrast.R"


rule deseq2_interaction:
    input:
        counts = "results/counts.csv",
        coldata = config["coldata"],
        gene_info = config["gene_info"]
    output:
        "results/de_analysis/interaction/{name}.deseq2_interaction.csv"
    log:
        "results/de_analysis/interaction/{name}.deseq2_interaction.log"
    params:
        levels = lambda w: config["de_analysis_jobs"]["interaction"][w.name]
    script:
        "../scripts/deseq2_interaction.R"


rule volcano_plot:
    input:
        "results/de_analysis/{type}/{name}.deseq2_{type}.csv"
    output:
        "results/de_analysis/{type}/{name}.deseq2_{type}.{p_column}-{p_th}.fc-{fc_th}.volcano.svg"
    log:
        "results/de_analysis/{type}/{name}.deseq2_{type}.{p_column}-{p_th}.fc-{fc_th}.volcano.log"
    script:
        "../scripts/volcano_plot.R"


rule enricher:
    input:
        de_res = "results/de_analysis/{type}/{name}.deseq2_{type}.csv",
        enrichment_data = [dataset.values() for dataset in config["enrichment_data"].values()]
    output:
        xlsx = "results/de_analysis/{type}/{name}.deseq2_{type}.{p_column}-{p_th}.fc-{fc_th}.enricher.xlsx",
        pdf = "results/de_analysis/{type}/{name}.deseq2_{type}.{p_column}-{p_th}.fc-{fc_th}.enricher.pdf"
    params:
        enrichment_data = config["enrichment_data"]
    log:
        "results/de_analysis/{type}/{name}.deseq2_{type}.{p_column}-{p_th}.fc-{fc_th}.enricher.log"
    script:
        "../scripts/enricher.R"
