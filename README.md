This is a pipeline for rice and arabidopsis RNA-seq analysis.

## How to run this pipeline

Step 1. Filling `config/reads2star-config.yaml` with STAR index, GTF, raw FASTQ paths.

Step 2. Filling `config/sample_check-config.yaml` with sample group information.

Step 3. Filling `config/de_analysis-config.yaml` with contrasts for DESeq2.

Step 4-1. For arabidopsis annotation, filling `config/at_anno-config.yaml` with anntation file path
(which can be downloaded from TAIR), and uncomment these line in `workflow/Snakefile`.

```
include: "rules/at_anno.smk"
```

```
# at_anno
expand("results/at_anno/{name}/{name}.deseq2_results_anno.csv", name=config["levels"]),
expand("results/at_anno/{name}/{name}.clusterProfiler.padj-0.01.fc-1.5.xlsx", name=config["levels"])
```

Else, for rice annotation, comment out arabidopsis lines above, filling `config/os_anno-config.yaml`
with annotation file path (which can be downloaded from RAP-DB) and oryzabase annotation xlsx path (see
<https://github.com/panyq357/onto-scripts>), and uncomment these lines in `workflow/Snakefile`.

```
include: "rules/os_anno.smk"
```

```
# os_anno
expand("results/os_anno/{name}/{name}.oryzabase-enricher.padj-0.01.fc-1.5.xlsx", name=config["levels"]),
expand("results/os_anno/{name}/{name}.kegg-enrich.padj-0.01.fc-1.5.csv", name=config["levels"])
```

Finally run this.

```bash
snakemake --cores 20
```

## Results explained

```
results/
├── at_anno                                                              # Annotation of results.
│   ├── MOCK_WT_vs_JA_WT
│   └── MOCK_WT_vs_MOCK_lhp1
├── collect_star_counts                                                  # Counts matrix.
│   └── counts.csv
├── deseq2_pipeline                                                      # DESeq2 results.
│   ├── MOCK_WT_vs_JA_WT.deseq2_results.csv
│   ├── MOCK_WT_vs_JA_WT.deseq2_results.rds
│   ├── MOCK_WT_vs_MOCK_lhp1.deseq2_results.csv
│   └── MOCK_WT_vs_MOCK_lhp1.deseq2_results.rds
├── fastp                                                                # raw FASTQ QC reports.
│   ├── JA_WT_1.html
│   ├── JA_WT_1.json
│   ├── JA_WT_2.html
│   ├── JA_WT_2.json
│   ├── JA_WT_3.html
│   ├── JA_WT_3.json
│   ├── MOCK_lhp1_1.html
│   ├── MOCK_lhp1_1.json
│   ├── MOCK_lhp1_2.html
│   ├── MOCK_lhp1_2.json
│   ├── MOCK_lhp1_3.html
│   ├── MOCK_lhp1_3.json
│   ├── MOCK_WT_1.html
│   ├── MOCK_WT_1.json
│   ├── MOCK_WT_2.html
│   ├── MOCK_WT_2.json
│   ├── MOCK_WT_3.html
│   └── MOCK_WT_3.json
├── sample_check                                                         # Sample check plots.
│   └── sample_dist_heatmap.svg
├── star_mapping                                                         # STAR mapping results.
│   ├── JA_WT_1.Aligned.sortedByCoord.out.bam
│   ├── JA_WT_1.Aligned.sortedByCoord.out.bam.bai
│   ├── JA_WT_1.Log.final.out
│   ├── JA_WT_1.Log.out
│   ├── JA_WT_1.Log.progress.out
│   ├── JA_WT_1.ReadsPerGene.out.tab
│   ├── JA_WT_1.SJ.out.tab
│   ├── JA_WT_2.Aligned.sortedByCoord.out.bam
│   ├── JA_WT_2.Aligned.sortedByCoord.out.bam.bai
│   ├── JA_WT_2.Log.final.out
│   ├── JA_WT_2.Log.out
│   ├── JA_WT_2.Log.progress.out
│   ├── JA_WT_2.ReadsPerGene.out.tab
│   ├── JA_WT_2.SJ.out.tab
│   ├── JA_WT_3.Aligned.sortedByCoord.out.bam
│   ├── JA_WT_3.Aligned.sortedByCoord.out.bam.bai
│   ├── JA_WT_3.Log.final.out
│   ├── JA_WT_3.Log.out
│   ├── JA_WT_3.Log.progress.out
│   ├── JA_WT_3.ReadsPerGene.out.tab
│   ├── JA_WT_3.SJ.out.tab
│   ├── MOCK_lhp1_1.Aligned.sortedByCoord.out.bam
│   ├── MOCK_lhp1_1.Aligned.sortedByCoord.out.bam.bai
│   ├── MOCK_lhp1_1.Log.final.out
│   ├── MOCK_lhp1_1.Log.out
│   ├── MOCK_lhp1_1.Log.progress.out
│   ├── MOCK_lhp1_1.ReadsPerGene.out.tab
│   ├── MOCK_lhp1_1.SJ.out.tab
│   ├── MOCK_lhp1_2.Aligned.sortedByCoord.out.bam
│   ├── MOCK_lhp1_2.Aligned.sortedByCoord.out.bam.bai
│   ├── MOCK_lhp1_2.Log.final.out
│   ├── MOCK_lhp1_2.Log.out
│   ├── MOCK_lhp1_2.Log.progress.out
│   ├── MOCK_lhp1_2.ReadsPerGene.out.tab
│   ├── MOCK_lhp1_2.SJ.out.tab
│   ├── MOCK_lhp1_3.Aligned.sortedByCoord.out.bam
│   ├── MOCK_lhp1_3.Aligned.sortedByCoord.out.bam.bai
│   ├── MOCK_lhp1_3.Log.final.out
│   ├── MOCK_lhp1_3.Log.out
│   ├── MOCK_lhp1_3.Log.progress.out
│   ├── MOCK_lhp1_3.ReadsPerGene.out.tab
│   ├── MOCK_lhp1_3.SJ.out.tab
│   ├── MOCK_WT_1.Aligned.sortedByCoord.out.bam
│   ├── MOCK_WT_1.Aligned.sortedByCoord.out.bam.bai
│   ├── MOCK_WT_1.Log.final.out
│   ├── MOCK_WT_1.Log.out
│   ├── MOCK_WT_1.Log.progress.out
│   ├── MOCK_WT_1.ReadsPerGene.out.tab
│   ├── MOCK_WT_1.SJ.out.tab
│   ├── MOCK_WT_2.Aligned.sortedByCoord.out.bam
│   ├── MOCK_WT_2.Aligned.sortedByCoord.out.bam.bai
│   ├── MOCK_WT_2.Log.final.out
│   ├── MOCK_WT_2.Log.out
│   ├── MOCK_WT_2.Log.progress.out
│   ├── MOCK_WT_2.ReadsPerGene.out.tab
│   ├── MOCK_WT_2.SJ.out.tab
│   ├── MOCK_WT_3.Aligned.sortedByCoord.out.bam
│   ├── MOCK_WT_3.Aligned.sortedByCoord.out.bam.bai
│   ├── MOCK_WT_3.Log.final.out
│   ├── MOCK_WT_3.Log.out
│   ├── MOCK_WT_3.Log.progress.out
│   ├── MOCK_WT_3.ReadsPerGene.out.tab
│   └── MOCK_WT_3.SJ.out.tab
└── volcano_plot                                                         # Volcano plots.
    ├── MOCK_WT_vs_JA_WT.padj-0.01.fc-1.5.volcano.pdf
    └── MOCK_WT_vs_MOCK_lhp1.padj-0.01.fc-1.5.volcano.pdf

```

