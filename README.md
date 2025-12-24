A pipeline for RNA-seq analysis.

## How to run this pipeline

First, Install follwing dependencies:

- [fastp](https://github.com/OpenGene/fastp)
- [STAR](https://github.com/alexdobin/STAR)
  - or [salmon](https://github.com/COMBINE-lab/salmon)
- [samtools](https://github.com/samtools/samtools)
- R packages:
  - from CRAN
    - `tidyverse`
    - `writexl`
  - from Bioconductor
    - `rtracklayer`
    - `GenomicRanges`
    - `tximport`
    - `ComplexHeatmap`
    - `DESeq2`
    - `clusterProfiler`

Dependencies can also be installed by running follwing scripts:

```
bash dependencies.bash
Rscript dependencies.R
```

Then run these to copy config templates.

```bash
cp workflow/Snakefile_template workflow/Snakefile
cp -r config_template config
```

Then filling config files in `config` directory.

Finally, run this.

```bash
snakemake --cores 20 --resources io=100
```
