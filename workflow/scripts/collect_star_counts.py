# description: collect counts in `*.ReadsPerGene.out.tab` files generated by STAR, concat them into a matrix.
# author: panyq357
# date: 2023-08-23

from pathlib import Path

import pandas as pd

count_matrix = pd.DataFrame()
for sample_id, file in zip(snakemake.params["sample_id_list"], snakemake.input):
    path = Path(file)
    df = pd.read_table(path, header=None, index_col=0)
    df = df.iloc[4:, ]  # Remove leading 4 statistical rows.
    sample_counts = df.iloc[:, 0]
    sample_counts.name = sample_id
    count_matrix = pd.concat([count_matrix, sample_counts], axis=1)

count_matrix.to_csv(snakemake.output[0])
