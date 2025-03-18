#!/bin/bash
# Group proteins in sample and count occurances
# Write to new parquet file

duckdb "${snakemake_output[1]}"  << EOF
SET threads = ${snakemake[threads]}; 
COPY (
    SELECT
        sseqid,
        count(sseqid) as count_proteins,
        sum(proportion) AS sum_proportion
    FROM 
        read_parquet('${snakemake_input[sample_topx]}/*.parquet')
    GROUP BY
        sseqid
    )
TO  "${snakemake_output[sample_grouped]}" (FORMAT PARQUET, COMPRESSION 'zstd', PER_THREAD_OUTPUT false)
EOF
