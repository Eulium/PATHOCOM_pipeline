#!/bin/bash
# Group proteins in sample and count occurances
# Write to new parquet file

duckdb "${snakemake_output[1]}"  << EOF
SET threads = ${snakemake[threads]}; 
COPY (
    SELECT
        sseqid_id,
	count(sseqid_id) as count_protein_keys,
	sum(proportion) AS sum_proportion
    FROM 
        read_parquet('${snakemake_input[0]}/*.parquet')
    GROUP BY
        sseqid_id
    )
TO  "${snakemake_output[0]}" (FORMAT PARQUET, COMPRESSION 'zstd', PER_THREAD_OUTPUT false)
EOF
