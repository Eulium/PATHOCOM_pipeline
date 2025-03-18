#!/bin/bash
# Loads all samples protein counst into a single DuckDB table for easy query
duckdb "${snakemake_output[0]}"  << EOF
.timer on
.echo off
.bail off
SET threads = ${snakemake[threads]}; 
CREATE TABLE Protein_abundance
	(
	Sample_id varchar, 
	sseqid varchar, 
	sseqid_count int, 
	sseqid_proportion_sum real
	)
;

INSERT INTO Protein_abundance
	SELECT
		sample AS Sample_id,
		sseqid as sseqid_id, 
		count_proteins AS sseqid_count, 
		sum_proportion AS sseqid_proportion_sum	 
	FROM
		read_parquet('data/grouped_parquets/sample=*/data_*.parquet"', compression = 'zstd', hive_partitioning = true)
; 
EOF
