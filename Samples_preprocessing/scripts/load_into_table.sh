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
	sseqid_id varchar, 
	sseqid_count int, 
	sseqid_proportion_sum real
	)
;

INSERT INTO Protein_abundance
	SELECT
		parse_filename(filename,true, 'system') AS Sample_id,
		sseqid_id as sseqid_id, 
		count_protein_keys AS sseqid_count, 
		sum_proportion AS sseqid_proportion_sum	 
	FROM
		read_parquet('data/grouped_parquets/*.parquet', compression = 'zstd', filename = true)
; 
EOF
