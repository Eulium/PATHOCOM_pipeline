#!/bin/bash


echo ${snakemake_input[protein_abundance]}
duckdb ${snakemake_output[db]} << EOF 
-- SUMMARY:
-- SQL code use to pre calc distribtuions of GO terms and proteins clusters from hits.
-- Creats pivoted tables with go terms, one based on assembled data and one based on reads data.
-- Pivoted tables based on all observed protein clusters are not viable due to the large number of clusters. 
-- Therefore we first clalc distribution and then filter for the most significant (X% of all counts) clusters. 

-- READ DATA
.timer on
.echo on
.bail on 
CREATE TABLE IF NOT EXISTS protein_abundance AS
    SELECT
        *
    FROM
        read_parquet('data_top90/protein_abundance/type=*/isolate=*/data_*.parquet', hive_partitioning = true)
; 

CREATE TABLE IF NOT EXISTS go_abundance AS
    SELECT 
        *
    FROM
        read_parquet('data_top90/go_present/type=*/isolate=*/data_*.parquet', hive_partitioning = true)
;

-- NCBI MAPPING 
CREATE TABLE IF NOT EXISTS NCBI_mapping AS
    SELECT
        *
        REPLACE (string_split(GO, ';') AS GO)
        RENAME (UniRef90 AS sseqid)
    FROM
        read_csv("/ptmp/pboppert/Pathocom/diamond-postprocess-playground/input/idmapping_selected.tab.gz",
                header = false,
                names = [UniProtKB_AC, UniProtKB_ID, GeneID, RefSeq, GI, PDB, GO, UniRef100, UniRef90, UniRef50, UniParc, PIR, NCBI_taxon, MIM, UniGene, PubMed, EMBL, EMBL_CDS, Ensembl, Ensembl_TRS, Ensembl_PRO, Additional_PubMed],
                delim = '\t',
                compression = 'gzip')
;


-- VARIABLES:
--Protein clusters total counts
SET VARIABLE total_assembled_counts = 
(
    SELECT 
        sum(sum_proportion)
    FROM
        protein_abundance
    WHERE
        type = 'assembled'
)
;

SET VARIABLE total_reads_counts = 
(
    SELECT 
        sum(sum_proportion)
    FROM
        protein_abundance
    WHERE
        type = 'reads'
)
;

SET VARIABLE total_assembled_counts_go = 
(
    SELECT 
        sum(sum_proportion)
    FROM
        go_abundance
    WHERE
        type = 'assembled'
)
;

SET VARIABLE total_reads_counts_go = 
(
    SELECT 
        sum(sum_proportion)
    FROM
        go_abundance
    WHERE
        type = 'reads'
)
;


-- CLUSTER COUNTS
-- Total counts of protein clusters 
CREATE TABLE IF NOT EXISTS protein_clusters_counts AS
SELECT 
    sseqid,
    sum_proportion_assembled, 
    sum_proportion_reads
FROM
    (SELECT 
        sseqid,
        sum(assembled.sum_proportion) AS sum_proportion_assembled
    FROM
        protein_abundance AS assembled
    WHERE
        type = 'assembled'
    GROUP BY
        sseqid
    ) AS assembled
FULL JOIN
    (SELECT 
        sseqid,
        sum(reads.sum_proportion) AS sum_proportion_reads
    FROM
        protein_abundance AS reads
    WHERE
        type = 'reads'
    GROUP BY
        sseqid
    ) AS reads
USING
    (sseqid)
;

-- Distribution of protein clusters counts based on assmbled data
CREATE TABLE IF NOT EXISTS protein_clusters_dist_assembled AS
SELECT 
    sseqid,
    sum_proportion_assembled AS sum_proportion,
    sum(sum_proportion) over cumulative_counts AS cumulative_sum, 
    cumulative_sum / getvariable('total_assembled_counts') AS cumulative_sum_perc
FROM
    (
    SELECT
        sseqid, 
        sum_proportion_assembled
    FROM
        protein_clusters_counts
    WHERE
        sum_proportion_assembled is not NULL
    ORDER BY
        sum_proportion_assembled DESC
    )
WINDOW 
    cumulative_counts AS (ROWS BETWEEN UNBOUNDED PRECEDING AND CURRENT ROW)
;

-- Distribution of protein clusters counts based on assmbled data
CREATE TABLE IF NOT EXISTS protein_clusters_dist_reads AS
SELECT 
    sseqid,
    sum_proportion_reads AS sum_proportion,
    sum(sum_proportion) over cumulative_counts AS cumulative_sum, 
    cumulative_sum / getvariable('total_reads_counts') AS cumulative_sum_perc
FROM
    (
    SELECT
        sseqid, 
        sum_proportion_reads
    FROM
        protein_clusters_counts
    WHERE
        sum_proportion_reads is not NULL
    ORDER BY
        sum_proportion_reads DESC
    )
WINDOW 
    cumulative_counts AS (ROWS BETWEEN UNBOUNDED PRECEDING AND CURRENT ROW)
;


-- GO TERMS
CREATE TABLE IF NOT EXISTS go_counts AS
SELECT 
    go,
    sum_proportion_assembled, 
    sum_proportion_reads
FROM
    (SELECT 
        go,
        sum(assembled.sum_proportion) AS sum_proportion_assembled
    FROM
        go_abundance AS assembled
    WHERE
        type = 'assembled'
    GROUP BY
        go
    ) AS assembled
FULL JOIN
    (SELECT 
        go,
        sum(reads.sum_proportion) AS sum_proportion_reads
    FROM
        go_abundance AS reads
    WHERE
        type = 'reads'
    GROUP BY
        go
    ) AS reads
USING
    (go)
;

-- Distribution of GO terms counts based on assembled data
CREATE TABLE IF NOT EXISTS go_dist_assembled AS
SELECT 
    go,
    sum_proportion_assembled AS sum_proportion,
    sum(sum_proportion) over cumulative_counts AS cumulative_sum, 
    cumulative_sum / getvariable('total_assembled_counts_go') AS cumulative_sum_perc
FROM
    (
    SELECT
        go, 
        sum_proportion_assembled
    FROM
        go_counts
    WHERE
        sum_proportion_assembled is not NULL
    ORDER BY
        sum_proportion_assembled DESC
    )
WINDOW 
    cumulative_counts AS (ROWS BETWEEN UNBOUNDED PRECEDING AND CURRENT ROW)
;

-- Distribution of GO terms counts based on read data
CREATE TABLE IF NOT EXISTS go_dist_reads AS
SELECT 
    go,
    sum_proportion_reads AS sum_proportion,
    sum(sum_proportion) over cumulative_counts AS cumulative_sum, 
    cumulative_sum / getvariable('total_reads_counts_go') AS cumulative_sum_perc
FROM
    (
    SELECT
        go, 
        sum_proportion_reads
    FROM
        go_counts
    WHERE
        sum_proportion_reads is not NULL
    ORDER BY
        sum_proportion_reads DESC
    )
WINDOW 
    cumulative_counts AS (ROWS BETWEEN UNBOUNDED PRECEDING AND CURRENT ROW)
;

CREATE TABLE IF NOT EXISTS protein_clusters_counts_UNREF50 AS
    SELECT 
        NCBI_mapping.UniRef50,
        sum(sum_proportion_assembled) AS sum_proportion_assembled,
        sum(sum_proportion_reads) AS sum_proportion_reads
    FROM
        (SELECT 
            sseqid,
            sum(assembled.sum_proportion) AS sum_proportion_assembled
        FROM
            protein_abundance AS assembled
        WHERE
            type = 'assembled'
        GROUP BY
            sseqid
        ) AS assembled
    FULL JOIN
        (SELECT 
            sseqid,
            sum(reads.sum_proportion) AS sum_proportion_reads
        FROM
            protein_abundance AS reads
        WHERE
            type = 'reads'
        GROUP BY
            sseqid
        ) AS reads
    USING
        (sseqid)
    JOIN 
        NCBI_mapping
    ON
        sseqid = NCBI_mapping.sseqid
    GROUP BY 
        NCBI_mapping.UniRef50
;

SET VARIABLE total_assembled_counts_UNREF50 = 
(
    SELECT 
        sum(sum_proportion_assembled)
    FROM
        protein_clusters_counts_UNREF50
)
;

SET VARIABLE total_reads_counts_UNREF50 = 
(
    SELECT 
        sum(sum_proportion_reads)
    FROM
        protein_clusters_counts_UNREF50
)
;

.print 'Exporting'
-- EXPORT data

SET pivot_limit=1000000; 

-- TEST WITH UNREF50 clusters assmbled
COPY 
(
    SELECT 
        UniRef50 AS sseqid,
        sum_proportion_assembled AS sum_proportion,
        sum(sum_proportion) over cumulative_counts AS cumulative_sum, 
        cumulative_sum / getvariable('total_assembled_counts_UNREF50') AS cumulative_sum_perc
    FROM
        (
        SELECT
            protein_clusters_counts_UNREF50.UniRef50,
            sum_proportion_assembled
        FROM
            protein_clusters_counts_UNREF50
        WHERE
            sum_proportion_assembled is not NULL
        ORDER BY
            sum_proportion_assembled DESC
        )
    WINDOW 
        cumulative_counts AS (ROWS BETWEEN UNBOUNDED PRECEDING AND CURRENT ROW)
)
TO 
    "data_top90/protein_clusters_dist_UNREF50_assembled.pq"
(FORMAT parquet)
;

-- TEST WITH UNREF50 clusters reads
COPY 
(
    SELECT 
        UniRef50 AS sseqid,
        sum_proportion_reads AS sum_proportion,
        sum(sum_proportion) over cumulative_counts AS cumulative_sum, 
        cumulative_sum / getvariable('total_reads_counts_UNREF50') AS cumulative_sum_perc
    FROM
        (
        SELECT
            protein_clusters_counts_UNREF50.UniRef50,
            sum_proportion_reads
        FROM
            protein_clusters_counts_UNREF50
        WHERE
            sum_proportion_reads is not NULL
        ORDER BY
            sum_proportion_reads DESC
        )
    WINDOW 
        cumulative_counts AS (ROWS BETWEEN UNBOUNDED PRECEDING AND CURRENT ROW)
)
TO 
    "data_top90/protein_clusters_dist_UNREF50_reads.pq"
(FORMAT parquet)
;


-- -- UNref50 Cluster heatmap assmbled
-- COPY 
-- (
--     PIVOT 
--     (
--         SELECT
--             isolate,
--             NCBI_mapping.UniRef50 as UniRef50,
--             sum(assembled.sum_proportion) AS sum_proportion_assembled
--         FROM
--             protein_abundance AS assembled
--         JOIN
--             NCBI_mapping
--         ON  
--             assembled.sseqid = NCBI_mapping.sseqid
--         WHERE
--             type = 'assembled'
--         GROUP BY
--             isolate,NCBI_mapping.UniRef50
--     )
--     ON 
--         UniRef50
--     USING
--         sum(sum_proportion_assembled)
--     ORDER BY
--         isolate
-- )
-- TO 
--     "data/assmbled_cluster_matrix_UNREF50.pq"
-- (FORMAT PARQUET)
-- ;


-- COPY 
-- (
--     PIVOT 
--     (
--         SELECT
--             isolate,
--             NCBI_mapping.UniRef50 as UniRef50,
--             sum(assembled.sum_proportion) AS sum_proportion_assembled
--         FROM
--             protein_abundance AS assembled
--         JOIN
--             NCBI_mapping
--         ON  
--             assembled.sseqid = NCBI_mapping.sseqid
--         WHERE
--             type = 'reads'
--         GROUP BY
--             isolate,NCBI_mapping.UniRef50
--     )
--     ON 
--         UniRef50
--     USING
--         sum(sum_proportion_assembled)
--     ORDER BY
--         isolate
-- )
-- TO 
--     "data/reads_cluster_matrix_UNREF50.pq"
-- (FORMAT PARQUET)
;


-- GO terms matrix on assebmled data
COPY 
(
    PIVOT 
    (
        SELECT
            isolate,
            go,
            go_abundance.sum_proportion as sum_proportion
        FROM
            go_abundance
        JOIN
            go_dist_assembled
        USING
            (go)
        WHERE
            go_dist_assembled.cumulative_sum_perc < 0.95 and go_abundance.type = 'assembled'
    )
    ON 
        go
    USING
        sum(sum_proportion)
    ORDER BY
        isolate
)
TO 
    '${snakemake_output[assembled_matrix]}'
(FORMAT PARQUET)
;


-- GO terms matrix on read data
-- Only the most significant (0.95) GO terms are included in the matrix
COPY
(
    PIVOT 
    (
    SELECT
        isolate,
        go,
        go_abundance.sum_proportion as sum_proportion
    FROM
        go_abundance
    JOIN
        go_dist_reads
    USING
        (go)
    WHERE
        go_dist_reads.cumulative_sum_perc < 0.95 and go_abundance.type = 'reads'
    )
    ON 
        go
    USING
        sum(sum_proportion)
    ORDER BY
        isolate
)
TO 
    '${snakemake_output[read_matrix]}'
(FORMAT PARQUET)
;


COPY (
    SELECT 
        *
    FROM
        protein_clusters_dist_assembled 
)
TO 
    '${snakemake_output[assembled_protein_distribution]}'
(FORMAT PARQUET)
;

COPY (
    SELECT 
        *
    FROM
        protein_clusters_dist_reads 
)
TO 
    '${snakemake_output[reads_protein_distribution]}'
(FORMAT PARQUET)
;

COPY (
    SELECT 
        *
    FROM
        go_dist_assembled 
)
TO 
    '${snakemake_output[assembled_go_distribution]}'
(FORMAT PARQUET)
;   

COPY (
    SELECT 
        *
    FROM
        go_dist_reads 
)
TO 
    '${snakemake_output[reads_go_distribution]}'
(FORMAT PARQUET)
;
EOF
