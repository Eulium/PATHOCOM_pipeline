#!/bin/bash
# Apply the DIAMOND top_x setting to a sample set
# Final result does not contain unmapped reads

duckdb "${snakemake_output[1]}"  << EOF
.timer on
.echo off
.bail on
SET threads = ${snakemake[threads]}; 

-- Table containing reads
.print 'LOAD DATA'
CREATE OR REPLACE TEMP TABLE inserts (
    sample_id varchar,
    qseqid_id varchar,
    qseqid_key int,
    sseqid_id varchar,
    pident REAL,
    length SMALLINT,
    mismatch SMALLINT,
    gapopen SMALLINT,
    qstart SMALLINT,
    qend SMALLINT,
    sstart INT,
    send INT,
    evalue REAL,
    bitscore REAL,
    mapped BOOL)
;

-- Insert mapped reads into table
INSERT INTO inserts BY NAME
    (
    SELECT
        'SID' AS sample_id,
        qseqid_id,
        qseqid_key,
        sseqid_id,
        pident,
        length,
        mismatch,
        gapopen,
        qstart,
        qend,
        sstart,
        send,
        evalue,
        bitscore,
        mapped
    FROM
        read_parquet('${snakemake_input[0]}/*.parquet', compression = 'zstd', filename = true)
    WHERE mapped = true
    )
;


-- Calulates the percentage (perc) of each entry's bitscore with respect to the reads highest score
.print 'CREATE tmp table TOP_PERC'
CREATE TEMP TABLE
    TOP_PERC AS
    (
    WITH TOP_SCORE AS                                          -- Per read highest bitscore
        (                                    
        SELECT
            qseqid_key,
            max(bitscore) AS top_score                 
        FROM
            inserts
        GROUP BY
            qseqid_key
        )                               
    SELECT
        inserts.qseqid_key AS qseqid_key,                      -- Query and Protein keys needed to enable later mapping of
        inserts.sseqid_id AS sseqid_id,                        -- perc to propper read/protein entry
        inserts.bitscore / TOP_SCORE.top_score * 100 AS perc   -- Percentage of hit's A bitscore in reads highest bitscore B
        
    FROM
        inserts
    JOIN
        TOP_SCORE
    ON
        TOP_SCORE.qseqid_key = inserts.qseqid_key               -- JOINS read specific top_score to each entry, via read key
    ) 
;

-- Mapps the calculated percentage (perc) to each entry via query and subject id
-- Removes all rows with perc smaller than X  
.print 'Inserts top percentage, filter low perc' 
.print ${snakemake_params[topx]} 
CREATE OR REPLACE TEMP TABLE inserts AS 
    SELECT
        reads.sample_id,
        reads.qseqid_key,
        reads.qseqid_id,
        reads.sseqid_id,
        reads.pident,
        reads.length,
        reads.mismatch,
        reads.gapopen,
        reads.qstart,
        reads.qend,
        reads.sstart,
        reads.send,
        reads.evalue,
        reads.bitscore,
        TOP_PERC.perc AS perc
    FROM
        inserts AS reads
    JOIN
        TOP_PERC
    ON
        reads.qseqid_key = TOP_PERC.qseqid_key                
        AND 
        reads.sseqid_id = TOP_PERC.sseqid_id
    WHERE
        TOP_PERC.perc >= ${snakemake_params[topx]} 

-- Table no longer needed
DROP TABLE TOP_PERC;



-- Calculates the proportion of each reads entry in respect to the total number of good matches per read.
-- Write to parquet file, partitioned by sample id. 
-- By 'partitioning' we create a hive partiton readable directory 
.print "Creating hive partitioned parquet"
.print "Retain all hits that have bitscores at least TOPX precent as great as a reads best hits " 
COPY (
        SELECT
        reads.sample_id,
        reads.qseqid_key,
        reads.qseqid_id,
        reads.sseqid_id,
        reads.pident,
        reads.length,
        reads.mismatch,
        reads.gapopen,
        reads.qstart,
        reads.qend,
        reads.sstart,
        reads.send,
        reads.evalue,
        reads.bitscore,
        reads.perc,
        1 / NEW.n_hits AS proportion           -- proportion of entry to toal number of hits per read
    FROM
        inserts AS reads
    JOIN
        (   
        SELECT
            OLD.qseqid_key AS qseqid_key,
            count(OLD.qseqid_key) AS n_hits    -- Count number of hits per read
        FROM 
            inserts AS OLD
        GROUP BY
            OLD.qseqid_key
        ) AS NEW
    ON
        reads.qseqid_key = NEW.qseqid_key
    ORDER BY 
        reads.sseqid_id                        -- Order by protein key to enable faser querying of proteins 
    )
TO "${snakemake_output[0]}" (FORMAT PARQUET, OVERWRITE_OR_IGNORE true, PER_THREAD_OUTPUT true);
EOF
