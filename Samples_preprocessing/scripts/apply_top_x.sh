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
    qseqid varchar,
    sseqid varchar,
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
        qseqid,
        sseqid,
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
        read_parquet('${snakemake_input[sample_parquet]}/*.parquet', compression = 'zstd')
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
            qseqid,
            max(bitscore) AS top_score                 
        FROM
            inserts
        GROUP BY
            qseqid
        )                               
    SELECT
        inserts.qseqid AS qseqid,                      -- Query and Protein keys needed to enable later mapping of
        inserts.sseqid AS sseqid,                        -- perc to propper read/protein entry
        inserts.bitscore / TOP_SCORE.top_score * 100 AS perc   -- Percentage of hit's A bitscore in reads highest bitscore B
        
    FROM
        inserts
    JOIN
        TOP_SCORE
    USING 
        (qseqid)               -- JOINS read specific top_score to each entry, via read key
    ) 
;

-- Mapps the calculated percentage (perc) to each entry via query and subject id
-- Removes all rows with perc smaller than X  
.print 'Inserts top percentage, filter low perc' 
.print ${snakemake_params[topx]} 
CREATE OR REPLACE TEMP TABLE inserts AS 
    SELECT
        reads.qseqid,
        reads.sseqid,
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
        reads.qseqid = TOP_PERC.qseqid                
        AND 
        reads.sseqid = TOP_PERC.sseqid
    WHERE
        TOP_PERC.perc >=  ${snakemake_params[topx]} 
;

-- Table no longer needed
DROP TABLE TOP_PERC;



-- Calculates the proportion of each reads entry in respect to the total number of good matches per read.
-- Write to parquet file, partitioned by sample id. 
-- By 'partitioning' we create a hive partiton readable directory 
.print "Creating hive partitioned parquet"
.print "Retain all hits that have bitscores at least TOPX precent as great as a reads best hits " 
COPY (
        SELECT
        reads.qseqid,
        reads.qseqid_id,
        reads.sseqid,
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
            OLD.qseqid AS qseqid,
            count(OLD.qseqid) AS n_hits    -- Count number of hits per read
        FROM 
            inserts AS OLD
        GROUP BY
            OLD.qseqid
        ) AS NEW
    ON
        reads.qseqid = NEW.qseqid
    ORDER BY 
        reads.sseqid                        -- Order by protein key to enable faser querying of proteins 
    )
TO "${snakemake_output[sample_topx]}" (FORMAT PARQUET, OVERWRITE_OR_IGNORE true, PER_THREAD_OUTPUT true);
EOF
