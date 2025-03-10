#!/bin/bash
# Transform a given tsv.zstd file into a parquet file
# Mapps each read id to a unique int number 

duckdb "${snakemake_output[1]}" << EOF
    SET threads = ${snakemake[threads]};     
    CREATE sequence qseqid_key minvalue 0 START WITH 0;

    CREATE temp TABLE qseqid_mapping AS 
    SELECT 
            DISTINCT ON(qseqid_id)
            qseqid_id,
            nextval('qseqid_key') AS qseqid_key
    FROM 
            read_csv(
            '${snakemake_input[0]}',
            delim = '\t', 
            header = false,
            compression = 'zstd',
            columns = {
            'qseqid_id':'varchar', 'sseqid_id':'varchar', 'pident':'real',
            'length':'SMALLINT', 'mismatch':'SMALLINT', 'gapopen':'SMALLINT', 'qstart':'SMALLINT',
            'qend':'SMALLINT', 'sstart':'INT', 'send':'INT', 'evalue':'real', 'bitscore':'real'}
            );

    COPY 
    (
    SELECT
        qseqid_id,
        qseqid_mapping.qseqid_key AS qseqid_key,
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
        IF(length = -1, false, true) AS mapped
    FROM            -- Read list of FILES, add filenames 
        read_csv(
        '${snakemake_input[0]}',
        delim = '\t', 
        header = false,
        compression = 'zstd',
        columns = {
        'qseqid_id':'varchar', 'sseqid_id':'varchar', 'pident':'real',
        'length':'SMALLINT', 'mismatch':'SMALLINT', 'gapopen':'SMALLINT', 'qstart':'SMALLINT',
        'qend':'SMALLINT', 'sstart':'INT', 'send':'INT', 'evalue':'real', 'bitscore':'real'}
        )
        JOIN 
                qseqid_mapping
        USING   (qseqid_id)
    )
    TO '${snakemake_output[0]}' (FORMAT parquet, COMPRESSION 'zstd', COMPRESSION_LEVEL 10, PER_THREAD_OUTPUT true);
EOF
