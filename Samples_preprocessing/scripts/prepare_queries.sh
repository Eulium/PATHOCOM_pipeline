
# Loads ncbi mapping table as well as prepares queries for fast search by GO term and taxon id.
# Prepared quries can be effieciently mass executed via DuckDB API's in e.g. Python.
duckdb ${snakemake_input[DB]}  << EOF
        
-- LOAD ncbi_mapping table into Database, this allows us to connect our sample result to taxids, GO terms and other information.
CREATE TABLE ncbi_mapping AS
    SELECT 
        *
    FROM
        read_csv("${snakemake_input[NCBI_mapping]}",
            header = false,
            names = [UniProtKB_AC, UniProtKB_ID, GeneID, RefSeq, GI, PDB, GO, UniRef100, UniRef90, UniRef50, UniParc, PIR, NCBI_taxon, MIM, UniGene, PubMed, EMBL, EMBL_CDS, Ensembl, Ensembl_TRS, Ensembl_PRO, Additional_PubMed],
            delim = '\t',
            compression = 'gzip')
;
-- Split GO term entries from single string into list for faster search later on
ALTER TABLE ncbi_mapping ALTER GO SET DATA TYPE VARCHAR[] using string_split(GO, ';');

-- Load a table containing taxon linages as list, by grouping a Kraken based metagneome classifictaion results table into records based on tax id and a counter i.
-- The Kraken table was previously manualy compiled calculated and had the counters extened.
-- This allows search by linge information for a given taxon id.
CREATE TABLE NCBI_Taxon_lineages AS 
    SELECT  
        taxid AS NCBI_taxon,
        list(lineage_taxid) AS Lineage
    FROM 
        (
        SELECT 
            *
        FROM
            read_parquet(${snakemake_input[Linage_file)]}",
        ORDER BY taxid, i
        )
    GROUP BY
        taxid
    ; 

-- Create repared qyery for fast search by GO term, prepared queries can be mass executed efficiently via DuckDB API's in e.g. Python

-- Prepare query for search by a single GO term and a single taxon id.
-- This query will return all protein abundances that are grouped with the given GO term and the given NCBI_taxon_id is part of the linage table.
PREPARE Query_BY_GO_Taxa AS 
    WITH SUBSET AS
        (
            SELECT
                DISTINCT ON(UniRef90) UniRef90
            FROM
                ncbi_mapping
            JOIN
                (
                SELECT
                    UniProtKB_AC
                FROM
                    ncbi_mapping
                JOIN
                    NCBI_Taxon_lineages
                USING	
                    (NCBI_taxon)
                WHERE
                    contains(ncbi_mapping.GO, $GO_term)
                    AND
                    list_contains(NCBI_Taxon_lineages.Lineage, $NCBI_taxon_id)
                ) AS A
            ON
                    ncbi_mapping.UniProtKB_AC = A.UniProtKB_AC
        )
    SELECT
            *
    FROM
            Protein_abundance
    JOIN
            SUBSET
    ON
            Protein_abundance.sseqid = SUBSET.UniRef90
;

-- Prepare query for protein abundance search by a single GO term.
PREPARE grep_GO_hits AS
    WITH SUBSET AS
    (
        SELECT
            DISTINCT ON(UniRef90) UniRef90
        FROM
            ncbi_mapping
        JOIN
            (
            SELECT
                UniProtKB_AC
            FROM
                ncbi_mapping
            WHERE
                contains(GO, $GO_term)
            )
        USING
            (UniProtKB_AC)
        )
    SELECT
        *
    FROM
        Protein_abundance
    JOIN
        SUBSET
    ON
        Protein_abundance.sseqid = SUBSET.UniRef90 

EOF