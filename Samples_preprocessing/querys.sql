.timer on
.echo on

-- A collection of SQL prepared querys for showcasing different DB query possibilities as well as DuckDB querying speed.
-- Prepared statments can also be called from DuckDB's python connection



.print 'NCBI_Taxon_lineages'

CREATE TABLE NCBI_Taxon_lineages
AS 
SELECT  
	taxid AS NCBI_taxon,
	list(lineage_taxid) AS Lineage
FROM 	(
	SELECT 
		*
	FROM
		read_parquet('/ptmp/pboppert/Pathocom/diamond-postprocess-playground/taxonomy/lineage.pqt')
	ORDER BY taxid, i
	)
GROUP BY
	taxid
; 


.print 'prep'
PREPARE Query_BY_GO_Taxa AS 
COPY (
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
                                 contains(ncbi_mapping.GO, $1)
                         AND
                                 list_contains(NCBI_Taxon_lineages.Lineage, $2)
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
         Protein_abundance.sseqid_id = SUBSET.UniRef90
) 
TO 'GO:Taxa_0006826_135621' (FORMAT PARQUET);  



PREPARE grep_GO_hits AS
WITH SUBSET AS
(
SELECT
        UniRef90_key
FROM
        ncbi_mapping_new
WHERE 
        ncbi_mapping_new.UniProtKB_AC_key IN 
        (
        SELECT
                UniProtKB_AC_key
        FROM
                ncbi_mapping_new
        WHERE
                list_contains(GO, ?)
        )
)
SELECT
        *
FROM
        Protein_abundance_mapped
JOIN
        SUBSET
ON
        Protein_abundance_mapped.UniRef90_key = SUBSET.UniRef90_key
;


.print 'EXECUTE Query GO:0006826,135621'
Execute Query_BY_GO_Taxa('GO:0006826',135621); 



.print 'Unmapped'

COPY (
WITH SUBSET AS
(
        SELECT
                 distinct on(UniRef90) UniRef90
        FROM
                 ncbi_mapping
        JOIN
                 (
                         SELECT
                                 UniProtKB_AC
                         FROM
                                 ncbi_mapping
                         WHERE
                                 contains(GO, 'GO:0006826')
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
         Protein_abundance.sseqid_id = SUBSET.UniRef90 

) TO 'GO:0006826_unmapped' (Format parquet); 


COPY (
 WITH subset AS
  (
        SELECT
                UniProtKB_AC_key,
                UniRef90_key
        FROM
                ncbi_mapping_new
        WHERE
                UniRef50 = 'UniRef50_Q887C9'
 )
        SELECT
                *

        FROM
                Protein_abundance
        JOIN
                subset
        USING 
                (UniRef90_key)
)
TO './Pseudomonas_effector_avrE_hits.csv' (FORMAT CSV, header true, delim '\t')
;
