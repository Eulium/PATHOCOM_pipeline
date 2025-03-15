# PATHOCOM_pipeline
PATHOCOM is an ongoing ERC SYNERGY project of the Max-Plank Institute for Biology Tübingen Germany,  Laboratoire des interactions plantes - microbes - environnement, CNRS Tolouse, France and the New York University, NY, USA that aims to understand and predict pathogen communities on plants to better combat infectious diseases.
To achieve this goal it combines large-scale field observations of microbes in Arabidopsis thaliana with ultra-high-throughput experimental tests of host-dependent microbial interactions.
Two pipelines are used to transform the collected metagenome shotgun data into a usable database: **Samples preprocessing** and **Bacteria to uniref** 
(The project is still work in process, so these pipelines are neither finished nor currently in a clean and 100% up to date format) 
### Samples preprocessing:

Expects a folder with sample blasted aganist the clustered uniref database, each tsv file should at least contain qseqid, sseqid and bitscore. For each sample it will create a parquet file apply diamonds topx filtering (if not needed set topx = 0 in config file) and calculates per sample protein aboundance tables with counts.
At the end a single DuckDB database is created, containing all per sample proteins abundance counts in an single table.
The file querys.sql can be used to load the NCBI mapping table and pepare SQL queries, by GO term and  UniProtKB_AC.
Pepared sql queries can be accessed by DuckDB programm api's and allow for faster repeated queries. 

### Bacteria to uniref: 

Used to mapp bacterial isolates against sample 


## Instalation
To run it you will have to install [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) (full or minimal) and [DuckDB CLI](https://duckdb.org/docs/installation/?version=stable&environment=cli&platform=linux&download_method=direct&architecture=x86_64). DuckDb requires no installation and can be run from the downloded binary (you may have to add it to your path variable), databases are either in temporay and in memory or a writen to a file given when launching DuckDB

**Snakemake (full):**

`conda create -c conda-forge -c bioconda -n snakemake snakemake`

**DuckDB:**

`curl https://install.duckdb.org | sh`


