# PATHOCOM pipeline

### PIPELINE
This Repo contains two of the **Snakemake** pipeline developed for the PATHCOM project's data analysis.
Both pipelines utilise DIAMOND, DuckDB, hive partition parquet files and some Python to process the 50TB collected by PATHOCOM.
Currently it is only made with PATHOCOM in mind but can be applied to any set of shotgun sequencing reads and will be generalized once it is completed.

### PATHOCOM
PATHOCOM is an ongoing ERC SYNERGY project of the Max-Plank Institute for Biology Tübingen Germany,  Laboratoire des interactions plantes - microbes - environnement, CNRS Toulouse, France and the New York University, NY, USA that aims to understand and predict pathogen communities on plants to better combat infectious diseases.
To achieve this goal it combines large-scale field observations of microbes in Arabidopsis thaliana with ultra-high-throughput experimental tests of host-dependent microbial interactions.
Two pipelines are used to transform the collected metagenome shotgun data into a usable database: **Samples preprocessing** and **Bacteria to uniref** 
(The project is still work in process, so these pipelines are neither finished nor currently in a clean and 100% up to date format) 

### Samples preprocessing:

Expects a folder with sample blasted against the clustered uniref database, each tsv file should at least contain qseqid, sseqid and bitscore.
For each sample it will create a parquet file, apply diamonds topx filtering (if not needed set topx = 0 in config file), calculates weighted target hits, and prepares per sample tables containing hit counts.
It then loads the NCBI mapping table, a custom table containing linage information and pepare SQL queries.
The result is a single DuckDB database, containing all per sample proteins abundance counts in an single table and prepared queries for faster data retrieval.
Prepared queries link single hits to gene, gene function, taxonomy etc. allowing near instant access to the metagenomic data found in every samples,  e.g select all proteins and their counts belonging to a single given GO term.
DuckDB databases and prepared queries can also be accessed by DuckDB's Python, C++ and Rust API's which allows for automated querying and plotting via external programs.  

### Bacteria to uniref: 

Used to compare assembles and raw shotgun sequencing reads from isolated bacterial samples to each other by aligning them against the clustered uniref50 reference database.
After DIAMOND blastx/blastp resulting hits are used to create protein abundance tables, use the NCBI mapping table to assign per hit GO terms and map the isolate results against the resulting tables from the 'sample preprocessing" table.
Extracts for the "read" and "assembled" data such as per sample Go term counts, per sample UNIREF counts and count value distribution.
Extracted data is then used to plot heatmaps and value distributions, see jupyter notebook.
Expects nucleotide fastq reads for raw read and amino acid fasta for assembled sequences. 

### Runtime
We ran both pipelines on the MPCDF Viper cluster, using 8 nodes the total runtime for 6000 file (~50TB of data) with samples preprocessing is less than 24 hours while the Bacteria to uniref pipeline needs around 18 hours to run the pipeline on 200 files containing assembled sequences and 200 files with raw reads.

## Installation
To run it you will have to install [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) (full or minimal) and [DuckDB CLI](https://duckdb.org/docs/installation/?version=stable&environment=cli&platform=linux&download_method=direct&architecture=x86_64). DuckDb requires no installation and can be run from the downloaded binary (you may have to add it to your path variable), databases are either in temporary and in memory or a written to a file given when launching DuckDB

**Snakemake (full):**

`conda create -c conda-forge -c bioconda -n snakemake snakemake`

**DuckDB:**

`curl https://install.duckdb.org | sh`

### Executing
You can run this pipeline via Snakemake on a single node/computer or run in on multiple nodes via the build in Slurm executor.

