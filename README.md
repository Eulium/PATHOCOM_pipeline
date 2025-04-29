# PATHOCOM_pipeline
PATHOCOM is an ongoing ERC SYNERGY project of the Max-Plank Institute for Biology Tübingen Germany,  Laboratoire des interactions plantes - microbes - environnement, CNRS Tolouse, France and the New York University, NY, USA that aims to understand and predict pathogen communities on plants to better combat infectious diseases.
To achieve this goal it combines large-scale field observations of microbes in Arabidopsis thaliana with ultra-high-throughput experimental tests of host-dependent microbial interactions.
Two pipelines are used to transform the collected metagenome shotgun data into a usable database: **Samples preprocessing** and **Bacteria to uniref** 
(The project is still work in process, so these pipelines are neither finished nor currently in a clean and 100% up to date format) 

### Samples preprocessing:

Expects a folder with sample blasted aganist the clustered uniref database, each tsv file should at least contain qseqid, sseqid and bitscore.
For each sample it will create a parquet file apply diamonds topx filtering (if not needed set topx = 0 in config file), calculates weighted target hits, and prepares per sample tables containing hit counts.
The result is a single DuckDB database is created, containing all per sample proteins abundance counts in an single table.
It then loads the NCBI mapping table, a custom table containing linage information and pepare SQL queries.
The prepared queries link single hits to gene, gene funtion, taxonomy etc. 

Pepared sql queries can be accessed by DuckDB programm api's and allow for faster repeated queries e.g select all proteins and their counts belonging to a single given GO term. 

### Bacteria to uniref: 

Used to compare assembles and raw shotgun sequencing reads from isolated bacterial samples to each other by aligning them against the clustered uniref50 reference database.
After DIAMOND blastx/blastp resulting hits are used to create protein abundance tabels, use the ncbi mapping table to asigne per hit GO terms and map the isolate results against the resulting tables from the sample preprocessing table.
Extracts for the "read" and "assmbled" data such as per sample Go term counts, per sample UNIREF counts and count value distribution.
Extracted data is then plotted has comperative heatmaps etc.
Expects nucleotide fastq reads for raw read and amio acid fasta for assmbled sequences. 

### Runtime
We ran both pipelines on the MPCDF Viper cluster, using 8 nodes the total runtime for 6000 file (~50TB of data) with samples preprocessing is less than 24 hours while the Bacteria to uniref pipeline needs around 18 hours to run the pipeline on 200 files containing assmbled sequences and 200 files with raw reads.

## Instalation
To run it you will have to install [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) (full or minimal) and [DuckDB CLI](https://duckdb.org/docs/installation/?version=stable&environment=cli&platform=linux&download_method=direct&architecture=x86_64). DuckDb requires no installation and can be run from the downloded binary (you may have to add it to your path variable), databases are either in temporay and in memory or a writen to a file given when launching DuckDB

**Snakemake (full):**

`conda create -c conda-forge -c bioconda -n snakemake snakemake`

**DuckDB:**

`curl https://install.duckdb.org | sh`

### Executing
You can run this pipline via Snakemake on a single node, use the Snakemake"s build in Slurm executor or a custom slurm executor fort faster execution on multiple nodes. 


