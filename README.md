**EVEREST (pipEline for Viral assEmbly and chaRactEriSaTion)**
**This version, meta, is adapted for pair-end and single-end reads.**

EVEREST is a snakemake pipeline for virus discovery, its main purpose is to characterise phage genomes isolates but can be also used for all the virome.

**Running EVEREST:**

1. Clone EVEREST repository.
```
(base)$ git clone --recursive https://github.com/agudeloromero/EVEREST.git
```

2. Create everest environment
Creating conda environment for Snakemake and EVEREST from the file everest.yml, provided in this repository.
```
(base)$ conda env -n everest create -f everest.yml
(base)$ conda activate everest
(everest)$
```

3. Databases for EVEREST are available here: [https://zenodo.org/records/8404860](https://zenodo.org/records/8404860).

**How to run EVEREST:**

* Edit conf/config.yml file to point the input, output and databases directories.
* The databases include the human genome as a reference to get rid off those reads. However, this link can be changed by an alternative genome. Example of how to download other genomes of interest can be seen here (https://github.com/agudeloromero/Reference_Genomes).
```
$ conda activate everest
(everest)$ snakemake --use-conda -j 2 --keep-going
```
Input directory should contain the .fastq.gz files to analyse ( -j option have to be set depending on number of available cores).
Files are expected to be named as "name_R1" and "name_R2" plus extension. In case you need to rename then, please see this example (https://github.com/agudeloromero/rename_fastq_files).
