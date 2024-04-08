## Introduction
This repository contains all code for performing the analyses and generating the figures described in the associated manuscript.

## Requirements
All jupyter notebooks were originally run in a GCP VM environment with 64 GB of memory and 8 CPU cores. These notebooks were adapted to run on a machine with 64 GB of memory and a Ryzen 7 series CPU (without parallel computing across cores). These notebooks are untested in environments with lower compute resources, but it is expected that most notebooks will fully execute with only 8 GB of memory and a single CPU core. A subset of notebooks require a large memory overhead and are not expected to execute with less than 32 GB of memory.

The shell scripts in this repository were originally run in a compute environment with 32 CPU cores and 4 GB of memory per core. These scripts assume such resources are available. These scripts can be modified to run on a single CPU core, but will take a long time to fully execute. 

## Input files
Download the input files and name them as indicated. All input files must be in the ./data directory of this repository.

'''

From the DepMap web portal: \
21Q4 Achilles Guide Map -> ./data/21q4_Achilles_guide_map.csv \
21Q4 Achilles Logfold Change -> .data/21q4_Achilles_logfold_change.csv \
21Q4 Achilles Replicate Map -> ./data/21q4_Achilles_replicate_map.csv \
21Q4 CRISPR Gene Effect -> ./data/21q4_crispr_gene_effect.csv \
22q1 Achilles Gene Effect -> ./data/22q1_Achilles_gene_effect.csv \
22q1 Achilles Guide Map -> ./data/22q1_Achilles_guide_map.csv \
22q1 Achilles Logfold Change -> ./data/22q1_Achilles_logfold_change.csv \

'''








code chunk for downloading data:
```

```
