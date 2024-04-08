## Introduction
This repository contains all code for performing the analyses and generating the figures described in the associated manuscript.

## Requirements
All jupyter notebooks were originally run in a GCP VM environment with 64 GB of memory and 8 CPU cores. These notebooks were adapted to run on a machine with 64 GB of memory and a Ryzen 7 series CPU. These notebooks are untested in environments with lower compute resources, but it is expected that most notebooks will fully execute with only 8 GB of memory and a single CPU. A subset of notebooks require a large memory overhead and are not expected to execute with less than 32 GB of memory.

The shell scripts in this repository were originally run in a compute environment with 32 CPU cores and 128 GB of memory. These scripts assume such resources are available. These scripts can be modified to run on a single CPU core, but will take a long time to fully execute. 

## Generating input files
Download or create the input files and name them as indicated. All input files must be in the ./data directory of this repository.
Some previously generated input files are restricted access. You must apply for access to get download permissions.

#### From the DepMap web portal: 
```
21Q4 Achilles Guide Map -> ./data/21q4_Achilles_guide_map.csv 
21Q4 Achilles Logfold Change -> .data/21q4_Achilles_logfold_change.csv 
21Q4 Achilles Replicate Map -> ./data/21q4_Achilles_replicate_map.csv 
21Q4 CRISPR Gene Effect -> ./data/21q4_crispr_gene_effect.csv 
22q1 Achilles Gene Effect -> ./data/22q1_Achilles_gene_effect.csv 
22q1 Achilles Guide Map -> ./data/22q1_Achilles_guide_map.csv 
22q1 Achilles Logfold Change -> ./data/22q1_Achilles_logfold_change.csv
22q1 Achilles Replicate Map -> ./data/22q1_Achilles_replicate_map.csv
22q1 CCLE Gene CN -> ./data/22q1_CCLE_gene_cn.csv
22q1 CRISPR Gene Effect -> ./data/22q1_crispr_gene_effect.csv
22q1 Expression -> ./data/22q1_expression.csv
22q1 Sample Info -> ./data/22q1_sample_info.csv
23q4 Omics Signatures -> ./data/23q4_omics_signatures.csv
22q2 CRISPR Gene Effect -> ./data/22q2_crispr_gene_effect.csv
23q4 OmicsSignatures -> ./data/23q4_omics_signatures.csv
```

#### From jupyter notebooks:
```
figure_1d.ipynb -> ./data/lm_ancestry_associated_dependency_pvals.txt
```


#### From shell scripts:
```
```

#### From UCSC genome browser:
```
#Download hg38 gtf
cd ./data
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz
gzip -d hg38.refGene.gtf.gz

```


## Notebook outputs
This section describes the expected output files from running each notebook. All output files are written to ./output unless otherwise noted.

```
figure_1b.ipynb -> fibure_1b_data.txt
figure_1b.ipynb -> figure_1b.pdf
figure_1b.ipynb -> collapsed_snp6_ancestry_calls.txt
figure_1c.ipynb -> figure_1c_data.txt
figure_1c.ipynb -> figure_1c.pdf
figure_1d.ipynb -> ./data/lm_ancestry_Associated_dependency_pvals.txt
figure_1d.ipynb -> figure_1d.pdf
figure_1d.ipynb -> figure_1d_figure_df.txt
figure_1e.ipynb -> figure_1e_data.txt
figure_1e.ipynb -> figure_1e.pdf
figure_2a_2b.ipynb -> figure_2a.pdf
figure_2a_2b.ipynb -> figure_2b.pdf
figure_2c.ipynb -> figure_2c_distance_to_tss_df.txt
figure_2c.ipynb -> figure_2c.pdf
figure_2d.ipynb -> figure_2d.pdf
figure_2d.ipynb -> figure_2d_eqtl_summary_df.txt
figure_2e.ipynb -> figure_2e.pdf
figure_3a.ipynb -> figure_3a_compiled_difference.txt
figure_3a.ipynb -> figure_3a.pdf
figure_3b.ipynb -> avana_filtering_bed_file.bed
figure_3b.ipynb -> figure_3b.pdf
figure_3b.ipynb -> figure_3b_distribution_table.txt
figure_3c.ipynb -> figure_3c.pdf
figure_3c.ipynb -> figure_3c_avana_affected_rate.txt
figure_3d.ipynb -> figure_3d.pdf
figure_3d.ipynb -> figure_3d.pdf
figure_3d.ipynb -> figure_3d_germline_somatic.txt
figure_3e_supplementalfigure_8.ipynb -> figure_3e_snp_in_guide.txt
figure_3e_supplementalfigure_8.ipynb -> figure_3e.pdf
figure_3e_supplementalfigure_8.ipynb -> figure3_plotting_df.txt
figure_3e_supplementalfigure_8.ipynb -> supplemental_figure_8.pdf
figure_3e_supplementalfigure_8.ipynb -> supplemental_figure_8_ours_doench_merged.txt
figure_4a_4b.ipynb -> figure_4a_affected_guides_df.txt
figure_4a_4b.ipynb -> figure_4a_top.pdf
figure_4a_4b.ipynb -> figure_4a_bottom.pdf
figure_4a_4b.ipynb -> figure_4a_complete_counts.txt
figure_4a_4b.ipynb -> figure_4b.pdf
figure_4a_4b.ipynb -> figure_4b_complete_counts.txt
figure_4c_4d.ipynb -> figure_4c.pdf
figure_4c_4d.ipynb -> figure_4c_affected_genes_per_person.txt
figure_4c_4d.ipynb -> figure_4d.pdf
figure_4c_4d.ipynb -> figure_4d_cosmic_matrix.txt
figure_4e.ipynb -> figure_4e.pdf
figure_4e.ipynb -> figure_4e_scatterplot_df.txt
supplemental_figure_1.ipynb -> supplemental_fig1_achilles_only_ancestry_pvals.txt
supplemental_figure_1.ipynb -> supplemental_fig1.pdf
supplemental_figure_1.ipynb -> supplemental_figure1_data_table.txt
supplemental_figure_2.ipynb -> supplemental_figure_2.pdf
supplemental_figure_2.ipynb -> supplemental_figure_2_merged_cn_snp.txt
supplemental_figure_3.ipynb -> supplemental_figure_3_left.pdf
supplemental_figure_3.ipynb -> supplemental_figure_3_right.pdf
supplemental_figure_3.ipynb -> supplemental_figure_3_df.txt
supplemental_figure_4.ipynb -> supplemental_figure_4.pdf
supplemental_figure_4.ipynb -> supplemental_figure_4_differential_df.txt
supplemental_figure_5_6.ipynb -> supplemental_figure_5_6_a.pdf
supplemental_figure_5_6.ipynb -> supplemental_figure_5_6_b.pdf
supplemental_figure_5_6.ipynb -> supplemental_figure_5_6_c.pdf
supplemental_figure_5_6.ipynb -> supplemental_figure_5_6_df.txt
supplemental_figure_7.ipynb -> supplemental_figure_7_snp_in_guide_df.txt
supplemental_figure_7.ipynb -> supplemental_figure_7_df.txt
supplemental_figure_7.ipynb -> supplemental_figure_7.pdf
supplemental_figure_11.ipynb -> supplemental_figure_11.pdf
supplemental_figure_11.ipynb -> supplemental_figure_11.df.txt

```

