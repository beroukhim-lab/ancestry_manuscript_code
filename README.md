## Introduction
This repository contains all code for performing the analyses and generating the figures described in the associated manuscript.

## Requirements
All jupyter notebooks were originally run in a GCP VM environment with 64 GB of memory and 8 CPU cores. These notebooks were adapted to run on a machine with 64 GB of memory and a Ryzen 7 series CPU. These notebooks are untested in environments with lower compute resources, but it is expected that most notebooks will fully execute with only 8 GB of memory and a single CPU. A subset of notebooks require a large memory overhead and are not expected to execute with less than 32 GB of memory.

The shell scripts in this repository were originally run in a compute environment with 32 CPU cores and 128 GB of memory. These scripts have been modified to run on a single CPU core (at a significant time cost), but can be modified to run in parallel.


## Simple Startup
For simplicity, critical input and intermediate files are provided in this github repository or are hosted on figshare. Many of these files have significant compute requirements to generate. Full descriptions on how to generate these inputs are below, but using the pre-generated inputs is strongly recommended.

To get started, clone this github repository and download inputs from figshare as follows:
```
#From DepMap figshare
wget -O ./data/21q4_Achilles_guide_map.csv https://figshare.com/ndownloader/files/31315819
wget -O .data/21q4_Achilles_logfold_change.csv https://figshare.com/ndownloader/files/31315903 
wget -O ./data/21q4_Achilles_replicate_map.csv https://figshare.com/ndownloader/files/31315876 
wget -O ./data/21q4_crispr_gene_effect.csv https://figshare.com/ndownloader/files/31315996 
wget -O ./data/22q1_Achilles_gene_effect.csv https://figshare.com/ndownloader/files/34008383 
wget -O ./data/22q1_Achilles_guide_map.csv https://figshare.com/ndownloader/files/34008362 
wget -O ./data/22q1_Achilles_logfold_change.csv https://figshare.com/ndownloader/files/34008443 
wget -O ./data/22q1_Achilles_replicate_map.csv https://figshare.com/ndownloader/files/34008398 
wget -O ./data/22q1_CCLE_gene_cn.csv https://figshare.com/ndownloader/files/34008428 
wget -O ./data/22q1_crispr_gene_effect.csv https://figshare.com/ndownloader/files/34008491 
wget -O ./data/22q1_expression.csv https://figshare.com/ndownloader/files/34008404 
wget -O ./data/22q1_sample_info.csv https://figshare.com/ndownloader/files/34008503 
wget -O ./data/22q2_crispr_gene_effect.csv https://figshare.com/ndownloader/files/34990036 
wget -O ./data/OmicsGuideMutationsBinaryKY.csv https://plus.figshare.com/ndownloader/files/43347465 
wget -O ./data/OmicsGuideMutationsBinaryHumagne.csv https://plus.figshare.com/ndownloader/files/43347432 
wget -O ./data/OmicsGuideMutationsBinaryAvana.csv https://plus.figshare.com/ndownloader/files/43347378 
wget -O ./data/AvanaGuideMap.csv https://plus.figshare.com/ndownloader/files/43346391 
wget -O ./data/HumagneGuideMap.csv https://plus.figshare.com/ndownloader/files/43346748 
wget -O ./data/KYGuideMap.csv https://plus.figshare.com/ndownloader/files/43346769 
wget -O ./data/common_essentials.csv https://plus.figshare.com/ndownloader/files/43346361 

#From figshare for this paper
wget -O ./data/final_avana.txt https://figshare.com/ndownloader/files/45747036 
wget -O ./data/final_moffat.txt https://figshare.com/ndownloader/files/45747045 
wget -O ./data/genetic_map_hg38_withX.txt https://figshare.com/ndownloader/files/45746991 
wget -O ./data/final_custom_library.txt https://figshare.com/ndownloader/files/45747039 
wget -O ./data/final_dolcetto.txt https://figshare.com/ndownloader/files/45747054 
wget -O ./data/gene.block.matrix.txt https://figshare.com/ndownloader/files/45747027 
wget -O ./data/avana.filtered.ccle.variant.calls.vcf.gz https://figshare.com/ndownloader/files/45747048 
wget -O ./data/final_sanger.txt https://figshare.com/ndownloader/files/45747003 
wget -O ./data/final_gecko.txt https://figshare.com/ndownloader/files/45747006 
wget -O ./data/final_minlibcas.txt https://figshare.com/ndownloader/files/45747024 
wget -O ./data/final_calabrese.txt https://figshare.com/ndownloader/files/45746943 
```



## Generating input files
Download or create the input files and name them as indicated. All input files must be in the ./data directory of this repository unless otherwise indicated.
Some previously generated input files are restricted access and require the user to apply for download permissions from the respective study. The following code will generate all inputs necessary for running the jupyter notebooks. Many processing steps are compute-intensive, and when possible pre-computed input files are provided in the ./data directory.

#### From the DepMap web portal (https://depmap.org/portal/download/all/): 
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
CCLE_SNP.Birdseed.Calls_2013-07-29.tar.gz -> ./snp_array_data
OmicsGuideMutationsBinaryKY -> ./data/OmicsGuideMutationsBinaryKY.csv
OmicsGuideMutationsBinaryHumagne -> ./data/OmicsGuideMutationsBinaryHumagne.csv
OmicsGuideMutationsBinaryAvana -> ./data/OmicsGuideMutationsBinaryAvana.csv
AvanaGuideMap -> ./data/AvanaGuideMap.csv
HumagneGuideMap -> ./data/HumagneGUideMap.csv
KYGUideMap -> ./data/KYGuideMap.csv
AchillesCommonEssentialControls -> ./data/common_essentials.csv
```

#### CCLE WGS (Note: Controlled access)
```
#Get access to the CCLE WGS controlled access data
#Download the files listed in ./ccle_wgs/wgs_file_names.txt
#Store all files in ./ccle_wgs/

#Create num_shared.txt
#Create num_snp6_only.txt
#Create num_wgs_only.txt
#Requires Bcftools (in PATH)
bash compute_intersections.sh

```

#### From UCSC genome browser:
```
#Download hg38 gtf
cd ./data
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz
gzip -d hg38.refGene.gtf.gz
```

#### TCGA somatic data from gdc portal (portal.gdc.cancer.gov)
```
#Download the somatic mutation calls listed in ./tcga_somatic/vcf_sample_sheet.tsv
#Store all files in ./tcga_somatic/

#Filter the files to only include variants that map to avana guides
cd ./code
bash create_avana_filtered_tcga_somatic.sh
```

#### TCGA germline variant calls (note: restricted access)
```
#TCGA germline variant calls are included as part of the following manuscript:
#Pathogenic germline variants in 10,389 adult cancers
#PMID: 29625052
#Download the variant call file (name = PCA.r1.TCGAbarcode.merge.tnSwapCorrected.10389.vcf.gz). Store in ./tcga_germline/

#Extract the sample names
#Requires Bcftools (in PATH)
cd ./code
bash extract_tcga_germline_sample_names.sh

#Filter the germline variant calls to include only variants in avana guides
#Requires Bcftools (in PATH)
cd ./code
bash extract_tcga_germline_variants_in_avana_guides.sh
```


#### Processing CCLE SNP6 genotyping
```
#Requires Python-3.6 (must be in PATH)
#Requires R-4.0 (must be in PATH)
#Requires Bcftools (must be in PATH)
#Requires Tabix (must be in PATH)

#Download the SNP6 annotation file
cd ./snp_array_data
wget https://software.broadinstitute.org/cancer/cga/sites/default/files/data/tools/contest/GenomeWideSNP_6.na30.annot.hg19.csv.pickle.gz

#Download the hg19 fasta
cd ./snp_array_data
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
gzip -d human_g1k_v37.fasta.gz
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.fai


#Process the vcf files. This script will take a long time to run, but can easily be modified to process vcf files in parallel
cd ./snp_array_data
bash process_ccle_genotyping.sh

############
#Imputation#
############
#Step 1) Perform genotype phasing and imputation using the Michigan Imputation Server
#Step 2) Download and unzip the phased/imputed vcf files. The zipped file is password protected, follow unpacking instructions.
#Step 3) Add the 'pe_' prefix to the phased/imputed vcf files.


#Cat the phased/imputed vcfs, pass filter, then zip and index
cd ./snp_array_data
bcftools concat -o ccle_snp6_phased_imputed.vcf pe_*
bcftools view -f PASS ccle_snp6_phased_imputed.vcf > ccle_all_called.vcf
bgzip ccle_all_called.vcf
tabix -p vcf ccle_all_called.vcf.gz
cp ccle_all_called.vcf.gz ../data
cp ccle_all_called.vcf.gz.tbi ../data


#Extract the sample headers, which is a useful input for many scripts
cd ./snp_array_data
bcftools query -l ccle_all_called.vcf.gz > ../data/ccle.vcf.sample.names.txt

```

#### Guide library filtering
```
#Filter the CCLE variant calls (ccle_all_called.vcf.gz) to retain only SNPS in Avana guides
#This script will create 'snps.in.all.avana.guides.vcf.gz'
#Requires Bcftools (in PATH)
cd ./code
bash create_snps_in_all_avana_guides_vcf.sh
```

#### Running RFMix-v2
```
#Requires RFMixv2 (https://github.com/slowkoni/rfmix)
#Requres Samtools-v1.9 (in PATH)
#Requires Bcftools (in PATH)
#Requires R-4.0 (in PATH)


#Download the RFMixv2 reference panel
cd ./data
for i in {1..22};
do
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr${i}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr${i}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz.tbi
done

#Download the genetic map (https://alkesgroup.broadinstitute.org/Eagle/#x1-250005.1.2)
#Download genetic_map_hg38_withX.txt (retain this name) and put in ./data
#Format the genetic map for our analysis
cd ./code
Rscript format_genetic_map.R


#Run RFMix
cd ./code
bash run_rfmix.sh
```

#### From R scripts:
```
#Create depmap_cell_lineage.csv
#Requires R-4.0 (or higher)
cd ./code
Rscript create_depmap_cell_lineage_file.R


#Create gene.block.matrix.txt
#Requires R-4.0 (or higher)
cd ./code
Rscript create_gene_block_matrix.R


#Create ancestry_top_snp_df.txt and merged.pvals.txt
#Requires R-4.0 (or higher)
#Requires PLINK 2.0 (in PATH)
cd ./code
Rscript create_ancestry_top_snp_df.txt


#Create merged_frequency_dataset.txt and ccle.ancestry.snps.vcf.gz
#Requires R-4.0 (or higher)
#Requires Bcftools (in PATH)
cd ./code
Rscript create_merged_frequency_dataset.R


#Create collapsed.ancestry.information.txt
#Requires R-4.0 (or higher)
cd ./code
Rscript create_collapsed_ancestry_information.R


#Create snv_position_single_guide_finaldf.txt
#Requires R-4.0 (or higher)
cd ./code
Rscript create_snv_position_single_guide_finaldf.R


#Create chronos_22q1_ancestry_associated_dependency_pvals.txt
#Create chronos_22q2_ancestry_associated_dependency_pvals.txt
#Requires R-4.0 (or higher)
cd ./code
Rscript create_chronos_22q1_22q2_ancestry_associated_dependency_pvals.R


#Create ccle_snp6_ancestry_calls.txt
#Requires R-4.0 (or higher)
cd ./code
Rscript create_ccle_snp6_ancestry_calls.R


#Create affected_guides_per_gnomad_sample.txt
#Require R-4.0 (or higher)
cd ./code
Rscript create_affected_guides_per_gnomad_sample.R


#Create top_snp_for_extraction.txt
#Requires R-4.0 (or higher)
cd ./code
Rscript create_top_snp_for_extraction.R


#Create top.snp.fdr.txt
#Create merged.pvals.txt
#Requires R-4.0 (or higher)
cd ./code
Rscript create_top_snp_fdr.R
```

#### From Doench et. al. 2016
```
Mismatch Tab from Supplemental Table 19 -> ./data/Doench_Data.txt
```


#### From COSMIC website
```
Download and unpack Cosmic_CancerGeneCensus_Tsv_v97_GRCh38.tar -> ./data/cosmic_genes.csv
```

#### From gnomAD. Note: Total file size is >5 TB
```
#Download the hgdp+1kg subset info file
gsutil cp https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_1kg_subset_sample_meta.tsv.bgz ./data

#Download the individualized hgdb+1kg vcf and index files
for i in {1..22};
do
gsutil cp https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr${i}.vcf.bgz ./data
gsutil cp https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr${i}.vcf.bgz.tbi ./data
done

#Download the aggregated vcf and index files
for i in {1..22};
do
gsutil cp https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr${i}.vcf.bgz ./data
gsutil cp https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr${i}.vcf.bgz.tbi ./data
done
```

#### From jupyter notebooks:
```
figure_1d.ipynb -> ./data/lm_ancestry_associated_dependency_pvals.txt
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
figure_3e_supplementalfigure_8.ipynb -> snv_position_single_guide_finaldf.txt
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

