#Extract the sample names from the vcf file
cd ../tcga_germline
bcftools query -l PCA.r1.TCGAbarcode.merge.tnSwapCorrected.10389.vcf.gz > tcga_germline_sample_names.txt
