#Filter the TCGA germline data so that it only includes variants that map to guides 
cd ../tcga_germline/
bcftools filter -R ../data/Avana14_filtering_hg19.bed PCA.r1.TCGAbarcode.merge.tnSwapCorrected.10389.vcf.gz -o tcga_germline_variants_in_avana_guides.vcf.gz -Oz
