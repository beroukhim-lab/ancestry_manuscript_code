#Now extract those SNPs with bcftools
cd ../data
bcftools filter --include 'ID=@top_snp_for_extraction.txt' ccle_all_called.vcf.gz > extracted.snps