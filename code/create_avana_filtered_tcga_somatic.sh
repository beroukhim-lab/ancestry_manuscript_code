#Loop through all of the somatic vcf files and filter them so that they only include variants in the avana hg38 bed file
#If you are using a big VM, consider running this in parallel (which requires significant compute resources)

#Loop through all of the files and filter
cd ../tcga_somatic/
for vcffile in *.vcf.gz;
do
bcftools filter -R ../data/Avana14_filtering.bed $vcffile -o avana.filtered.$vcffile -Oz
done
