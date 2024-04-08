#Get the directory of the script, which should be in the 'data' subdirectory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
cd "$SCRIPT_DIR"/..
BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"



#########
#Somatic#
#########

#Loop through all of the somatic vcf files and filter them so that they only include variants in the avana hg38 bed file
#Note: This takes a lot of compute to run because it's going to try to run everything in parallel.
cd ${BASE_DIR}/tcga_somatic
for vcffile in *.vcf.gz;
do
bcftools filter -R ${BASE_DIR}/data/Avana14_filtering.bed $vcffile -o avana.filtered.$vcffile -O z &
done





##########
#Germline#
##########

#Filter the TCGA germline data
cd ${BASE_DIR}/tcga_germline
bcftools filter -R Avana14_filtering_hg19.bed PCA.r1.TCGAbarcode.merge.tnSwapCorrected.10389.vcf.gz -o tcga_germline_variants_in_avana_guides.vcf.gz -O z

#Extract the sample names from the germline vcf file
cd ${BASE_DIR}/tcga_germline
bcftools query -l PCA.r1.TCGAbarcode.merge.tnSwapCorrected.10389.vcf.gz > tcga_germline_sample_names.txt

