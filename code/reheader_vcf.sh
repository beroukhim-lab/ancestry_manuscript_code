#The purpose of this script is to add a modified header to the ccle vcf files
#If this header isn't added, downstream steps will error because the vcf is not recognized correctly
#A 'small.header' file is included in this repository that contains the necessary information that must be added
#This script requires Bcftools and Tabix to be in PATH


#All vcf are stored here. The vcf files we want are those with the 'renamed' prefix
cd ../snp_array_data



for vcffile in renamed*;
do

#Make and add in the new header
head -4 $vcffile | tail -1 > fileheader.txt
cat smallheader.txt fileheader.txt > newheader.txt
bcftools reheader -h newheader.txt $vcffile -o $vcffile.newheader

#sort the file
cat $vcffile.newheader | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > $vcffile.sorted

#bgzip the file and generate the index
bgzip $vcffile.sorted
tabix -p vcf $vcffile.sorted.gz

done