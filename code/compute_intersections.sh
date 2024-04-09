#Calculate the intersecting set of samples
cd ../ccle_wgs/
file1="wgs_file_names.txt"  
file2="../data/ccle.vcf.sample.names" 
comm -12 <(sort "$file1") <(sort "$file2") > intersecting_wgs_snp6_samples.txt


#Create a file named directory_names.txt.
#It is easier to do this than it is to modify the other script
#Also put it in the ../data directory
cp intersecting_wgs_snp6_samples.txt ../data/directory_names.txt


#Loop through all of the files and calculate the intersection
cd ../ccle_wgs/
while IFS= read -r line; do
    echo "Processing line: $line"
	bcftools isec ${line}.ccle_all_called.vcf.gz ${line}.wgs.pass.vcf.gz -p $line -Oz &
done < intersecting_wgs_snp6_samples.txt



#Create places to store the data
cd ../ccle_wgs/
touch num_snp6_only.txt
touch num_wgs_only.txt
touch num_shared.txt


#Now loop through all of the output directories and record the number of variants for each file
for directory in ../ccle_wgs/ACH*/; do
    base_directory=$(pwd)
    cd "$directory"
	
	#Count the number of variants
	num_snp6_only=$(zcat 0000.vcf.gz | grep -v "#" | wc -l)
	num_wgs_only=$(zcat 0001.vcf.gz | grep -v "#" | wc -l)
	num_shared=$(zcat 0002.vcf.gz | grep -v "#" | wc -l)
	
	#Add this information to some new files
	echo ${num_snp6_only} >> ../num_snp6_only.txt
	echo ${num_wgs_only} >> ../num_wgs_only.txt
	echo ${num_shared} >> ../num_shared.txt
	
	#Leave the directory
    cd ${base_directory}
	
done