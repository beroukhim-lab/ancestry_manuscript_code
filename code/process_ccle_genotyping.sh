#Unpack the tar.gz file
cd ../snp_array_data
gzip -d CCLE_SNP.Birdseed.Calls_2013-07-29.tar.gz
tar -xvf CCLE_SNP.Birdseed.Calls_2013-07-29.tar

#Move the files up one directory
cd ../snp_array_data/birdseed
for file in *.txt;
do
mv $file ..
done
cd ..


#Convert birdseed to vcf
cd ../snp_array_data

for file in ../snp_array_data/*.txt
do
fbase=$(basename $file | cut -f1 -d'.')

python3 ../code/BirdseedToVCF.py \
--birdseed "${file##*/}" \
--output_vcf $fbase.vcf \
--snp_annotation_file ../snp_array_data/GenomeWideSNP_6.na30.annot.hg19.csv.pickle \
--array_sample $fbase \
--vcf_sample $fbase \
--fasta ../snp_array_data/human_g1k_v37.fasta
done



#Rename the birdseed files
cd ../snp_array_data 
RScript file_renaming.R



#Update the header in the vcf files
cd ../snp_array_data 
bash reheader_vcf.sh



#Merge the vcf files and then index
cd ../snp_array_data
bcftools merge *.vcf.sorted.gz -o ccle_snp6.vcf
bgzip ccle_snp6.vcf
tabix -p vcf ccle_snp6.vcf.gz



#Split the vcf into individual chromosomes (necessary for imputation)
cd ../snp_array_data
for i in {1..22};
do
tabix -h ccle_snp6.vcf.gz ${i} > ccle.all.called.chr${i}.vcf
bgzip ccle.all.called.chr${i}.vcf
done

#Also do the same for chrX
tabix -h ccle_snp6.vcf.gz X > ccle.all.called.chrX.vcf
bgzip ccle.all.called.chrX.vcf
