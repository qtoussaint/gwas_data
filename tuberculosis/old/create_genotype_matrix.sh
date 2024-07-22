# First get the merged VCF from SNP calling. Then run the following command to extract the VCF: 

fg
nice bcftools query -f '[%GT,]\n' merged_files.vcf.gz | sed -e 's/\.\/\./0/g' -e 's/1\/1/1/g' -e 's/,$//' > genotype_matrix.csv
