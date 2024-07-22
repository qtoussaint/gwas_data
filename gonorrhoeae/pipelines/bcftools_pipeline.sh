# Pipeline to concatonate individual VCFs (previously generated using snippy) using the bcftools package

bcftools merge /nfs/research/jlees/jacqueline/gwas_datasets/gonorrhoeae/genotype/variant_calls/results*/snps.vcf.gz -o merged_files.vcf.gz


