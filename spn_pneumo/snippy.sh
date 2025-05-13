#!/usr/bin/env bash

#SBATCH --time=48:00:00             # walltime
#SBATCH --ntasks=1                  # number of tasks
#SBATCH --cpus-per-task=48          # number of CPUs Per Task
#SBATCH --nodes=1                   # number of nodes
#SBATCH --mem=150GB                  # memory per node
#SBATCH -J "snippy"             # job name
#SBATCH -o "bcftools_merge.out"         # job output file
#SBATCH -e "bcftools_merge.err"         # job error file

#################################################################################

outdir=/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype

cd $outdir

# Make a list of all assemblies
#find /nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype/sequences/*.fa -type f > sequence_paths.txt
#cat sequence_paths.txt | cut -f10 -d"/" > sequence_ids.txt
#cat sequence_ids.txt | cut -f1 -d"." > sequence_ids_new.txt
#mv sequence_ids_new.txt sequence_ids.txt # clean up intermediary

# Use snippy to run samples
cat $outdir/sequence_ids.txt | while read line || [[ -n $line ]]; do snippy --cpus 47 --outdir "$outdir/snippy/results_$line" --ref /nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype/ATCC_700669_genbank.gbff --ctgs "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype/sequences/$line.contigs_velvet.fa" --force; done

# merge VCFs together
#cat $outdir/sequence_ids.txt | while read line || [[ -n $line ]]; do echo "$outdir/snippy/results_$line" >> $outdir/merge.txt; done
#awk 'NF{print $0 "/snps.vcf.gz"}' $outdir/merge.txt > $outdir/new.txt
#mv $outdir/new.txt $outdir/merge.txt
#bcftools merge -l $outdir/merge.txt -Oz --missing-to-ref -o $outdir/spn_snippy_merged.vcf.gz --threads 47

# fix formatting and filter snps
#bcftools +fixploidy spn_snippy_merged.vcf.gz -- -f 1 > spn_snippy_merged_haploid.vcf.gz
#bcftools filter -e'INFO/MAF<0.05' spn_snippy_merged_haploid.vcf.gz > spn_snippy_final.vcf.gz

# rename chroms to match snpeff database
bcftools annotate --rename-chrs chr.txt spn_snippy_final.vcf.gz | bgzip > snps_renamed_chr.vcf.gz

# run snpeff to get gene annotations & effect annotations
snpEff eff -v Streptococcus_pneumoniae_atcc_700669 -lof snps_renamed_chr.vcf.gz | bgzip > snps_effannotated.vcf.gz

# delete intermediary files
rm spn_snippy_merged.vcf.gz
rm spn_snippy_merged_haploid.vcf.gz
rm spn_snippy_final.vcf.gz
rm snps_renamed_chr.vcf.gz
