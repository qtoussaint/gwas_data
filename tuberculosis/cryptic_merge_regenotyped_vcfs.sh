#!/usr/bin/env bash

#SBATCH --time=10:00:00   # walltime
#SBATCH --ntasks=1   # number of tasks
#SBATCH --cpus-per-task=20   # number of CPUs Per Task
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem=100GB   # memory per node
#SBATCH --partition=datamover
#SBATCH --propagate=STACK # prevents bcftools 'too many open files' error
#SBATCH -J "merge_cryptic_vcfs"   # job name
#SBATCH -o "/nfs/research/jlees/jacqueline/gwas_data/code/gwas_data/tuberculosis/merge_cryptic_vcfs.out"   # job output file
#SBATCH -e "/nfs/research/jlees/jacqueline/gwas_data/code/gwas_data/tuberculosis/merge_cryptic_vcfs.err"   # job error file

#####################

# datamover partition allows access to cryptic files on /nfs/ftp

# specify path
outdir=/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/cryptic_regenotyped_merge
data_dir=/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/cryptic_regeno_snpeff

# check that there are SNPs in the rpoB region
#bcftools query -f'%CHROM  %POS  %REF  %ALT{0}\n'  -r NC_000962.3:763370-767320 site.02.iso.1.subject.0010.lab_id.2014222040.seq_reps.14222040.regeno.vcf.gz

#bcftools merge --file-list $outdir/../file.txt -m none -O z -o $outdir/cryptic_merged_regenotyped.vcf.gz
#bcftools index $outdir/cryptic_merged_regenotyped.vcf.gz

cd $data_dir

# rename chromosome to match snpeff database
#bcftools annotate --rename-chrs chr.txt ../cryptic_regenotyped_merge/cryptic_merged_regenotyped_bgzip.vcf.gz --threads 20 --output $data_dir/cryptic_regeno_renamed.vcf.gz
#bcftools index cryptic_regeno_renamed.vcf.gz

# check number of samples in VCF
#bcftools head cryptic_regeno_renamed.vcf.gz | grep -o site | wc -l

# fix ploidy (diploid -> haploid)
#bcftools +fixploidy cryptic_regeno_renamed.vcf.gz -- -f 1 > $data_dir/cryptic_regeno.vcf
#bgzip $data_dir/cryptic_regeno.vcf
#bcftools index cryptic_regeno.vcf.gz

# check number of samples in VCF
#bcftools head cryptic_regeno.vcf.gz | grep -o site | wc -l

# add MAF column
#bcftools +fill-tags cryptic_regeno.vcf.gz --threads 20 -- -t MAF > $data_dir/cryptic_regeno_addtags.vcf
#bgzip $data_dir/cryptic_regeno_addtags.vcf
#bcftools index $data_dir/cryptic_regeno_addtags.vcf.gz 

# filter for >0.005 MAF
#bcftools filter -i'INFO/MAF>0.0004' cryptic_regeno_addtags.vcf.gz --threads 20 --output $data_dir/cryptic_regeno_MAF_0.0004.vcf.gz
#bcftools index cryptic_regeno_MAF_0.0004.vcf.gz # 87 mutations in rpoB, at least 5 samples per variant

# fix ref/alt calls (missings to 0, multiallelics to 1)
bcftools +setGT cryptic_regeno_MAF_0.0004_effannotated.vcf.gz  -Oz -o cryptic_regeno_MAF_0.0004_effannotated_tmp.vcf.gz -- -t q -n c:"0" -i 'GT="."'
bcftools index cryptic_regeno_MAF_0.0004_effannotated_tmp.vcf.gz

bcftools +setGT cryptic_regeno_MAF_0.0004_effannotated_tmp.vcf.gz  -Oz -o cryptic_regeno_MAF_0.0004_effannotated_binary.vcf.gz -- -t q -n c:"1" -i 'GT~"[1-9]"'
bcftools index cryptic_regeno_MAF_0.0004_effannotated_binary.vcf.gz

rm cryptic_regeno_MAF_0.0004_effannotated_tmp.vcf.gz*

