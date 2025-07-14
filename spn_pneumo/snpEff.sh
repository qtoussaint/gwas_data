#!/usr/bin/env bash

#SBATCH --time=48:00:00   # walltime
#SBATCH --ntasks=1   # number of tasks
#SBATCH --cpus-per-task=17   # number of CPUs Per Task
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem=150GB   # memory per node
#SBATCH -J "snpEff"   # job name
#SBATCH -o "snpEff.out"   # job output file
#SBATCH -e "snpEff.err"   # job error file

# specify path
outdir=/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/snpEff
data_dir=/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/pyseer_data/pyseer_tutorial

############

cd $outdir

# after looking at database, discovered "CHROM" column contains "Chromosome", not "26":
#snpEff download Streptococcus_pneumoniae_atcc_700669 -dataDir $outdir
#cd Streptococcus_pneumoniae_atcc_700669
#gunzip -c snpEffectPredictor.bin | head

# renamed "26" column in VCF to "Chromosome" like it is in the database:
#echo "26 Chromosome" > chr.txt
#bcftools annotate --rename-chrs chr.txt $data_dir/snps.vcf.gz | bgzip > snps_renamed.vcf.gz

# in bcftools env:

snpEff eff Streptococcus_pneumoniae_atcc_700669 -lof -t snps_renamed.vcf.gz > snps_effannotated.vcf.gz
# (-lof adds loss of function tags)

# extract fields with one line per effect
cat snps_effannotated.vcf | /hps/software/users/jlees/jacqueline/pkgs/snpsift-5.2-hdfd78af_0/share/snpsift-5.2-0/scripts/vcfEffOnePerLine.pl | SnpSift extractFields -  POS REF ALT "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "LOF[*].GENE" > fields_all.txt

cat snps_effannotated.vcf |  /hps/software/users/jlees/jacqueline/pkgs/snpsift-5.2-hdfd78af_0/share/snpsift-5.2-0/scripts/vcfEffOnePerLine.pl | SnpSift extractFields -  POS REF ALT "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "LOF[*].GENE" | sed '/MODIFIER/d' > fields_filtered.txt

# realized I need to split multiallelics:
bcftools norm --multiallelics -any snps_effannotated.vcf > snps_biallelic.vcf
gzip -c snps_biallelic.vcf > snps_biallelic.vcf.gz
