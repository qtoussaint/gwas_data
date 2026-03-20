#!/usr/bin/env bash

#SBATCH --time=48:00:00   # walltime
#SBATCH --ntasks=1   # number of tasks
#SBATCH --cpus-per-task=4   # number of CPUs Per Task
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem=20GB   # memory per node
#SBATCH --partition=research
#SBATCH -J "snpEff"   # job name
#SBATCH -o "snpEff.out"   # job output file
#SBATCH -e "snpEff.err"   # job error file

# specify path
outdir=/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/genotype/snpEff
data_dir=/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/genotype

############

cd $outdir

# after looking at database, discovered "CHROM" column contains "Chromosome", not "26":
#snpEff download Mycobacterium_tuberculosis_h37rv -dataDir $outdir
#cd Mycobacterium_tuberculosis_h37rv
#gunzip -c snpEffectPredictor.bin | head
#cd ..

# renamed "NC_000962.3" column in VCF to "Chromosome" like it is in the database:
#bcftools view $/data_dir/tb_mhunt_filtered.vcf.gz | cut -f 1  | head -40 # check contents of existing chromosome column
#echo "NC_000962.3 Chromosome" > chr.txt
#bcftools annotate --rename-chrs chr.txt $data_dir/tb_mhunt_filtered.vcf.gz | bgzip > tb_snps_renamed.vcf.gz
#gunzip -c tb_snps_renamed.vcf.gz > tb_snps_renamed.vcf

# in bcftools env:

snpEff eff -v Mycobacterium_tuberculosis_h37rv -lof $outdir/.vcf | bgzip > $outdir/snps_effannotated.vcf.gz # (-lof adds loss of function tags)

# extract fields with one line per effect

gunzip -c snps_effannotated.vcf.gz > snps_effannotated.vcf

cat snps_effannotated.vcf | /hps/software/users/jlees/jacqueline/pkgs/snpsift-5.2-hdfd78af_0/share/snpsift-5.2-0/scripts/vcfEffOnePerLine.pl | SnpSift extractFields -  POS REF ALT "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "LOF[*].GENE" > fields_all.txt

cat snps_effannotated.vcf |  /hps/software/users/jlees/jacqueline/pkgs/snpsift-5.2-hdfd78af_0/share/snpsift-5.2-0/scripts/vcfEffOnePerLine.pl | SnpSift extractFields -  POS REF ALT "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "LOF[*].GENE" | sed '/MODIFIER/d' > fields_filtered.txt
