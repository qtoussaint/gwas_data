#!/usr/bin/env bash

#SBATCH --time=24:00:00             # walltime
#SBATCH --ntasks=1                  # number of tasks
#SBATCH --cpus-per-task=5          # number of CPUs Per Task
#SBATCH --nodes=1                   # number of nodes
#SBATCH --mem=5GB                  # memory per node
#SBATCH -J "snippy"             # job name
#SBATCH -o "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype/snippy_bcftools_snpeff.out"         # job output file
#SBATCH -e "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype/snippy_bcftools_snpeff.err"         # job error file

#################################################################################

outdir=/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype

cd $outdir

# Make a list of all assemblies
#find /nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype/sequences/*.fa -type f > sequence_paths.txt
#cat sequence_paths.txt | cut -f10 -d"/" > sequence_ids.txt
#cat sequence_ids.txt | cut -f1 -d"." > sequence_ids_new.txt
#mv sequence_ids_new.txt sequence_ids.txt # clean up intermediary

# Use snippy to run samples
#cat $outdir/sequence_ids.txt | while read line || [[ -n $line ]]; do snippy --cpus 47 --outdir "$outdir/snippy/results_$line" --ref /nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype/ATCC_700669_genbank.gbff --ctgs "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype/sequences/$line.contigs_velvet.fa" --force; done

# merge VCFs together
#cat $outdir/sequence_ids.txt | while read line || [[ -n $line ]]; do echo "$outdir/snippy/results_$line" >> $outdir/merge.txt; done
#awk 'NF{print $0 "/snps.vcf.gz"}' $outdir/merge.txt > $outdir/new.txt
#mv $outdir/new.txt $outdir/merge.txt
#bcftools merge -l $outdir/merge.txt -Oz --missing-to-ref -o $outdir/spn_snippy_merged.vcf.gz --threads 47

# fix formatting and filter snps
#bcftools +fixploidy spn_snippy_merged.vcf.gz -- -f 1 > spn_snippy_merged_haploid.vcf.gz
#bcftools filter -e'INFO/MAF<0.01' spn_snippy_merged_haploid.vcf.gz > spn_snippy_final.vcf.gz

# rename chroms to match snpeff database
#bcftools annotate --rename-chrs chr.txt spn_snippy_final.vcf.gz | bgzip > snps_renamed_chr.vcf.gz

### KEEP MULTIALLELIC

# unzip
#gunzip -c snps_renamed_chr.vcf.gz > snps_renamed_chr_unzip.vcf

# run snpeff to get gene annotations & effect annotations
#snpEff eff -v Streptococcus_pneumoniae_atcc_700669 -lof snps_renamed_chr_unzip.vcf | bgzip > snps_snippy_maf01_multiallelic.vcf.gz

# unzip
#gunzip -c snps_snippy_maf01_multiallelic.vcf.gz > snps_snippy_maf01_multiallelic_unzip.vcf

# field extraction with snpEff
#cat snps_snippy_maf01_multiallelic_unzip.vcf | /hps/software/users/jlees/jacqueline/pkgs/snpsift-5.2-hdfd78af_0/share/snpsift-5.2-0/scripts/vcfEffOnePerLine.pl | SnpSift extractFields -  POS REF ALT "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "LOF[*].GENE" > fields_all_multiallelic.txt
#cat snps_snippy_maf01_multiallelic_unzip.vcf |  /hps/software/users/jlees/jacqueline/pkgs/snpsift-5.2-hdfd78af_0/share/snpsift-5.2-0/scripts/vcfEffOnePerLine.pl | SnpSift extractFields -  POS REF ALT "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "LOF[*].GENE" | sed '/MODIFIER/d' > fields_filtered_multiallelic.txt

# extract presence/absence matrix with biallelic (no 2,3,4+ allele markings, just 0 or 1)
# probably best in case snippy annotations aren't as good and end up with many low-frequency variants
#bcftools query -f '%CHROM\_%POS\_%REF\_%ALT[\t%GT]\n' snps_snippy_maf01_multiallelic.vcf.gz > genotype_multiallelic_tmp.tsv

#awk '{
#    for (i=2; i<=NF; i++) {
#        if ($i == "0" || $i == "0") {
#            $i = 0;
#        } else {
#            $i = 1;
#        }
#    }
#    print
#}' genotype_multiallelic_tmp.tsv > genotype_multiallelic_tmp2.tsv

# get header with samples
#bcftools query -l snps_snippy_maf01_multiallelic.vcf.gz > samples.txt
#echo -e "CHROM\_POS\_REF\_ALT$(awk '{printf "\t%s", $0}' samples.txt)" > header.tsv
#cat header.tsv genotype_multiallelic_tmp2.tsv > genotype_multiallelic.tsv

### MULTIALLELIC WITH MAF=0.05 INCL MODIFIER-ONLY SITES, MULTIALLELIC WITH MAF=0.05 AND !!! NO !!! MODIFIER-ONLY SITES

bcftools filter -e'INFO/MAF<0.05' snps_snippy_maf01_multiallelic.vcf.gz > snps_snippy_maf05_multiallelic.vcf.gz

# unzip
gunzip -c snps_snippy_maf05_multiallelic.vcf.gz > snps_snippy_maf05_multiallelic_unzip.vcf

# field extraction with snpEff
cat snps_snippy_maf05_multiallelic_unzip.vcf | /hps/software/users/jlees/jacqueline/pkgs/snpsift-5.2-hdfd78af_0/share/snpsift-5.2-0/scripts/vcfEffOnePerLine.pl | SnpSift extractFields -  POS REF ALT "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "LOF[*].GENE" > fields_all_maf05_multiallelic.txt
cat snps_snippy_maf05_multiallelic_unzip.vcf |  /hps/software/users/jlees/jacqueline/pkgs/snpsift-5.2-hdfd78af_0/share/snpsift-5.2-0/scripts/vcfEffOnePerLine.pl | SnpSift extractFields -  POS REF ALT "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "LOF[*].GENE" | sed '/MODIFIER/d' > fields_filtered_maf05_multiallelic.txt

# extract presence/absence matrix with biallelic (no 2,3,4+ allele markings, just 0 or 1)
# probably best in case snippy annotations aren't as good and end up with many low-frequency variants
bcftools query -f '%CHROM\_%POS\_%REF\_%ALT[\t%GT]\n' snps_snippy_maf05_multiallelic.vcf.gz > genotype_maf05_multiallelic_tmp.tsv

awk '{
    for (i=2; i<=NF; i++) {
        if ($i == "0" || $i == "0") {
            $i = 0;
        } else {
            $i = 1;
        }
    }
    print
}' genotype_maf05_multiallelic_tmp.tsv > genotype_maf05_multiallelic_tmp2.tsv

# get header with samples
bcftools query -l snps_snippy_maf05_multiallelic.vcf.gz > samples.txt
echo -e "CHROM\_POS\_REF\_ALT$(awk '{printf "\t%s", $0}' samples.txt)" > header.tsv
cat header.tsv genotype_maf05_multiallelic_tmp2.tsv > genotype_maf05_multiallelic.tsv

# NOW FILTER OUT MODIFIER-ONLY VARIANTS

# make bed file that only contains positions with at least one non-modifier variant 
awk 'NR>1 {print "Chromosome", $1}' fields_filtered_maf05_multiallelic.txt > positions_filtered_maf05_multiallelic.bed

# filter vcf based on this bed
bcftools view -T positions_filtered_maf05_multiallelic.bed snps_snippy_maf05_multiallelic.vcf.gz > snps_snippy_maf05_NOMODIFIERS_multiallelic.vcf.gz

# make genotype matrix as before
bcftools query -f '%CHROM\_%POS\_%REF\_%ALT[\t%GT]\n' snps_snippy_maf05_NOMODIFIERS_multiallelic.vcf.gz > genotype_maf05_NOMODIFIERS_multiallelic_tmp.tsv

awk '{
    for (i=2; i<=NF; i++) {
        if ($i == "0" || $i == "0") {
            $i = 0;
        } else {
            $i = 1;
        }
    }
    print
}' genotype_maf05_NOMODIFIERS_multiallelic_tmp.tsv > genotype_maf05_NOMODIFIERS_multiallelic_tmp2.tsv

# get header with samples
bcftools query -l snps_snippy_maf05_NOMODIFIERS_multiallelic.vcf.gz > samples.txt
echo -e "CHROM\_POS\_REF\_ALT$(awk '{printf "\t%s", $0}' samples.txt)" > header.tsv
cat header.tsv genotype_maf05_NOMODIFIERS_multiallelic_tmp2.tsv > genotype_maf05_NOMODIFIERS_multiallelic.tsv

# FILTER 'LOW' IMPACT (I THINK THIS ONLY MEANS SYNONYMOUS?) VARIANTS AS WELL



### SPLIT MULTIALLELIC

# split multiallelics:
#bcftools norm --multiallelics -any snps_renamed_chr.vcf.gz > snps_biallelic.vcf

# run snpeff to get gene annotations & effect annotations
#snpEff eff -v Streptococcus_pneumoniae_atcc_700669 -lof snps_biallelic.vcf | bgzip > snps_snippy_maf01_split_multialleic.vcf.gz

# unzip
#gunzip -c snps_snippy_maf01_split_multialleic.vcf.gz > snps_snippy_maf01_split_multialleic_unzipped.vcf

# field extraction with snpEff
#cat snps_snippy_maf01_split_multialleic_unzipped.vcf | /hps/software/users/jlees/jacqueline/pkgs/snpsift-5.2-hdfd78af_0/share/snpsift-5.2-0/scripts/vcfEffOnePerLine.pl | SnpSift extractFields -  POS REF ALT "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "LOF[*].GENE" > fields_all_split_multiallelic.txt
#cat snps_snippy_maf01_split_multialleic_unzipped.vcf |  /hps/software/users/jlees/jacqueline/pkgs/snpsift-5.2-hdfd78af_0/share/snpsift-5.2-0/scripts/vcfEffOnePerLine.pl | SnpSift extractFields -  POS REF ALT "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "LOF[*].GENE" | sed '/MODIFIER/d' > fields_filtered_split_multiallelic.txt

# extract presence/absence matrix with one line per allele
# more granular but could result in many low-frequency variants at a site and reduce power
#bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' snps_snippy_maf01_multiallelic.vcf.gz > genotype_split_multiallelic.tsv

### CLEAN UP

# delete intermediary files
rm spn_snippy_merged.vcf.gz
rm spn_snippy_merged_haploid.vcf.gz
rm spn_snippy_final.vcf.gz
rm snps_renamed_chr_unzip.vcf
rm snps_renamed_chr.vcf.gz
rm snps_biallelic.vcf
rm snps_snippy_maf01_split_multialleic.vcf
rm snps_snippy_maf01_multialleic.vcf
rm snps_snippy_maf01_split_multialleic_unzipped.vcf
rm snps_snippy_maf01_multiallelic_unzip.vcf
rm genotype_multiallelic_tmp.tsv
rm samples.txt
rm header.tsv
rm genotype_multiallelic_tmp2.tsv
rm snps_snippy_maf05_multiallelic_unzip.vcf
rm genotype_maf05_multiallelic_tmp.tsv
rm genotype_maf05_multiallelic_tmp2.tsv
rm positions_filtered_maf05_multiallelic.bed
rm genotype_maf05_NOMODIFIERS_multiallelic_tmp.tsv
rm genotype_maf05_NOMODIFIERS_multiallelic_tmp2.tsv
