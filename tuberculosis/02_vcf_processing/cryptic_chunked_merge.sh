#!/usr/bin/env bash

#SBATCH --time=48:00:00   # walltime
#SBATCH --cpus-per-task=48   # number of CPUs Per Task
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem=40GB   # memory per node
#SBATCH --partition=datamover
#SBATCH --propagate=STACK # prevents bcftools 'too many open files' error
#SBATCH -J "merge_cryptic_vcfs"   # job name
#SBATCH -o "/nfs/research/jlees/jacqueline/gwas_data/code/gwas_data/tuberculosis/merge_cryptic_vcfs_%A_%a.out"   # job output file
#SBATCH -e "/nfs/research/jlees/jacqueline/gwas_data/code/gwas_data/tuberculosis/merge_cryptic_vcfs_%A_%a.err"   # job error file


# specify path
outdir=/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/cryptic_regenotyped_merge
data_dir=/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/cryptic_regeno_snpeff

filelist="${outdir}/../file.txt"

chrom="NC_000962.3"
genome_len=4411532
n_chunks=20

threads=${SLURM_CPUS_PER_TASK}
outprefix="${outdir}/cryptic_merged_regenotyped"

task_id="${SLURM_ARRAY_TASK_ID}"

chunk_size=$(( genome_len / n_chunks ))

# 1-based inclusive start
start=$(( (task_id - 1) * chunk_size + 1 ))

# inclusive end; last chunk takes the remainder
if (( task_id == n_chunks )); then
  end=$genome_len
else
  end=$(( task_id * chunk_size ))
fi

region="${chrom}:${start}-${end}"
outvcf="${outprefix}.chunk${task_id}.vcf.gz"

#mkdir -p "${outdir}" logs

#bcftools merge \
#  --file-list "${filelist}" \
#  -m none \
#  -r "${region}" \
#  -Oz \
#  --threads "${threads}" \
#  -o "${outvcf}"

#bcftools index -f "${outvcf}"

cd $data_dir

# rename chromosome to match snpeff database
#bcftools annotate --rename-chrs chr.txt ../cryptic_regenotyped_merge/cryptic_merged_regenotyped.chunk${task_id}.vcf.gz --threads $threads -Oz --output $data_dir/cryptic_regeno_renamed.chunk${task_id}.vcf.gz 
#bcftools index cryptic_regeno_renamed.chunk${task_id}.vcf.gz

# check number of samples in VCF
#bcftools head cryptic_regeno_renamed.chunk${task_id}.vcf.gz | grep -o site | wc -l

# fix ploidy (diploid -> haploid)
#bcftools +fixploidy cryptic_regeno_renamed.chunk${task_id}.vcf.gz -- -f 1 > $data_dir/cryptic_regeno.chunk${task_id}.vcf
#bgzip $data_dir/cryptic_regeno.chunk${task_id}.vcf
#bcftools index cryptic_regeno.chunk${task_id}.vcf.gz

# check number of samples in VCF
#bcftools head cryptic_regeno.chunk${task_id}.vcf.gz | grep -o site | wc -l

# add MAF column
#bcftools +fill-tags cryptic_regeno.chunk${task_id}.vcf.gz --threads $threads -- -t MAF > $data_dir/cryptic_regeno_addtags.chunk${task_id}.vcf
#bgzip $data_dir/cryptic_regeno_addtags.chunk${task_id}.vcf
#bcftools index $data_dir/cryptic_regeno_addtags.chunk${task_id}.vcf.gz

# filter for >0.005 MAF
#bcftools filter -i'INFO/MAF>0.0004' cryptic_regeno_addtags.chunk${task_id}.vcf.gz --threads $threads -Oz --output $data_dir/cryptic_regeno_MAF_0.0004.chunk${task_id}.vcf.gz
#bcftools index $data_dir/cryptic_regeno_MAF_0.0004.chunk${task_id}.vcf.gz # 87 mutations in rpoB, at least 5 samples per variant

# merge MAF-filtered (much smaller) VCFs together

#ls -v "${data_dir}/cryptic_regeno_MAF_0.0004.chunk*.vcf.gz" > chunks_to_merge.txt


#bcftools concat \
#  --file-list $data_dir/chunks_to_merge.txt \
#  -Oz \
#  --threads "${threads}" \
#  -o "${data_dir}/cryptic_regeno_MAF_0.0004.allchunks.vcf.gz"

#bcftools index -f "${data_dir}/cryptic_regeno_MAF_0.0004.allchunks.vcf.gz"

#gunzip -c $data_dir/cryptic_regeno_MAF_0.0004.allchunks.vcf.gz > $data_dir/cryptic_regeno_MAF_0.0004.allchunks.vcf

# in bcftools env:

#snpEff eff -v Mycobacterium_tuberculosis_h37rv -lof $data_dir/cryptic_regeno_MAF_0.0004.allchunks.vcf | bgzip > $data_dir/cryptic_regeno_MAF_0.0004.allchunks_effannotated.vcf.gz # (-lof adds loss of function tags)

# extract fields with one line per effect

#gunzip -c $data_dir/cryptic_regeno_MAF_0.0004.allchunks_effannotated.vcf.gz > $data_dir/cryptic_regeno_MAF_0.0004.allchunks_effannotated.vcf

cat $data_dir/cryptic_regeno_MAF_0.0004.allchunks_effannotated.vcf | /hps/software/users/jlees/jacqueline/pkgs/snpsift-5.2-hdfd78af_0/share/snpsift-5.2-0/scripts/vcfEffOnePerLine.pl | SnpSift extractFields -  POS REF ALT "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "LOF[*].GENE" > fields_all.txt

cat $data_dir/cryptic_regeno_MAF_0.0004.allchunks_effannotated.vcf |  /hps/software/users/jlees/jacqueline/pkgs/snpsift-5.2-hdfd78af_0/share/snpsift-5.2-0/scripts/vcfEffOnePerLine.pl | SnpSift extractFields -  POS REF ALT "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "LOF[*].GENE" | sed '/MODIFIER/d' > fields_filtered.txt

# fix ref/alt calls (missings to 0, multiallelics to 1)
#bcftools +setGT $data_dir/cryptic_regeno_MAF_0.0004.allchunks_effannotated.vcf.gz -Oz -o $data_dir/cryptic_regeno_MAF_0.0004.allchunks_effannotated_tmp.vcf.gz -- -t q -n c:"0" -i 'GT="."'
#bcftools index $data_dir/cryptic_regeno_MAF_0.0004.allchunks_effannotated_tmp.vcf.gz

#bcftools +setGT $data_dir/cryptic_regeno_MAF_0.0004.allchunks_effannotated_tmp.vcf.gz -Oz -o $data_dir/cryptic_regeno_MAF_0.0004.allchunks_effannotated_binary.vcf.gz -- -t q -n c:"1" -i 'GT~"[1-9]"'
#bcftools index $data_dir/cryptic_regeno_MAF_0.0004.allchunks_effannotated_binary.vcf.gz 

#rm $data_dir/cryptic_regeno_MAF_0.0004.allchunks_effannotated_tmp.vcf.gz

