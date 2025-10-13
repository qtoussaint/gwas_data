################################################### LOAD LIBRARIES

library("dplyr")
library("data.table")

###################################################

# read in variants
genotype <- fread("/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/cryptic_regeno_snpeff/cryptic_regeno_MAF_0.0004_effannotated.vcf.gz")

# remove variant metadata
variants <- genotype[,c(1:9)] # retain variant ids
genotype <- genotype[,-c(1:9)]
samples <- colnames(genotype) # retain sample ids

# match regenotyped sample IDs to generic sample IDs
phenotype <- read.csv(file = "/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/CRyPTIC_reuse_table_20240917.csv")
phenotype_samples <- basename(phenotype$REGENOTYPED_VCF)
phenotype_samples <- sub("\\.regeno\\.vcf\\.gz$", "", phenotype_samples)
match(samples, phenotype_samples)

# above confirms the order was preserved, so can copy IDs w/o matching
#colnames(genotype) <- phenotype$ENA_RUN

# remove VCF formatting from variant calls
genotype <- as.matrix(genotype[, 1:ncol(genotype)])
genotype <- sub(":.*", "", genotype)

# fix ref/alt calls
genotype[genotype == "."] <- "0"
num_vals <- as.numeric(genotype)
genotype[!is.na(num_vals) & num_vals >= 1 & num_vals <= 1000] <- "1"
genotype[] <- lapply(genotype, as.numeric)

# now add variants to rownames
variant_id <- paste0(variants$`#CHROM`, "_", variants$POS, "_", variants$REF, "_", variants$ALT)
genotype <- cbind(variant_id, genotype)

# write to TSV file (for faster reading in next time)
fwrite(genotype,
       "/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/cryptic_regeno_snpeff/cryptic_regeno_MAF_0.0004_effannotated_presabs.tsv",
       sep = "\t", col.names = F, row.names = F, quote = F)
