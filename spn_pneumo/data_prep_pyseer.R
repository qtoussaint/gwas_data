################################################### LOAD LIBRARIES

library("vcfR")
library("dplyr")
library("data.table")

################################################### GENERATE DATA FOR PYSEER

# read in genotype 
sparc_snps <- read.vcfR("/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype/sparc_snps.vcf.gz", verbose = TRUE)
genotype <- extract.gt(sparc_snps, as.numeric = TRUE, IDtoRowNames = T)

# read in lineages
lineages <- read.csv("/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/lineages/database_clusters.csv")

# read in phenotype data
phenotypes <- read.csv("/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/phenotype/sparc_metadata.csv")

# select desired phenotype (here it's benzylpenicillin MICs)
penicillin_mics <- as.data.frame(phenotypes$Benzylpenicillin.MIC..ug.mL.)
penicillin_mics <- data.frame(MIC=penicillin_mics[,1], ID=phenotypes$Lane.name)

# ensure order/number of samples matches between these 3 data sets
# sort/filter those that dont
index <- match(dimnames(genotype)[[2]], penicillin_mics$ID)
penicillin_mics <- penicillin_mics[index,]

index <- match(dimnames(genotype)[[2]], lineages$Taxon)
lineages <- lineages[index,]

# transpose genotype data (need samples in rows and variants in columns)
genotype <- t(genotype)

# find and replace missing data in genotype file
# will use change all missing SNPs to 0
genotype[is.na(genotype)] <- 0

# filter for MAF>0.05 (if desired)
small_maf <- which((colSums(genotype)/nrow(genotype))<0.05)
genotype <- genotype[,-small_maf]

# Clean up the MIC values

# show original values
print("True MIC values:")
print(sort(unique(penicillin_mics$MIC)))

# clean duplicate values and remove symbols
penicillin_mics$MIC <- replace(penicillin_mics$MIC, penicillin_mics$MIC==">4.00", 4)
penicillin_mics$MIC <- replace(penicillin_mics$MIC, penicillin_mics$MIC==">4", 4)
penicillin_mics$MIC <- replace(penicillin_mics$MIC, penicillin_mics$MIC=="<=0.03", 0.03)
penicillin_mics$MIC <- replace(penicillin_mics$MIC, penicillin_mics$MIC=="<0.03", 0.03)
penicillin_mics$MIC <- replace(penicillin_mics$MIC, penicillin_mics$MIC=="<0.016", 0.016)
penicillin_mics$MIC <- replace(penicillin_mics$MIC, penicillin_mics$MIC=="<.016", 0.016)
penicillin_mics$MIC <- replace(penicillin_mics$MIC, penicillin_mics$MIC=="< 0.01", 0.01)

# make numeric
penicillin_mics$MIC <- as.numeric(penicillin_mics$MIC)

# test different bin breaks and make the cuts
# 13 bins (from 26 total)
MICs <- sort(unique(penicillin_mics$MIC))
print(unique(cut(MICs, breaks=c(0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 1, 2, 4, 6))))

# cut according to give breakpoints and add as separate column
penicillin_complete <- penicillin_mics %>% mutate(MIC_factor = cut(penicillin_mics$MIC, breaks=c(0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 1, 2, 4, 6)))

# replace factors with characters
penicillin_complete$MIC_factor <- print(penicillin_complete$MIC_factor)

# create sorted vector of all MIC bins
MIC_intervals <- as.character(sort(unique(penicillin_complete$MIC_factor)))

# convert MIC to sample,s ordinal category (1:24) in resistances
for (i in 1:nrow(penicillin_complete)) {
  penicillin_complete$MIC_ord_category[i] <- which(MIC_intervals == penicillin_complete$MIC_factor[i])
}

# switch order of columns and rename
penicillin_mics <- data.frame(ID=penicillin_complete$ID, K=penicillin_complete$MIC_ord_category)

# change from category to center of break (for pyseer)

penicillin_mics$K <- replace(penicillin_mics$K, penicillin_mics$K==1, 0.01)
penicillin_mics$K <- replace(penicillin_mics$K, penicillin_mics$K==2, 0.03)
penicillin_mics$K <- replace(penicillin_mics$K, penicillin_mics$K==3, 0.05)
penicillin_mics$K <- replace(penicillin_mics$K, penicillin_mics$K==4, 0.07)
penicillin_mics$K <- replace(penicillin_mics$K, penicillin_mics$K==5, 0.09)
penicillin_mics$K <- replace(penicillin_mics$K, penicillin_mics$K==6, 0.15)
penicillin_mics$K <- replace(penicillin_mics$K, penicillin_mics$K==7, 0.25)
penicillin_mics$K <- replace(penicillin_mics$K, penicillin_mics$K==8, 0.35)
penicillin_mics$K <- replace(penicillin_mics$K, penicillin_mics$K==9, 0.45)
penicillin_mics$K <- replace(penicillin_mics$K, penicillin_mics$K==10, 0.75)
penicillin_mics$K <- replace(penicillin_mics$K, penicillin_mics$K==11, 1.5)
penicillin_mics$K <- replace(penicillin_mics$K, penicillin_mics$K==12, 3)
penicillin_mics$K <- replace(penicillin_mics$K, penicillin_mics$K==13, 5)

# change phenotype header
colnames(penicillin_mics) <- c("sample", "continuous")

# export genotypes as Rtab
write.table(genotype, file = "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/pyseer_data/genotypes.rtab", 
            sep = "\t", col.names = T, row.names = T, quote = F)

# export phenotypes as tab separated
write.table(penicillin_mics, file = "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/pyseer_data/phenotypes_pen_mic.txt", 
            sep = "\t",
            col.names = T, row.names = F,
            quote = F)

# export sample list
write.table(penicillin_mics[,1], file = "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/pyseer_data/sample_list.txt", 
            sep = "\t", 
            col.names = F, row.names = F, quote = F)
