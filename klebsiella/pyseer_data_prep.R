################################################ LOAD LIBRARIES

library(dplyr)
library(data.table)

################################################ GENERATE DATA FOR PYSEER (single EF=100 continuous)

# read in genotype
genotype <- read.delim("/nfs/research/jlees/jacqueline/gwas_data/klebsiella/genotype/genotype_pres_abs.txt")
rownames(genotype) <- genotype[,1]
genotype <- genotype[,-1]
genotype <- t(genotype)

# read in lineages
lineages <- read.delim("/nfs/research/jlees/jacqueline/gwas_data/klebsiella/lineages/klebsiella_poppunk_clusters.txt", sep = "\t", header = FALSE)

# read in phenotype
phenotype <- read.delim("/nfs/research/jlees/jacqueline/gwas_data/klebsiella/phenotypes/simul_homoplasic_EF_100.0_herit_0.9_traitprev_0.1.phen", sep = " ", header = FALSE)

# Reorder samples in lineages file to match genotype and phenotype
sample_index <- match(phenotype[,2], lineages[,1])
lineages <- lineages[sample_index,]

# transpose back to normal format
genotype <- t(genotype)

# change phenotype header
phenotype <- phenotype[,2:3]
colnames(phenotype) <- c("sample", "continuous")

# turn row and column names to new row and column
genotype <- as.data.frame(genotype)
genotype <- cbind(rownames(genotype), genotype)
genotype <- rbind(colnames(genotype), genotype)
genotype[,1] <- c("variants", genotype[-1, 1])

# export genotypes as Rtab
write.table(genotype, file = "/nfs/research/jlees/jacqueline/gwas_data/klebsiella/pyseer_data/single_genotypes.rtab", 
            sep = "\t", col.names = F, row.names = F, quote = F)

# export phenotypes as tab separated
write.table(phenotype, file = "/nfs/research/jlees/jacqueline/gwas_data/klebsiella/pyseer_data/single_100_phenotypes.txt", 
            sep = "\t",
            col.names = T, row.names = F,
            quote = F)

# export sample list
write.table(phenotype[,1], file = "/nfs/research/jlees/jacqueline/gwas_data/klebsiella/pyseer_data/single_100_sample_list.txt", 
            sep = "\t", 
            col.names = F, row.names = F,
            quote = F)

################################################ GENERATE DATA FOR PYSEER (single EF=30 continuous)

# read in genotype
genotype <- read.delim("/nfs/research/jlees/jacqueline/gwas_data/klebsiella/genotype/genotype_pres_abs.txt")
rownames(genotype) <- genotype[,1]
genotype <- genotype[,-1]
genotype <- t(genotype)

# read in lineages
lineages <- read.delim("/nfs/research/jlees/jacqueline/gwas_data/klebsiella/lineages/klebsiella_poppunk_clusters.txt", sep = "\t", header = FALSE)

# read in phenotype
phenotype <- read.delim("/nfs/research/jlees/jacqueline/gwas_data/klebsiella/phenotypes/simul_homoplasic_EF_30.0_herit_0.9_traitprev_0.1.phen", sep = " ", header = FALSE)

# Reorder samples in lineages file to match genotype and phenotype
sample_index <- match(phenotype[,2], lineages[,1])
lineages <- lineages[sample_index,]

# change phenotype header
phenotype <- phenotype[,2:3]
colnames(phenotype) <- c("sample", "continuous")

# genotypes for single already written above (EF=100)

# export phenotypes as tab separated
write.table(phenotype, file = "/nfs/research/jlees/jacqueline/gwas_data/klebsiella/pyseer_data/single_30_phenotypes.txt", 
            sep = "\t",
            col.names = T, row.names = F,
            quote = F)

# export sample list
write.table(phenotype[,1], file = "/nfs/research/jlees/jacqueline/gwas_data/klebsiella/pyseer_data/single_30_sample_list.txt", 
            sep = "\t", 
            col.names = F, row.names = F,
            quote = F)

################################################ GENERATE DATA FOR PYSEER (oligogenic EF=100 continuous)

# read in genotype
genotype <- read.delim("/nfs/research/jlees/jacqueline/gwas_data/klebsiella/Kpneu_MIC_prediction/data/SNPs/matrix_and_lineages/pres_abs_matrix.txt")

# read in lineages
lineages <- read.delim("/nfs/research/jlees/jacqueline/gwas_data/klebsiella/lineages/klebsiella_poppunk_clusters.txt", sep = "\t", header = FALSE)

# read in phenotype
phenotype <- read.delim("/nfs/research/jlees/jacqueline/gwas_data/klebsiella/Kpneu_MIC_prediction/data/SNPs/simulated_traits/simul_oligogenic_EF_100.0_herit_0.9_traitprev_0.1.phen", sep = " ", header = FALSE)

# change phenotype header
phenotype <- phenotype[,2:3]
colnames(phenotype) <- c("sample", "continuous")

# turn row and column names to new row and column
genotype <- as.data.frame(genotype)
genotype <- rbind(colnames(genotype), genotype)
genotype[,1] <- c("variants", genotype[-1, 1])

# export genotypes as Rtab
write.table(genotype, file = "/nfs/research/jlees/jacqueline/gwas_data/klebsiella/pyseer_data/oligogenic_genotypes.rtab", 
            sep = "\t", col.names = F, row.names = F, quote = F)

# export phenotypes as tab separated
write.table(phenotype, file = "/nfs/research/jlees/jacqueline/gwas_data/klebsiella/pyseer_data/oligogenic_100_phenotypes.txt", 
            sep = "\t",
            col.names = T, row.names = F,
            quote = F)

# export sample list
write.table(phenotype[,1], file = "/nfs/research/jlees/jacqueline/gwas_data/klebsiella/pyseer_data/oligogenic_100_sample_list.txt", 
            sep = "\t", 
            col.names = F, row.names = F,
            quote = F)

################################################ GENERATE DATA FOR PYSEER (oligogenic EF=30 continuous)
# read in phenotype
phenotype <- read.delim("/nfs/research/jlees/jacqueline/gwas_data/klebsiella/Kpneu_MIC_prediction/data/SNPs/simulated_traits/simul_oligogenic_EF_30.0_herit_0.9_traitprev_0.1.phen", sep = " ", header = FALSE)

# change phenotype header
phenotype <- phenotype[,2:3]
colnames(phenotype) <- c("sample", "continuous")

# export phenotypes as tab separated
write.table(phenotype, file = "/nfs/research/jlees/jacqueline/gwas_data/klebsiella/pyseer_data/oligogenic_30_phenotypes.txt", 
            sep = "\t",
            col.names = T, row.names = F,
            quote = F)

# export sample list
write.table(phenotype[,1], file = "/nfs/research/jlees/jacqueline/gwas_data/klebsiella/pyseer_data/oligogenic_30_sample_list.txt", 
            sep = "\t", 
            col.names = F, row.names = F,
            quote = F)
