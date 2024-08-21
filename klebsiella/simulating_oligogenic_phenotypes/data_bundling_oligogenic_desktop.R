### Pipeline for continuous GWAS model

######## Oligogenic data bundling (original causal SNPs)

# load libraries
library(cmdstanr)

# read in genotype (SNPs)
genotype <- read.delim("~/Desktop/simulating_klebsiella_data/pres_abs_matrix.txt")
rownames(genotype) <- genotype[,1]
genotype <- genotype[,-1]
genotype <- t(genotype)

# read in lineages
lineages <- read.delim("~/Desktop/simulating_klebsiella_data/PP_clusters.txt", sep = "\t", header = FALSE)

# read in phenotype
phenotype <- read.delim("~/Desktop/simulating_klebsiella_data/sim_kleb_SNPs_EF_100_H_09.phen", sep = " ", header = FALSE)

# Reorder samples in lineages file to match genotype and phenotype
sample_index <- match(phenotype[,2], lineages[,1])
lineages <- lineages[sample_index,]

### Concat into one file (to check for correctness)

# define V = number of variants, L = number of lineages, N = number of samples
V <- as.numeric(ncol(genotype))
#L <- as.numeric(max(lineages[,2]))
L <- 347
N <- as.numeric(nrow(phenotype))

print(V)
print(L)
print(N)

checkme <- cbind(genotype[,1], phenotype, lineages)

# convert lineages into one-hot format
lineages <- model.matrix(lineages[,1]~as.factor(lineages[,2]), lineages)[,-1]

# make data into file
cleaned_gwas_data <- list(phenotype = phenotype[,3],
                          variant_matrix = genotype,
                          lineage_matrix = lineages,
                          N=N,
                          V=V,
                          L=L)

# save completed bundle, ready for gwas, to the appropriate folder
# we use stan json format as other formats (such as RDump) aren't supported by cmdstan

json_file <- tempfile(fileext = ".json", tmpdir = "~/Desktop/simulating_klebsiella_data/", pattern = "oligogenic_EF_100_H_0_9_associative_")
cmdstanr::write_stan_json(data = cleaned_gwas_data, file = json_file)

# create index of variant positions
index <- as.matrix(colnames(genotype))
colnames(index) <- c("variant")
write.csv(x = index,
          file = "~/Desktop/simulating_klebsiella_data/oligogenic_EF_100_H_0_9_associative_variant_index.csv")

############################################################################## OLIGO EF=100 N_causal_SNPs=8

# read in genotype (SNPs)
genotype <- read.delim("~/Desktop/simulating_klebsiella_data/pres_abs_matrix.txt")
rownames(genotype) <- genotype[,1]
genotype <- genotype[,-1]
genotype <- t(genotype)

# read in lineages
lineages <- read.delim("~/Desktop/simulating_klebsiella_data/PP_clusters.txt", sep = "\t", header = FALSE)

# read in phenotype
phenotype <- read.delim("~/Desktop/simulating_klebsiella_data/EF_100_ncausal_8/EF_100_ncausal_8.phen", sep = " ", header = FALSE)

# Reorder samples in lineages file to match genotype and phenotype
sample_index <- match(phenotype[,2], lineages[,1])
lineages <- lineages[sample_index,]

### Concat into one file (to check for correctness)

# define V = number of variants, L = number of lineages, N = number of samples
V <- as.numeric(ncol(genotype))
#L <- as.numeric(max(lineages[,2]))
L <- 347
N <- as.numeric(nrow(phenotype))

print(V)
print(L)
print(N)

checkme <- cbind(genotype[,1], phenotype, lineages)

# convert lineages into one-hot format
lineages <- model.matrix(lineages[,1]~as.factor(lineages[,2]), lineages)[,-1]

# make data into file
cleaned_gwas_data <- list(phenotype = phenotype[,3],
                          variant_matrix = genotype,
                          lineage_matrix = lineages,
                          N=N,
                          V=V,
                          L=L)

# save completed bundle, ready for gwas, to the appropriate folder
# we use stan json format as other formats (such as RDump) aren't supported by cmdstan

json_file <- tempfile(fileext = ".json", tmpdir = "~/Desktop/simulating_klebsiella_data/EF_100_ncausal_8/", pattern = "oligogenic_EF_100_H_0_9_ncausal_8_associative_")
cmdstanr::write_stan_json(data = cleaned_gwas_data, file = json_file)

# create index of variant positions
index <- as.matrix(colnames(genotype))
colnames(index) <- c("variant")
write.csv(x = index,
          file = "~/Desktop/simulating_klebsiella_data/EF_100_ncausal_8/oligogenic_EF_100_H_0_9_ncausal_8_associative_variant_index.csv")

############################################################################## OLIGO EF=100 N_causal_SNPs=20

# read in genotype (SNPs)
genotype <- read.delim("~/Desktop/simulating_klebsiella_data/pres_abs_matrix.txt")
rownames(genotype) <- genotype[,1]
genotype <- genotype[,-1]
genotype <- t(genotype)

# read in lineages
lineages <- read.delim("~/Desktop/simulating_klebsiella_data/PP_clusters.txt", sep = "\t", header = FALSE)

# read in phenotype
phenotype <- read.delim("~/Desktop/simulating_klebsiella_data/EF_100_ncausal_20/EF_100_ncausal_20.phen", sep = " ", header = FALSE)

# Reorder samples in lineages file to match genotype and phenotype
sample_index <- match(phenotype[,2], lineages[,1])
lineages <- lineages[sample_index,]

### Concat into one file (to check for correctness)

# define V = number of variants, L = number of lineages, N = number of samples
V <- as.numeric(ncol(genotype))
#L <- as.numeric(max(lineages[,2]))
L <- 347
N <- as.numeric(nrow(phenotype))

print(V)
print(L)
print(N)

checkme <- cbind(genotype[,1], phenotype, lineages)

# convert lineages into one-hot format
lineages <- model.matrix(lineages[,1]~as.factor(lineages[,2]), lineages)[,-1]

# make data into file
cleaned_gwas_data <- list(phenotype = phenotype[,3],
                          variant_matrix = genotype,
                          lineage_matrix = lineages,
                          N=N,
                          V=V,
                          L=L)

# save completed bundle, ready for gwas, to the appropriate folder
# we use stan json format as other formats (such as RDump) aren't supported by cmdstan

json_file <- tempfile(fileext = ".json", tmpdir = "~/Desktop/simulating_klebsiella_data/EF_100_ncausal_20/", pattern = "oligogenic_EF_100_H_0_9_ncausal_20_associative_")
cmdstanr::write_stan_json(data = cleaned_gwas_data, file = json_file)

# create index of variant positions
index <- as.matrix(colnames(genotype))
colnames(index) <- c("variant")
write.csv(x = index,
          file = "~/Desktop/simulating_klebsiella_data/EF_100_ncausal_20/oligogenic_EF_100_H_0_9_ncausal_20_associative_variant_index.csv")
