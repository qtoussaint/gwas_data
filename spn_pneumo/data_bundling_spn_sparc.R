################################################### LOAD LIBRARIES

library("vcfR")
library("dplyr")

################################################### BUNDLE DATA

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

### to retain all MIC categories:
# make equal values equal (eg 1.00 -> 1)
#penicillin_mics$MIC <- replace(penicillin_mics$MIC, penicillin_mics$MIC=="0.50", 0.5)
#penicillin_mics$MIC <- replace(penicillin_mics$MIC, penicillin_mics$MIC=="1.00", 1)
#penicillin_mics$MIC <- replace(penicillin_mics$MIC, penicillin_mics$MIC=="2.00", 2)
#penicillin_mics$MIC <- replace(penicillin_mics$MIC, penicillin_mics$MIC=="4.00", 4)
#penicillin_mics$MIC <- replace(penicillin_mics$MIC, penicillin_mics$MIC==">4.00", ">4")
#penicillin_mics$MIC <- replace(penicillin_mics$MIC, penicillin_mics$MIC=="<=0.03", "<0.03")
#penicillin_mics$MIC <- replace(penicillin_mics$MIC, penicillin_mics$MIC=="<.016", "<0.016")

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

# combine into one file (to check for correctness)
checkme <- cbind(genotype[,1], penicillin_mics, lineages)

# define V = number of variants, L = number of lineages, N = number of samples
V <- as.numeric(ncol(genotype))
L <- as.numeric(max(lineages$Cluster))
N <- as.numeric(nrow(penicillin_mics))
K <- as.numeric(length(unique(penicillin_mics$K)))

print(V)
print(L)
print(N)
print(K)

# convert lineages into one-hot format
lineages <- model.matrix(Taxon~as.factor(Cluster), lineages)[,-1]

# make data into file
cleaned_gwas_data <- list(phenotype = penicillin_mics$K,
                          variant_matrix = genotype,
                          lineage_matrix = lineages,
                          N=N,
                          V=V,
                          L=L,
                          K=K)

# save completed bundle, ready for gwas, to the appropriate folder
# we use stan json format as other formats (such as RDump) aren't supported by cmdstan
json_file <- tempfile(fileext = ".json", tmpdir = "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/gwas_bundles", pattern = "spn_sparc_gwas_bundle_K_13_maf_0.05_")
cmdstanr::write_stan_json(data = cleaned_gwas_data, file = json_file)

# create index of variant positions
index <- as.matrix(colnames(genotype))
colnames(index) <- c("variant")
write.csv(x = index,
          file = "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/gwas_bundles/spn_sparc_gwas_bundle_K_13_maf_0.05_associative_variant_index.csv")

################################################### CREATE PREDICTIVE BUNDLE

# set seed for reproducibility
set.seed(101010)

# Split data into training and testing data
N_train <- round(N*0.8)
N_test <- N-N_train

train_ind <- sample(c(1:N), size = N_train, replace = FALSE)

var_train <- genotype[train_ind,]
var_test <- genotype[-train_ind,]
lin_train <- lineages[train_ind,]
lin_test <- lineages[-train_ind,]
phen_train <- penicillin_mics$K[train_ind]
phen_test <- penicillin_mics$K[-train_ind]
phen_test_w_samples <- penicillin_mics[-train_ind,1:2]
colnames(phen_test_w_samples) <- c("sample_id", "true_phenotype")

# make data into file
cleaned_gwas_data <- list(training_phenotypes = phen_train,
                          training_variants = var_train,
                          training_lineages = lin_train,
                          test_variants = var_test,
                          test_lineages = lin_test,
                          N_train = N_train,
                          N_test = N_test,
                          V=V,
                          L=L,
                          K=K)

# save completed bundle, ready for gwas, to the appropriate folder
# we use stan json format as other formats (such as RDump) aren't supported by cmdstan
json_file <- tempfile(fileext = ".json", tmpdir = "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/gwas_bundles", pattern = "spn_sparc_gwas_bundle_K_13_maf_0.05_predictive_")
cmdstanr::write_stan_json(data = cleaned_gwas_data, file = json_file)

# write file with test phenotypes (ground truth)
write.csv(phen_test_w_samples, file = "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/gwas_bundles/spn_sparc_gwas_bundle_K_13_maf_0.05_predictive_test_phenotypes.csv")

################################################### SUBCLUSTERING VERSION

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

### to retain all MIC categories:
# make equal values equal (eg 1.00 -> 1)
#penicillin_mics$MIC <- replace(penicillin_mics$MIC, penicillin_mics$MIC=="0.50", 0.5)
#penicillin_mics$MIC <- replace(penicillin_mics$MIC, penicillin_mics$MIC=="1.00", 1)
#penicillin_mics$MIC <- replace(penicillin_mics$MIC, penicillin_mics$MIC=="2.00", 2)
#penicillin_mics$MIC <- replace(penicillin_mics$MIC, penicillin_mics$MIC=="4.00", 4)
#penicillin_mics$MIC <- replace(penicillin_mics$MIC, penicillin_mics$MIC==">4.00", ">4")
#penicillin_mics$MIC <- replace(penicillin_mics$MIC, penicillin_mics$MIC=="<=0.03", "<0.03")
#penicillin_mics$MIC <- replace(penicillin_mics$MIC, penicillin_mics$MIC=="<.016", "<0.016")

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

# combine into one file (to check for correctness)
checkme <- cbind(genotype[,1], penicillin_mics, lineages)

# define V = number of variants, L = number of lineages, N = number of samples
V <- as.numeric(ncol(genotype))
L <- as.numeric(max(lineages$Cluster))
S <- 
N <- as.numeric(nrow(penicillin_mics))
K <- as.numeric(length(unique(penicillin_mics$K)))

print(V)
print(L)
print(N)
print(K)

# convert lineages into one-hot format
lineages <- model.matrix(Taxon~as.factor(Cluster), lineages)[,-1]

# make data into file
cleaned_gwas_data <- list(phenotype = penicillin_mics$K,
                          variant_matrix = genotype,
                          lineage_matrix = lineages,
                          N=N,
                          V=V,
                          L=L,
                          K=K)

# save completed bundle, ready for gwas, to the appropriate folder
# we use stan json format as other formats (such as RDump) aren't supported by cmdstan
json_file <- tempfile(fileext = ".json", tmpdir = "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/gwas_bundles", pattern = "spn_sparc_gwas_bundle_K_13_maf_0.05_")
cmdstanr::write_stan_json(data = cleaned_gwas_data, file = json_file)

# create index of variant positions
index <- as.matrix(colnames(genotype))
colnames(index) <- c("variant")
write.csv(x = index,
          file = "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/gwas_bundles/spn_sparc_gwas_bundle_K_13_maf_0.05_associative_variant_index.csv")

################################################### CREATE PREDICTIVE BUNDLE

# set seed for reproducibility
set.seed(101010)

# Split data into training and testing data
N_train <- round(N*0.8)
N_test <- N-N_train

train_ind <- sample(c(1:N), size = N_train, replace = FALSE)

var_train <- genotype[train_ind,]
var_test <- genotype[-train_ind,]
lin_train <- lineages[train_ind,]
lin_test <- lineages[-train_ind,]
phen_train <- penicillin_mics$K[train_ind]
phen_test <- penicillin_mics$K[-train_ind]
phen_test_w_samples <- penicillin_mics[-train_ind,1:2]
colnames(phen_test_w_samples) <- c("sample_id", "true_phenotype")

# make data into file
cleaned_gwas_data <- list(training_phenotypes = phen_train,
                          training_variants = var_train,
                          training_lineages = lin_train,
                          test_variants = var_test,
                          test_lineages = lin_test,
                          N_train = N_train,
                          N_test = N_test,
                          V=V,
                          L=L,
                          K=K)

# save completed bundle, ready for gwas, to the appropriate folder
# we use stan json format as other formats (such as RDump) aren't supported by cmdstan
json_file <- tempfile(fileext = ".json", tmpdir = "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/gwas_bundles", pattern = "spn_sparc_gwas_bundle_K_13_maf_0.05_predictive_")
cmdstanr::write_stan_json(data = cleaned_gwas_data, file = json_file)

# write file with test phenotypes (ground truth)
write.csv(phen_test_w_samples, file = "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/gwas_bundles/spn_sparc_gwas_bundle_K_13_maf_0.05_predictive_test_phenotypes.csv")


################################################### LEAVE ONE LINEAGE OUT PREDICTION VERSION

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

### to retain all MIC categories:
# make equal values equal (eg 1.00 -> 1)
#penicillin_mics$MIC <- replace(penicillin_mics$MIC, penicillin_mics$MIC=="0.50", 0.5)
#penicillin_mics$MIC <- replace(penicillin_mics$MIC, penicillin_mics$MIC=="1.00", 1)
#penicillin_mics$MIC <- replace(penicillin_mics$MIC, penicillin_mics$MIC=="2.00", 2)
#penicillin_mics$MIC <- replace(penicillin_mics$MIC, penicillin_mics$MIC=="4.00", 4)
#penicillin_mics$MIC <- replace(penicillin_mics$MIC, penicillin_mics$MIC==">4.00", ">4")
#penicillin_mics$MIC <- replace(penicillin_mics$MIC, penicillin_mics$MIC=="<=0.03", "<0.03")
#penicillin_mics$MIC <- replace(penicillin_mics$MIC, penicillin_mics$MIC=="<.016", "<0.016")

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

# define V = number of variants, L = number of lineages, N = number of samples
V <- as.numeric(ncol(genotype))
L <- as.numeric(max(lineages$Cluster))
N <- as.numeric(nrow(penicillin_mics))
K <- as.numeric(length(unique(penicillin_mics$K)))

# Split data into training and testing data
# Pick a relatively big cluster with a decent number of samples (here 25 samples)
N_test <- as.numeric(length(which(lineages$Cluster==5)))
N_train <- N-N_test
train_ind <- which(lineages$Cluster!=5)

# convert lineages into one-hot format
lineages <- model.matrix(Taxon~as.factor(Cluster), lineages)[,-1]

var_train <- genotype[train_ind,]
var_test <- genotype[-train_ind,]
lin_train <- lineages[train_ind,]
lin_test <- lineages[-train_ind,]
phen_train <- penicillin_mics$K[train_ind]
phen_test <- penicillin_mics$K[-train_ind]
phen_test_w_samples <- penicillin_mics[-train_ind,1:2]
colnames(phen_test_w_samples) <- c("sample_id", "true_phenotype")

# make data into file
cleaned_gwas_data <- list(training_phenotypes = phen_train,
                          training_variants = var_train,
                          training_lineages = lin_train,
                          test_variants = var_test,
                          test_lineages = lin_test,
                          N_train = N_train,
                          N_test = N_test,
                          V=V,
                          L=L,
                          K=K)

# save completed bundle, ready for gwas, to the appropriate folder
# we use stan json format as other formats (such as RDump) aren't supported by cmdstan
json_file <- tempfile(fileext = ".json", tmpdir = "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/gwas_bundles", pattern = "spn_sparc_leaveoneout_K_13_maf_0.05_predictive_")
cmdstanr::write_stan_json(data = cleaned_gwas_data, file = json_file)

# write file with test phenotypes (ground truth)
write.csv(phen_test_w_samples, file = "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/gwas_bundles/spn_sparc_leaveoneout_K_13_maf_0.05_predictive_test_phenotypes.csv")
