################################################### LOAD LIBRARIES

library("vcfR")
library("dplyr")
library("cmdstanr")

################################################### BUNDLE DATA K=13

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

################################################### SUBCLUSTERING VERSION (remove samples in subclusters with =<3 samples)

# read in genotype
sparc_snps <- read.vcfR("/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype/sparc_snps.vcf.gz", verbose = TRUE)
genotype <- extract.gt(sparc_snps, as.numeric = TRUE, IDtoRowNames = T)

# read in lineages and sublineages from poppipe
lineages <- read.csv(file = "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/lineages/sublineages_from_poppipe/min_cluster_3/all_clusters.txt")
sublineages <- lineages[,c(1,3)]

# remove NA samples
NA_sublin <- which(is.na(sublineages[,2]))
sublineages <- sublineages[-NA_sublin,]
lineages <- lineages[-NA_sublin,c(1,2)]

# read in phenotype data
phenotypes <- read.csv("/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/phenotype/sparc_metadata.csv")

# select desired phenotype (here it's benzylpenicillin MICs)
penicillin_mics <- as.data.frame(phenotypes$Benzylpenicillin.MIC..ug.mL.)
penicillin_mics <- data.frame(MIC=penicillin_mics[,1], ID=phenotypes$Lane.name)

# ensure order/number of samples matches between these 4 data sets
# sort/filter those that don't
index <- match(lineages$Taxon, dimnames(genotype)[[2]])
index <- index[-which(is.na(index))]
genotype <- genotype[,index]

index <- match(dimnames(genotype)[[2]], penicillin_mics$ID)
penicillin_mics <- penicillin_mics[index,]

index <- match(penicillin_mics$ID, lineages$Taxon)
lineages <- lineages[index,]

index <- match(penicillin_mics$ID, sublineages$Taxon)
sublineages <- sublineages[index,]

# transpose genotype data (need samples in rows and variants in columns)
genotype <- t(genotype)

# find and replace missing data in genotype file
# will use change all missing SNPs to 0
genotype[is.na(genotype)] <- 0

# filter for MAF>0.05 (if desired)
small_maf <- which((colSums(genotype)/nrow(genotype))<0.05)
genotype <- genotype[,-small_maf]

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

# define V = number of variants, L = number of lineages, S = number of sublineages, N = number of samples
V <- as.numeric(ncol(genotype))
L <- as.numeric(length(unique(lineages$Strain)))
S <- as.numeric(length(unique(sublineages$Subcluster_1)))
N <- as.numeric(nrow(penicillin_mics))
K <- as.numeric(length(unique(penicillin_mics$K)))

# convert coarse lineages into one-hot format
lineages <- model.matrix(Taxon~as.factor(Strain), lineages)[,-1]

# convert fine lineages into one-hot format
sublineages <- model.matrix(Taxon~as.factor(Subcluster_1), sublineages)[,-1]

# make data into file
cleaned_gwas_data <- list(phenotype = penicillin_mics$K,
                          variant_matrix = genotype,
                          lineage_matrix = lineages,
                          sublineage_matrix = sublineages,
                          N=N,
                          V=V,
                          L=L,
                          S=S,
                          K=K)


# save completed bundle, ready for gwas, to the appropriate folder
# we use stan json format as other formats (such as RDump) aren't supported by cmdstan
json_file <- tempfile(fileext = ".json", tmpdir = "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/gwas_bundles", pattern = "spn_sparc_gwas_bundle_K_13_maf_0.05_sublineages_remove_")
cmdstanr::write_stan_json(data = cleaned_gwas_data, file = json_file)

# create index of variant positions
index <- as.matrix(colnames(genotype))
colnames(index) <- c("variant")
write.csv(x = index,
          file = "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/gwas_bundles/spn_sparc_gwas_bundle_K_13_maf_0.05_sublineages_remove_variant_index.csv")

################################################### SUBCLUSTERING VERSION (replace subclusters with =<3 samples with coarse clusters)

# read in genotype
sparc_snps <- read.vcfR("/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype/sparc_snps.vcf.gz", verbose = TRUE)
genotype <- extract.gt(sparc_snps, as.numeric = TRUE, IDtoRowNames = T)

# read in lineages and sublineages from poppipe
lineages <- read.csv(file = "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/lineages/sublineages_from_poppipe/min_cluster_3/all_clusters.txt")
sublineages <- lineages[,c(1,3)]
lineages <- lineages[,c(1,2)]

# replace NAs with coarse cluster designation
# appended zeros to avoid accidentally assigning a to a subcluster that already exists
NA_sublin <- which(is.na(sublineages[,2]))
replacement <- as.numeric(paste(lineages[NA_sublin,2], "0000", sep = ""))
sublineages[NA_sublin,2] <- replacement

# read in phenotype data
phenotypes <- read.csv("/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/phenotype/sparc_metadata.csv")

# select desired phenotype (here it's benzylpenicillin MICs)
penicillin_mics <- as.data.frame(phenotypes$Benzylpenicillin.MIC..ug.mL.)
penicillin_mics <- data.frame(MIC=penicillin_mics[,1], ID=phenotypes$Lane.name)

# ensure order/number of samples matches between these 4 data sets
# sort/filter those that don't
index <- match(lineages$Taxon, dimnames(genotype)[[2]])
index <- index[-which(is.na(index))]
genotype <- genotype[,index]

index <- match(dimnames(genotype)[[2]], penicillin_mics$ID)
penicillin_mics <- penicillin_mics[index,]

index <- match(penicillin_mics$ID, lineages$Taxon)
lineages <- lineages[index,]

index <- match(penicillin_mics$ID, sublineages$Taxon)
sublineages <- sublineages[index,]

# transpose genotype data (need samples in rows and variants in columns)
genotype <- t(genotype)

# find and replace missing data in genotype file
# will use change all missing SNPs to 0
genotype[is.na(genotype)] <- 0

# filter for MAF>0.05 (if desired)
small_maf <- which((colSums(genotype)/nrow(genotype))<0.05)
genotype <- genotype[,-small_maf]

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

# define V = number of variants, L = number of lineages, S = number of sublineages, N = number of samples
V <- as.numeric(ncol(genotype))
L <- as.numeric(length(unique(lineages$Strain)))
S <- as.numeric(length(unique(sublineages$Subcluster_1)))
N <- as.numeric(nrow(penicillin_mics))
K <- as.numeric(length(unique(penicillin_mics$K)))

# convert coarse lineages into one-hot format
lineages <- model.matrix(Taxon~as.factor(Strain), lineages)[,-1]

# convert fine lineages into one-hot format
sublineages <- model.matrix(Taxon~as.factor(Subcluster_1), sublineages)[,-1]

# make data into file
cleaned_gwas_data <- list(phenotype = penicillin_mics$K,
                          variant_matrix = genotype,
                          lineage_matrix = lineages,
                          sublineage_matrix = sublineages,
                          N=N,
                          V=V,
                          L=L,
                          S=S,
                          K=K)

# save completed bundle, ready for gwas, to the appropriate folder
# we use stan json format as other formats (such as RDump) aren't supported by cmdstan
json_file <- tempfile(fileext = ".json", tmpdir = "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/gwas_bundles", pattern = "spn_sparc_gwas_bundle_K_13_maf_0.05_sublineages_replace_")
cmdstanr::write_stan_json(data = cleaned_gwas_data, file = json_file)

# create index of variant positions
index <- as.matrix(colnames(genotype))
colnames(index) <- c("variant")
write.csv(x = index,
          file = "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/gwas_bundles/spn_sparc_gwas_bundle_K_13_maf_0.05_sublineages_replace_variant_index.csv")

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


##################################### 5 MIC CATEGORIES

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
# 4 bins (from 26 total)
MICs <- sort(unique(penicillin_mics$MIC))
print(unique(cut(MICs, breaks=c(0, 0.032, 0.065, 0.2, 2, 6))))

# cut according to give breakpoints and add as separate column
penicillin_complete <- penicillin_mics %>% mutate(MIC_factor = cut(penicillin_mics$MIC, breaks=c(0, 0.032, 0.065, 0.2, 2, 6)))

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

ggplot(data = penicillin_mics, aes(x = K)) +
  geom_bar()

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
json_file <- tempfile(fileext = ".json", tmpdir = "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/gwas_bundles", pattern = "spn_sparc_gwas_bundle_K_5_maf_0.05_")
cmdstanr::write_stan_json(data = cleaned_gwas_data, file = json_file)

# create index of variant positions
index <- as.matrix(colnames(genotype))
colnames(index) <- c("variant")
write.csv(x = index,
          file = "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/gwas_bundles/spn_sparc_gwas_bundle_K_5_maf_0.05_associative_variant_index.csv")

##################################### BINARY

# read in genotype 
sparc_snps <- read.vcfR("/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype/sparc_snps.vcf.gz", verbose = TRUE)
genotype <- extract.gt(sparc_snps, as.numeric = TRUE, IDtoRowNames = T)

# read in lineages
lineages <- read.csv("/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/lineages/database_clusters.csv")

# read in phenotype data
phenotypes <- read.delim("/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/pyseer_data/pyseer_tutorial/resistances.pheno")
colnames(phenotypes) <- c("samples", "binary")

# ensure order/number of samples matches between these 3 data sets
# sort/filter those that dont
index <- match(dimnames(genotype)[[2]], phenotypes$samples)
phenotypes <- phenotypes[index,]

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

# make numeric
phenotypes$binary <- as.numeric(phenotypes$binary)

# switch order of columns and rename
phenotypes <- data.frame(ID=phenotypes$sample, binary=phenotypes$binary)

# combine into one file (to check for correctness)
checkme <- cbind(genotype[,1], phenotypes, lineages)

# define V = number of variants, L = number of lineages, N = number of samples
V <- as.numeric(ncol(genotype))
L <- as.numeric(max(lineages$Cluster))
N <- as.numeric(nrow(phenotypes))

print(V)
print(L)
print(N)

# convert lineages into one-hot format
lineages <- model.matrix(Taxon~as.factor(Cluster), lineages)[,-1]

# make data into file
cleaned_gwas_data <- list(phenotype = phenotypes$binary,
                          variant_matrix = genotype,
                          lineage_matrix = lineages,
                          N=N,
                          V=V,
                          L=L)

# save completed bundle, ready for gwas, to the appropriate folder
# we use stan json format as other formats (such as RDump) aren't supported by cmdstan
json_file <- tempfile(fileext = ".json", tmpdir = "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/gwas_bundles/penicillin", pattern = "spn_sparc_binary_penicillin_maf_0.05_")
cmdstanr::write_stan_json(data = cleaned_gwas_data, file = json_file)

# create index of variant positions
index <- as.numeric(as.matrix(colnames(genotype)))
variant_pos <- as.numeric(sparc_snps@fix[,"POS"])[index] # REAL GENOME COORDINATES
colnames(variant_pos) <- c("variant")
write.csv(x = variant_pos,
          file = "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/gwas_bundles/penicillin/spn_sparc_binary_penicillin_maf_0.05_associative_variant_index.csv")

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
phen_train <- phenotypes$binary[train_ind]
phen_test <- phenotypes$binary[-train_ind]
phen_test_w_samples <- phenotypes[-train_ind,1:2]
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
                          L=L)

# save completed bundle, ready for gwas, to the appropriate folder
# we use stan json format as other formats (such as RDump) aren't supported by cmdstan
json_file <- tempfile(fileext = ".json", tmpdir = "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/gwas_bundles", pattern = "spn_sparc_gwas_binary_maf_0.05_predictive_")
cmdstanr::write_stan_json(data = cleaned_gwas_data, file = json_file)

# write file with test phenotypes (ground truth)
write.csv(phen_test_w_samples, file = "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/gwas_bundles/spn_sparc_gwas_binary_maf_0.05_predictive_test_phenotypes.csv")


##################################### BINARY TRIMETHOPRIM

# read in genotype 
sparc_snps <- read.vcfR("/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype/sparc_snps.vcf.gz", verbose = TRUE)
genotype <- extract.gt(sparc_snps, as.numeric = TRUE, IDtoRowNames = T)
variant_pos <- as.numeric(sparc_snps@fix[,"POS"])

# read in lineages
lineages <- read.csv("/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/lineages/database_clusters.csv")

# read in phenotype data
phenotypes <- read.delim("/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/phenotype/sparc_metadata.csv", sep = ",")
phenotypes <- phenotypes[,c("Lane.name", "Trimethoprim.MIC..ug.mL.")]
colnames(phenotypes) <- c("samples", "MIC")

print(unique(phenotypes$MIC))

### from EUCAST 2025 guidelines for trimethoprim in S. pnuemoniae (R is >=2, S is <=1, units mg/L)

# clean duplicate values and remove symbols
phenotypes$MIC <- replace(phenotypes$MIC, phenotypes$MIC==">2.00", 2)
phenotypes$MIC <- replace(phenotypes$MIC, phenotypes$MIC=="<=0.25", 0.25)
phenotypes$MIC <- replace(phenotypes$MIC, phenotypes$MIC==">32", 32)
### "/" format does not indicate a fraction, but rather two antibiotics given to a patient at once
# the first number is the MIC for the reported antibiotic (here, trimethoprim)
phenotypes$MIC <- replace(phenotypes$MIC, phenotypes$MIC=="4/76", 4)
phenotypes$MIC <- replace(phenotypes$MIC, phenotypes$MIC==">4/76", 4)
phenotypes$MIC <- replace(phenotypes$MIC, phenotypes$MIC=="0.25/4", 0.25)
phenotypes$MIC <- replace(phenotypes$MIC, phenotypes$MIC=="0.5/9.", 0.5)
phenotypes$MIC <- replace(phenotypes$MIC, phenotypes$MIC=="0.12/2", 0.12)
phenotypes$MIC <- replace(phenotypes$MIC, phenotypes$MIC=="1/19", 1)
phenotypes$MIC <- replace(phenotypes$MIC, phenotypes$MIC=="<0.06/", 0.06)
phenotypes$MIC <- replace(phenotypes$MIC, phenotypes$MIC=="2/38", 2)

sort(unique(phenotypes$MIC))

phenotypes$MIC <- as.numeric(phenotypes$MIC)

# make binary based on EUCAST
phenotypes$MIC[which(phenotypes$MIC <= 1)] <- 0
phenotypes$MIC[which(phenotypes$MIC >= 2)] <- 1
phenotypes <- phenotypes[-which(phenotypes$MIC %in% c(1,0) == FALSE),]
colnames(phenotypes) <- c("samples", "binary")

# ensure order/number of samples matches between these 3 data sets
# sort/filter those that dont

index <- match(dimnames(genotype)[[2]], phenotypes$samples)
phenotypes <- phenotypes[index,]

phenotypes <- phenotypes[-which(is.na(phenotypes$binary)),]

index <- match(phenotypes$samples, dimnames(genotype)[[2]])
genotype <- genotype[,index]

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

# make numeric
phenotypes$binary <- as.numeric(phenotypes$binary)

# switch order of columns and rename
phenotypes <- data.frame(ID=phenotypes$sample, binary=phenotypes$binary)

# combine into one file (to check for correctness)
checkme <- cbind(genotype[,1], phenotypes, lineages)

# define V = number of variants, L = number of lineages, N = number of samples
V <- as.numeric(ncol(genotype))
L <- as.numeric(max(lineages$Cluster))
N <- as.numeric(nrow(phenotypes))

print(V)
print(L)
print(N)

# convert lineages into one-hot format
lineages <- model.matrix(Taxon~as.factor(Cluster), lineages)[,-1]

# make data into file
cleaned_gwas_data <- list(phenotype = phenotypes$binary,
                          variant_matrix = genotype,
                          lineage_matrix = lineages,
                          N=N,
                          V=V,
                          L=L)

# save completed bundle, ready for gwas, to the appropriate folder
# we use stan json format as other formats (such as RDump) aren't supported by cmdstan
json_file <- tempfile(fileext = ".json", tmpdir = "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/gwas_bundles/trimethoprim/", pattern = "spn_sparc_binary_trimethoprim_maf_0.05_")
cmdstanr::write_stan_json(data = cleaned_gwas_data, file = json_file)

# create index of variant positions
index <- as.numeric(as.matrix(colnames(genotype)))
variant_pos <- as.numeric(sparc_snps@fix[,"POS"])[index] # REAL GENOME COORDINATES
colnames(variant_pos) <- c("variant")
write.csv(x = variant_pos,
          file = "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/gwas_bundles/trimethoprim/spn_sparc_binary_trimethoprim_maf_0.05_associative_variant_index.csv")

##################################### CONTINUOUS

# read in genotype 
sparc_snps <- read.vcfR("/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype/sparc_snps.vcf.gz", verbose = TRUE)
genotype <- extract.gt(sparc_snps, as.numeric = TRUE, IDtoRowNames = T)

# read in lineages
lineages <- read.csv("/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/lineages/database_clusters.csv")

# read in phenotype data
# same as benzylpenicillin MICs from sparc metadata but with >\<\= removed (and slightly dif numbers); from pyseer input data
phenotypes <- read.delim("/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/pyseer_data/phenotypes_pen_mic.txt")

# ensure order/number of samples matches between these 3 data sets
# sort/filter those that dont
index <- match(dimnames(genotype)[[2]], phenotypes$sample)
phenotypes <- phenotypes[index,]

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

# make numeric
phenotypes$continuous <- as.numeric(phenotypes$continuous)

# switch order of columns and rename
phenotypes <- data.frame(ID=phenotypes$sample, continuous=phenotypes$continuous)

# combine into one file (to check for correctness)
checkme <- cbind(genotype[,1], phenotypes, lineages)

# define V = number of variants, L = number of lineages, N = number of samples
V <- as.numeric(ncol(genotype))
L <- as.numeric(max(lineages$Cluster))
N <- as.numeric(nrow(phenotypes))

print(V)
print(L)
print(N)

# convert lineages into one-hot format
lineages <- model.matrix(Taxon~as.factor(Cluster), lineages)[,-1]

# make data into file
cleaned_gwas_data <- list(phenotype = phenotypes$continuous,
                          variant_matrix = genotype,
                          lineage_matrix = lineages,
                          N=N,
                          V=V,
                          L=L)

# save completed bundle, ready for gwas, to the appropriate folder
# we use stan json format as other formats (such as RDump) aren't supported by cmdstan
json_file <- tempfile(fileext = ".json", tmpdir = "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/gwas_bundles", pattern = "spn_sparc_gwas_bundle_continuous_maf_0.05_")
cmdstanr::write_stan_json(data = cleaned_gwas_data, file = json_file)

# create index of variant positions
index <- as.matrix(colnames(genotype))
colnames(index) <- c("variant")
write.csv(x = index,
          file = "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/gwas_bundles/spn_sparc_gwas_bundle_continuous_maf_0.05_associative_variant_index.csv")

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
phen_train <- phenotypes$continuous[train_ind]
phen_test <- phenotypes$continuous[-train_ind]
phen_test_w_samples <- phenotypes[-train_ind,1:2]
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
                          L=L)

# save completed bundle, ready for gwas, to the appropriate folder
# we use stan json format as other formats (such as RDump) aren't supported by cmdstan
json_file <- tempfile(fileext = ".json", tmpdir = "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/gwas_bundles", pattern = "spn_sparc_gwas_continuous_maf_0.05_predictive_")
cmdstanr::write_stan_json(data = cleaned_gwas_data, file = json_file)

# write file with test phenotypes (ground truth)
write.csv(phen_test_w_samples, file = "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/gwas_bundles/spn_sparc_gwas_continuous_maf_0.05_predictive_test_phenotypes.csv")
