################################################ LOAD LIBRARIES

library(dplyr)
library(data.table)

################################################ GENERATE DATA BUNDLE (binary)

# read in phenotype
# use newer (2023) table
phenotype <- read.csv(file = "/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/CRyPTIC_reuse_table_20231208.csv")

# pick phenotype
# here I'm using the binary rifabutin resistance (RFB) phenotype
phenotype <- phenotype[,c(1,15)]

# remove NAs
phenotype <- na.omit(phenotype)

# change S/R to binary 0/1
phenotype$RFB_BINARY_PHENOTYPE[phenotype$RFB_BINARY_PHENOTYPE == "S"] <- 0
phenotype$RFB_BINARY_PHENOTYPE[phenotype$RFB_BINARY_PHENOTYPE == "R"] <- 1
phenotype$RFB_BINARY_PHENOTYPE <- lapply(phenotype$RFB_BINARY_PHENOTYPE, as.numeric)

# read in variants
genotype <- fread("/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/genotype/tb_mhunt_genotypes_cleaned.txt",
                       sep = " ",
                       header=F)
genotype <- as.data.frame(genotype)

# transpose genotype (samples in rows, variants in columns)
genotype <- transpose(genotype)

# make variant positions into header
colnames(genotype) <- as.character(genotype[1,])
genotype <- genotype[-1,]

# clean genotype sample names to just be sample (for matching later on)
genotype[,1] <- sub("\\[.*\\]E", "E", genotype[,1])
genotype[,1] <- sub(":GT", "", genotype[,1])
row.names(genotype) <- genotype[,1]
genotype <- genotype[,-1]

# change to haploid format
# convert missing to reference
genotype[] <- lapply(genotype, function(x) {
  x[x == "0/0"] <- 0
  x[x == "./."] <- 0
  x[x == "1/1"] <- 1
  x[x == "0"] <- 0
  x[x == "1"] <- 1
  x[x == "2/2"] <- 1
  return(x)
})

# convert genotype matrix to numeric
genotype <- genotype %>%
  mutate(across(where(is.character), as.numeric))

# read in lineages
lineages <- read.delim(file = "/nfs/research/jlees/jacqueline/gwas_data/code/gwas_data/tuberculosis/old/output_fastlin.txt")

# remove samples with k-mers that have more than one lineage designation (about 300 samples, probably due to contamination)
# clean up lineage designations by removing the median k-mer occurrences in parentheses
# remove lineage designations of empty strings ("")
lineages <- lineages %>% 
  filter(!grepl(",", lineages)) %>%
  mutate(lineages = sub("\\s+\\(.*\\)$", "", lineages)) %>%
  filter(lineages != "")

# Reorder samples in lineages file to match genotype and phenotype
sample_index <- match(phenotype$ENA_RUN, lineages$X.sample)
sample_index <- na.omit(sample_index)
lineages <- lineages[sample_index,]

sample_index <- match(lineages[,1], row.names(genotype))
sample_index <- na.omit(sample_index)
genotype <- genotype[sample_index,]

sample_index <- match(lineages[,1], phenotype$ENA_RUN)
sample_index <- na.omit(sample_index)
phenotype <- phenotype[sample_index,]
phenotype$RFB_BINARY_PHENOTYPE <- unlist(phenotype$RFB_BINARY_PHENOTYPE)

# define V = number of variants, L = number of lineages, N = number of samples
V <- as.numeric(ncol(genotype))
L <- as.numeric(length(unique(lineages$lineages)))
N <- as.numeric(nrow(phenotype))

print(V)
print(L)
print(N)

# convert lineages into one-hot format
lineages <- lineages[,c(1,5)]
lineages <- model.matrix(X.sample~as.factor(lineages), lineages)[,-1]

# make data into file
cleaned_gwas_data <- list(phenotype = phenotype$RFB_BINARY_PHENOTYPE,
                          variant_matrix = genotype,
                          lineage_matrix = lineages,
                          N=N,
                          V=V,
                          L=L)

# save completed bundle, ready for gwas, to the appropriate folder
# we use stan json format as other formats (such as RDump) aren't supported by cmdstan

json_file <- tempfile(fileext = ".json", tmpdir = "/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/gwas_bundles/", pattern = "tb_mhunt_rfb_binary_")
cmdstanr::write_stan_json(data = cleaned_gwas_data, file = json_file)

# create index of variant positions
index <- as.matrix(colnames(genotype))
colnames(index) <- c("variant")
write.csv(x = index,
          file = "/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/gwas_bundles/tb_mhunt_rfb_binary_variant_index.csv")

################################################ GENERATE DATA BUNDLE (binary predictive)

# set seed for reproducibility
set.seed(42)

# Split data into training and testing data
N_train <- round(N*0.8)
N_test <- N-N_train

set.seed(123)
train_ind <- sample(c(1:N), size = N_train, replace = FALSE)

var_train <- genotype[train_ind,]
var_test <- genotype[-train_ind,]
lin_train <- lineages[train_ind,]
lin_test <- lineages[-train_ind,]
phen_train <- phenotype$RFB_BINARY_PHENOTYPE[train_ind]
phen_test <- phenotype$RFB_BINARY_PHENOTYPE[-train_ind]
phen_test_w_samples <- phenotype[-train_ind,]
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
json_file <- tempfile(fileext = ".json", tmpdir = "/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/gwas_bundles/", pattern = "tb_mhunt_rfb_binary_predictive_")
cmdstanr::write_stan_json(data = cleaned_gwas_data, file = json_file)

# write file with test phenotypes (ground truth)
write.csv(phen_test_w_samples, file = "/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/gwas_bundles/tb_mhunt_rfb_binary_predictive_test_phenotypes.csv")

################################################ GENERATE DATA BUNDLE (MICs)

# read in phenotype
# use newer (2023) table
phenotype <- read.csv(file = "/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/CRyPTIC_reuse_table_20231208.csv")

# pick phenotype
# here I'm using the rifabutin resistance (RFB) phenotype with MICs
phenotype <- phenotype[,c(1,28)]

# remove NAs
phenotype <- na.omit(phenotype)

# check unique values and put in order
order <- unique(phenotype$RFB_MIC)
desired_order <- order[c(1,3,6,4,7,5,2)]

# change MICs to ordered categories in new row
for (i in 1:nrow(phenotype)) {
  index <- match(phenotype$RFB_MIC[i], desired_order)
  phenotype$RFB_MIC_categorical[i] <- index
}

# read in variants
genotype <- fread("/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/genotype/tb_mhunt_genotypes_cleaned.txt",
                  sep = " ",
                  header=F)
genotype <- as.data.frame(genotype)

# transpose genotype (samples in rows, variants in columns)
genotype <- transpose(genotype)

# make variant positions into header
colnames(genotype) <- as.character(genotype[1,])
genotype <- genotype[-1,]

# clean genotype sample names to just be sample (for matching later on)
genotype[,1] <- sub("\\[.*\\]E", "E", genotype[,1])
genotype[,1] <- sub(":GT", "", genotype[,1])
row.names(genotype) <- genotype[,1]
genotype <- genotype[,-1]

# change to haploid format
# convert missing to reference
genotype[] <- lapply(genotype, function(x) {
  x[x == "0/0"] <- 0
  x[x == "./."] <- 0
  x[x == "1/1"] <- 1
  x[x == "0"] <- 0
  x[x == "1"] <- 1
  x[x == "2/2"] <- 1
  return(x)
})

# convert genotype matrix to numeric
genotype <- genotype %>%
  mutate(across(where(is.character), as.numeric))

# read in lineages
lineages <- read.delim(file = "/nfs/research/jlees/jacqueline/gwas_data/code/gwas_data/tuberculosis/old/output_fastlin.txt")

# remove samples with k-mers that have more than one lineage designation (about 300 samples, probably due to contamination)
# clean up lineage designations by removing the median k-mer occurrences in parentheses
# remove lineage designations of empty strings ("")
lineages <- lineages %>% 
  filter(!grepl(",", lineages)) %>%
  mutate(lineages = sub("\\s+\\(.*\\)$", "", lineages)) %>%
  filter(lineages != "")

# Reorder samples in lineages file to match genotype and phenotype
sample_index <- match(phenotype$ENA_RUN, lineages$X.sample)
sample_index <- na.omit(sample_index)
lineages <- lineages[sample_index,]

sample_index <- match(lineages[,1], row.names(genotype))
sample_index <- na.omit(sample_index)
genotype <- genotype[sample_index,]

sample_index <- match(lineages[,1], phenotype$ENA_RUN)
sample_index <- na.omit(sample_index)
phenotype <- phenotype[sample_index,]
phenotype$RFB_MIC_categorical <- unlist(phenotype$RFB_MIC_categorical)

# define V = number of variants, L = number of lineages, N = number of samples
V <- as.numeric(ncol(genotype))
L <- as.numeric(length(unique(lineages$lineages)))
N <- as.numeric(nrow(phenotype))
K <- as.numeric(length(unique(phenotype$RFB_MIC_categorical)))

# convert lineages into one-hot format
lineages <- lineages[,c(1,5)]
lineages <- model.matrix(X.sample~as.factor(lineages), lineages)[,-1]

# make data into file
cleaned_gwas_data <- list(phenotype = phenotype$RFB_MIC_categorical,
                          variant_matrix = genotype,
                          lineage_matrix = lineages,
                          N=N,
                          V=V,
                          L=L,
                          K=K)

# save completed bundle, ready for gwas, to the appropriate folder
# we use stan json format as other formats (such as RDump) aren't supported by cmdstan

json_file <- tempfile(fileext = ".json", tmpdir = "/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/gwas_bundles/", pattern = "tb_mhunt_rfb_MICs_")
cmdstanr::write_stan_json(data = cleaned_gwas_data, file = json_file)

# create index of variant positions
index <- as.matrix(colnames(genotype))
colnames(index) <- c("variant")
write.csv(x = index,
          file = "/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/gwas_bundles/tb_mhunt_rfb_MICs_variant_index.csv")

################################################ GENERATE DATA BUNDLE (MICs predictive)

# set seed for reproducibility
set.seed(42)

# Split data into training and testing data
N_train <- round(N*0.8)
N_test <- N-N_train

set.seed(123)
train_ind <- sample(c(1:N), size = N_train, replace = FALSE)

var_train <- genotype[train_ind,]
var_test <- genotype[-train_ind,]
lin_train <- lineages[train_ind,]
lin_test <- lineages[-train_ind,]
phen_train <- phenotype$RFB_MIC_categorical[train_ind]
phen_test <- phenotype$RFB_MIC_categorical[-train_ind]
phen_test_w_samples <- phenotype[-train_ind,]
colnames(phen_test_w_samples) <- c("sample_id", "true_phenotype (MIC)", "true phenotype (categorical MIC)")


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
json_file <- tempfile(fileext = ".json", tmpdir = "/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/gwas_bundles/", pattern = "tb_mhunt_rfb_MICs_predictive_")
cmdstanr::write_stan_json(data = cleaned_gwas_data, file = json_file)

# write file with test phenotypes (ground truth)
write.csv(phen_test_w_samples, file = "/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/gwas_bundles/tb_mhunt_rfb_MICs_predictive_test_phenotypes.csv")


################################################ GENERATE LEAVE ONE OUT BUNDLE WITH MICs

# read in phenotype
# use newer (2023) table
phenotype <- read.csv(file = "/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/CRyPTIC_reuse_table_20231208.csv")

# pick phenotype
# here I'm using the rifabutin resistance (RFB) phenotype with MICs
phenotype <- phenotype[,c(1,28)]

# remove NAs
phenotype <- na.omit(phenotype)

# check unique values and put in order
order <- unique(phenotype$RFB_MIC)
desired_order <- order[c(1,3,6,4,7,5,2)]

# change MICs to ordered categories in new row
for (i in 1:nrow(phenotype)) {
  index <- match(phenotype$RFB_MIC[i], desired_order)
  phenotype$RFB_MIC_categorical[i] <- index
}

# read in variants
genotype <- fread("/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/genotype/tb_mhunt_genotypes_cleaned.txt",
                  sep = " ",
                  header=F)
genotype <- as.data.frame(genotype)

# transpose genotype (samples in rows, variants in columns)
genotype <- transpose(genotype)

# make variant positions into header
colnames(genotype) <- as.character(genotype[1,])
genotype <- genotype[-1,]

# clean genotype sample names to just be sample (for matching later on)
genotype[,1] <- sub("\\[.*\\]E", "E", genotype[,1])
genotype[,1] <- sub(":GT", "", genotype[,1])
row.names(genotype) <- genotype[,1]
genotype <- genotype[,-1]

# change to haploid format
# convert missing to reference
genotype[] <- lapply(genotype, function(x) {
  x[x == "0/0"] <- 0
  x[x == "./."] <- 0
  x[x == "1/1"] <- 1
  x[x == "0"] <- 0
  x[x == "1"] <- 1
  x[x == "2/2"] <- 1
  return(x)
})

# convert genotype matrix to numeric
genotype <- genotype %>%
  mutate(across(where(is.character), as.numeric))

# read in lineages
lineages <- read.delim(file = "/nfs/research/jlees/jacqueline/gwas_data/code/gwas_data/tuberculosis/old/output_fastlin.txt")

# remove samples with k-mers that have more than one lineage designation (about 300 samples, probably due to contamination)
# clean up lineage designations by removing the median k-mer occurrences in parentheses
# remove lineage designations of empty strings ("")
lineages <- lineages %>% 
  filter(!grepl(",", lineages)) %>%
  mutate(lineages = sub("\\s+\\(.*\\)$", "", lineages)) %>%
  filter(lineages != "")

# Reorder samples in lineages file to match genotype and phenotype
sample_index <- match(phenotype$ENA_RUN, lineages$X.sample)
sample_index <- na.omit(sample_index)
lineages <- lineages[sample_index,]

sample_index <- match(lineages[,1], row.names(genotype))
sample_index <- na.omit(sample_index)
genotype <- genotype[sample_index,]

sample_index <- match(lineages[,1], phenotype$ENA_RUN)
sample_index <- na.omit(sample_index)
phenotype <- phenotype[sample_index,]
phenotype$RFB_MIC_categorical <- unlist(phenotype$RFB_MIC_categorical)

# define V = number of variants, L = number of lineages, N = number of samples
V <- as.numeric(ncol(genotype))
L <- as.numeric(length(unique(lineages$lineages)))
N <- as.numeric(nrow(phenotype))
K <- as.numeric(length(unique(phenotype$RFB_MIC_categorical)))

# Split data into training and testing data
N_test <- as.numeric(length(which(lineages$lineages!="4.5")))
N_train <- N-N_test
train_ind <- as.numeric(which(lineages$lineages=="4.5"))

# convert lineages into one-hot format
lineages <- lineages[,c(1,5)]
lineages <- model.matrix(X.sample~as.factor(lineages), lineages)[,-1]

var_train <- genotype[train_ind,]
var_test <- genotype[-train_ind,]
lin_train <- lineages[train_ind,]
lin_test <- lineages[-train_ind,]
phen_train <- phenotype$RFB_MIC_categorical[train_ind]
phen_test <- phenotype$RFB_MIC_categorical[-train_ind]
phen_test_w_samples <- phenotype[-train_ind,]
colnames(phen_test_w_samples) <- c("sample_id", "true_phenotype (MIC)", "true phenotype (categorical MIC)")

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
json_file <- tempfile(fileext = ".json", tmpdir = "/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/gwas_bundles/", pattern = "tb_mhunt_rfb_MICs_leaveoneout_")
cmdstanr::write_stan_json(data = cleaned_gwas_data, file = json_file)

# write file with test phenotypes (ground truth)
write.csv(phen_test_w_samples, file = "/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/gwas_bundles/tb_mhunt_rfb_MICs_leaveoneout_test_phenotypes.csv")
