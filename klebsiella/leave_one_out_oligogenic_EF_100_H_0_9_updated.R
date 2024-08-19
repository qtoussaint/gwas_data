### Pipeline for continuous GWAS model

# load libraries
library(cmdstanr)

# read in genotype
genotype <- read.delim("/nfs/research/jlees/jacqueline/gwas_data/klebsiella/Kpneu_MIC_prediction/data/genes/matrix_and_lineages/genes_presence_absence.txt")
rownames(genotype) <- genotype[,1]
genotype <- genotype[,-1]
genotype <- t(genotype)

# read in lineages
lineages <- read.delim("/nfs/research/jlees/jacqueline/gwas_data/klebsiella/lineages/klebsiella_poppunk_clusters.txt", sep = "\t", header = FALSE)

# read in phenotype
phenotype <- read.delim("/nfs/research/jlees/jacqueline/gwas_data/klebsiella/phenotypes/simul_oligogenic_EF_100.0_herit_0.9_traitprev_0.1.phen", sep = " ", header = FALSE)

# Reorder samples in lineages file to match genotype and phenotype
sample_index <- match(phenotype[,2], lineages[,1])
lineages <- lineages[sample_index,]

### Concat into one file (to check for correctness)

# define V = number of variants, L = number of lineages, N = number of samples
V <- as.numeric(ncol(genotype))
#L <- as.numeric(max(lineages[,2]))
L <- 347
N <- as.numeric(nrow(phenotype))

# Split data into training and testing data
# Pick a relatively big cluster with a decent number of samples (here 25 samples)
N_test <- as.numeric(length(which(lineages$V2=="PP3")))
N_train <- N-N_test
train_ind <- which(lineages$V2!="PP3")

# convert lineages into one-hot format
lineages <- model.matrix(lineages[,1]~as.factor(lineages[,2]), lineages)[,-1]

# set datasets
var_train <- genotype[train_ind,]
var_test <- genotype[-train_ind,]
lin_train <- lineages[train_ind,]
lin_test <- lineages[-train_ind,]
phen_train <- phenotype$V3[train_ind]
phen_test <- phenotype$V3[-train_ind]
phen_test_w_samples <- phenotype[-train_ind,2:3]
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

json_file <- tempfile(fileext = ".json", tmpdir = "/nfs/research/jlees/jacqueline/gwas_data/klebsiella/gwas_data_bundles", pattern = "oligogenic_EF_100_H_0_9_predictive_leave_one_out_")
cmdstanr::write_stan_json(data = cleaned_gwas_data, file = json_file)

# write file with test phenotypes (ground truth)
write.csv(phen_test_w_samples, file = "/nfs/research/jlees/jacqueline/gwas_data/klebsiella/gwas_data_bundles/oligogenic_EF_100_H_0_9_predictive_leave_one_out_test_phenotypes.csv")

