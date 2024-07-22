### Pipeline for ordinal GWAS model

library("vcfR")

#genotype

#USA
genotype <- read.vcfR("/nfs/research/jlees/jacqueline/gwas_datasets/gonorrhoeae/genotype/variant_calls/USA_all.vcf.gz")
genotype <- extract.gt(genotype, as.numeric = TRUE, IDtoRowNames = T)
genotype <- t(genotype)

#UK
#genotype <- read.vcfR("/nfs/research/jlees/jacqueline/gwas_datasets/gonorrhoeae/genotype/variant_calls/UK_all.vcf.gz")
#genotype <- extract.gt(genotype, as.numeric = TRUE, IDtoRowNames = T)
#genotype <- t(genotype)

#missing2ref
genotype[is.na(genotype)] <- 0

#rename rows for later data sorting
genotype_rownames <- rownames(genotype)
genotype_rownames <- substring(genotype_rownames, 9)
rownames(genotype) <- genotype_rownames


# lineages

#USA
lineages <- read.csv("/nfs/research/jlees/jacqueline/gwas_datasets/gonorrhoeae/lineages/poppunk_clusters_USA/poppunk_clusters_USA_clusters.csv")
#UK/Canada
#lineages <- read.csv("/nfs/research/jlees/jacqueline/gwas_datasets/gonorrhoeae/poppunk_clusters_UK/poppunk_clusters_UK_clusters.csv")


# phenotype
phenotype <- read.csv("/nfs/research/jlees/jacqueline/gwas_datasets/gonorrhoeae/pipelines/41564_2018_309_MOESM3_ESM.csv", header = TRUE)
colnames(phenotype) <- phenotype[1,]
phenotype <- phenotype[-1,]

#get only desired phenotype dataset

#index separate countries
USA <- which(phenotype[,1] == "USA")
#UK <- which(phenotype[,1] == "England")
#Canada <- which(phenotype[,1] == "Canada")
#UK_Canada <- c(UK, Canada)

#USA dataset
phenotype <- phenotype[USA,]

#UK/Canada dataset
#phenotype <- phenotype[UK_Canada,]

# choose which phenotype column (which antibiotic) you want to test
chosen_antibiotic <- 3

# Reorder samples in genotype & lineages file to match phenotype
# genotype samples from the non-desired phenotype dataset (USA or UK/Canada) will be removed

sample_index <- match(phenotype[,2], lineages[,1])
#there is one sample in the lineage data that's not in the phenotype data - remove it from lineages
removethese <- which(is.na(sample_index))
phenotype <- phenotype[-removethese,]
#redo sample index and sort as usual
sample_index <- match(phenotype[,2], lineages[,1])
lineages <- lineages[sample_index,]


# Define V = number of variants, L = number of lineages, N = number of samples
L <- as.numeric(length(unique(lineages$Cluster)))
V <- as.numeric(ncol(genotype))
N <- as.numeric(nrow(phenotype))
K <- as.numeric(length(unique(phenotype[,chosen_antibiotic])))

# print for future reference
print(V)
print(L)
print(N)
print(K)

# Check that the phenotypes are ordered correctly
ordered_categories <- sort(as.numeric(unique(phenotype[,chosen_antibiotic])))
print(ordered_categories)

# Convert the chosen phenotype to its ordinal category
categorical_phenotype <- as.numeric(phenotype[,chosen_antibiotic])
for (i in 1:N) {
  categorical_phenotype[i] <- which(ordered_categories == categorical_phenotype[i])
}

# Concat into one file (to check for correctness)
checkme <- cbind(rownames(genotype), genotype[,1], phenotype, categorical_phenotype, lineages, row.names = rownames(phenotype))

### Read in appropriate Stan model

library(cmdstanr)
mod <- cmdstan_model("/nfs/research/jlees/jacqueline/gwas_datasets/gonorrhoeae/pipelines/BBGWAS-ord-association_gon.stan")

# convert lineages into one-hot format
lineages <- model.matrix(lineages[,1]~as.factor(lineages[,2]), lineages)[,-1]

# make data into file
dat <- list(phenotype = categorical_phenotype,
            variant_matrix = genotype,
            lineage_matrix = lineages,
            N=N,
            V=V,
            L=L,
            K=K)

# Variational inference algorithm
fit <- mod$variational(dat,
                       seed = 3930,
                       refresh = 500)

### Save fitted model object
temp_rds_file <- tempfile(pattern = "gonorrhoeae_USA_PEN_", fileext = ".RDS", tmpdir = "/nfs/research/jlees/jacqueline/gwas_analyses/gonorrhoeae/")
fit$save_object(file = temp_rds_file)
