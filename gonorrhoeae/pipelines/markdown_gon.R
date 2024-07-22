### Pipeline for ordinal GWAS model
library("knitr")
library("vcfR")

#genotype

#USA
```{r USA}
genotype <- read.vcfR("/nfs/research/jlees/jacqueline/gwas_datasets/gonorrhoeae/genotype/variant_calls/USA_all.vcf.gz")
genotype <- extract.gt(genotype, as.numeric = TRUE, IDtoRowNames = T)
genotype <- t(genotype)

#missing2ref
genotype[is.na(genotype)] <- 0

#rename rows for later data sorting
genotype_rownames <- rownames(genotype)
genotype_rownames <- substring(genotype_rownames, 9)
rownames(genotype) <- genotype_rownames
```

#UK
genotype <- read.vcfR("/nfs/research/jlees/jacqueline/gwas_datasets/gonorrhoeae/genotype/variant_calls/USA_all.vcf.gz")
genotype <- extract.gt(genotype, as.numeric = TRUE, IDtoRowNames = T)
genotype <- t(genotype)

#missing2ref
genotype[is.na(genotype)] <- 0

#rename rows for later data sorting
genotype_rownames <- rownames(genotype)
genotype_rownames <- substring(genotype_rownames, 9)
rownames(genotype) <- genotype_rownames

# lineages

#USA
lineages <- read.csv("/nfs/research/jlees/jacqueline/gwas_datasets/gonorrhoeae/poppunk_clusters_USA/poppunk_clusters_USA_clusters.csv")
#UK/Canada
#lineages <- read.csv("/nfs/research/jlees/jacqueline/gwas_datasets/gonorrhoeae/poppunk_clusters_UK/poppunk_clusters_UK_clusters.csv")


# phenotype
phenotype <- read.csv("/nfs/research/jlees/jacqueline/gwas_datasets/gonorrhoeae/phenotype/41564_2018_309_MOESM3_ESM.csv", header = TRUE)
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


# Reorder samples in genotype & lineages file to match phenotype
# genotype samples from the non-desired phenotype dataset (USA or UK/Canada) will be removed

sample_index <- match(phenotype[,2], lineages[,1])
#there is one sample in the lineage data that's not in the phenotype data - remove it from lin
removethese <- which(is.na(sample_index))
phenotype <- phenotype[-removethese,]
#redo sample index and sort as usual
sample_index <- match(phenotype[,2], lineages[,1])
lineages <- lineages[sample_index,]


### Concat into one file (to check for correctness)

# define V = number of variants, L = number of lineages, N = number of samples
V <- as.numeric(ncol(genotype))
L <- as.numeric(length(unique(lineages$Cluster)))
N <- as.numeric(nrow(phenotype))

# DONT FORGET to change the phenotype column to the drug you're studying!
K <- as.numeric(length(unique(phenotype[,3])))

# check that ordinal categories make sense
ordered_categories <- sort(as.numeric(unique(phenotype[,3])))
print(ordered_categories)

ordered_phenotypes <- as.numeric(phenotype[,3])

ordered_phenotypes <- replace(ordered_phenotypes, which(ordered_phenotypes == 0.03), values = 1)
ordered_phenotypes <- replace(ordered_phenotypes, which(ordered_phenotypes == 0.06), values = 2)
ordered_phenotypes <- replace(ordered_phenotypes, which(ordered_phenotypes == 0.125), values = 3)
ordered_phenotypes <- replace(ordered_phenotypes, which(ordered_phenotypes == 0.250), values = 4)
ordered_phenotypes <- replace(ordered_phenotypes, which(ordered_phenotypes == 0.5), values = 5)
ordered_phenotypes <- replace(ordered_phenotypes, which(ordered_phenotypes == 1), values = 6)
ordered_phenotypes <- replace(ordered_phenotypes, which(ordered_phenotypes == 2), values = 7)
ordered_phenotypes <- replace(ordered_phenotypes, which(ordered_phenotypes == 4), values = 8)
ordered_phenotypes <- replace(ordered_phenotypes, which(ordered_phenotypes == 8), values = 9)
ordered_phenotypes <- replace(ordered_phenotypes, which(ordered_phenotypes == 16), values = 10)
ordered_phenotypes <- replace(ordered_phenotypes, which(ordered_phenotypes == 32), values = 11)
ordered_phenotypes <- replace(ordered_phenotypes, which(ordered_phenotypes == 64), values = 12)


print(V)
print(L)
print(N)
print(K)

checkme <- cbind(rownames(genotype), genotype[,1], phenotype, lineages, row.names = rownames(phenotype))
write.csv(checkme, "/nfs/research/jlees/jacqueline/gwas_datasets/gonorrhoeae/pipelines/checkme.csv", row.names=TRUE)


### Read in appropriate Stan model

library(cmdstanr)

mod <- cmdstan_model("/nfs/research/jlees/jacqueline/gwas_datasets/gonorrhoeae/pipelines/BBGWAS-ord-association_gon.stan")

#set_cmdstan_path("/hps/software/users/jlees/jacqueline/anaconda3/envs/cmdstan/bin/cmdstan")
#mod <- cmdstan_model("/nfs/research/jlees/jacqueline/gwas_code/BBGWAS-ord-association.stan") # for continuous association model

# convert lineages into one-hot format
lineages <- model.matrix(lineages[,1]~as.factor(lineages[,2]), lineages)[,-1]

# make data into file
dat <- list(phenotype = ordered_phenotypes,
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