################################################ LOAD LIBRARIES

library(dplyr)
library(data.table)

################################################ GENERATE DATA BUNDLE (MICs) -- LINEAGE SUBCLUSTER VERSION

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

# test lineage subclustering with Fastlin "coarse" lineage clusters

# get unique lineages
lineage_clusters <- unique(lineages$lineages)

# split clusters by delimiter (.)
split_clusters <- strsplit(lineage_clusters, split = '[.]')

# find max and min numbers of levels to each fastlin cluster designation
depths <- list()
for (i in 1:length(split_clusters)) {
  depths <- append(depths, length(split_clusters[[i]]))
}
max_lineage_levels <- as.numeric(max(unlist(depths)))
min_lineage_levels <- as.numeric(min(unlist(depths)))

# now count from the whole dataset
split_df <- strsplit(lineages$lineages, split = '[.]')

# create a dataframe with NAs for uneven lengths
split_df <- do.call(rbind, lapply(split_df, function(x) {
  length(x) <- max_lineage_levels
  return(x)
}))

# convert to dataframe 
split_df <- as.data.frame(split_df, stringsAsFactors = FALSE)

# find number of non-NAs in each column (coarseness level)
non_na_counts <- colSums(!is.na(split_df))
# find number of samples that will have the same designation in coarse (overall cluster) as in fine (subcluster)
non_na_counts[1]-non_na_counts[2]
# find out how many of each unique 'coarse' lineage there are in the dataset
table(split_df[,1])

# set first column to be lineage cluster and second+ columns to be subcluster
# we retain the first value on the subcluster labels because it doesn't matter (as they're categorical)
lineages$lineages_coarse <- split_df[,1]

# define V = number of variants, L = number of lineages, N = number of samples
V <- as.numeric(ncol(genotype))
L <- as.numeric(length(unique(lineages$lineages_coarse)))
S <- as.numeric(length(unique(lineages$lineages)))
N <- as.numeric(nrow(phenotype))
K <- as.numeric(length(unique(phenotype$RFB_MIC_categorical)))

# convert coarse lineages into one-hot format
lineages_coarse <- lineages[,c(1,8)]
lineages_coarse <- model.matrix(X.sample~as.factor(lineages_coarse), lineages_coarse)[,-1]

# convert lineage subclusters into one-hot format
lineages <- lineages[,c(1,5)]
lineages <- model.matrix(X.sample~as.factor(lineages), lineages)[,-1]

# make data into file
cleaned_gwas_data <- list(phenotype = phenotype$RFB_MIC_categorical,
                          variant_matrix = genotype,
                          lineage_matrix = lineages_coarse,
                          sublineage_matrix = lineages,
                          N=N,
                          V=V,
                          L=L,
                          S=S,
                          K=K)

# save completed bundle, ready for gwas, to the appropriate folder
# we use stan json format as other formats (such as RDump) aren't supported by cmdstan

json_file <- tempfile(fileext = ".json", tmpdir = "/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/gwas_bundles/", pattern = "tb_mhunt_rfb_MICs_sublineages_")
cmdstanr::write_stan_json(data = cleaned_gwas_data, file = json_file)

# create index of variant positions
index <- as.matrix(colnames(genotype))
colnames(index) <- c("variant")
write.csv(x = index,
          file = "/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/gwas_bundles/tb_mhunt_rfb_MICs_sublineages_variant_index.csv")

