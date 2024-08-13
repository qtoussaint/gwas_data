################################################ LOAD LIBRARIES

library(dplyr)
library(data.table)

################################################ GENERATE DATA FOR PYSEER 

# read in phenotype
# use newer (2023) table
phenotype <- read.csv(file = "/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/CRyPTIC_reuse_table_20231208.csv")

# pick phenotype
# here I'm using the rifabutin resistance (RFB) phenotype with MICs
phenotype <- phenotype[,c(1,28)]

# remove NAs
#phenotype <- na.omit(phenotype)

# check unique values and put in order
order <- unique(phenotype$continuous)
phenotype[which(phenotype[,2]=="<=0.06"),2] <- "0.06"
phenotype[which(phenotype[,2]==">2"),2] <- "2"
phenotype$continuous <- as.numeric(phenotype$continuous)

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

# transpose back to normal format
genotype <- t(genotype)

# ensure that the genotypes and phenotypes match
# including removing samples not in assemblies dir (due to not being able to be assembled by shovill)
sample_index <- match(colnames(genotype), phenotype$ENA_RUN)
sample_index <- na.omit(sample_index)
phenotype <- phenotype[sample_index,]

assemblies <- read.delim(file = "/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/genotype/assemblies.txt")
assemblies <- t(as.data.frame(strsplit(x = unlist(assemblies), split = ".fa")))

sample_index <- match(assemblies, phenotype$ENA_RUN)
sample_index <- na.omit(sample_index)
phenotype <- phenotype[sample_index,]

sample_index <- match(phenotype$ENA_RUN, colnames(genotype))
sample_index <- na.omit(sample_index)
genotype <- genotype[,sample_index]

# change phenotype header
colnames(phenotype) <- c("sample", "continuous")

# export genotypes as Rtab
write.table(genotype, file = "/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/pyseer_data/genotypes.Rtab", 
            sep = "\t", col.names = T, row.names = T)

# export phenotypes as tab separated
write.table(phenotype, file = "/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/pyseer_data/phenotypes_rfb_mic.txt", 
            sep = "\t",
            col.names = T, row.names = F,
            quote = F)

# export sample list
write.table(phenotype[,1], file = "/nfs/research/jlees/jacqueline/gwas_data/tuberculosis/pyseer_data/sample_list.txt", 
            sep = "\t", 
            col.names = F, row.names = F)
