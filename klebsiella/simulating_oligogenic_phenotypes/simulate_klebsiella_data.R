########################################## WRITE TPED & TFAM FROM SNP MATRIX ##########################################

# Function to convert SNP matrix to TPED and TFAM format in R
snp_matrix_to_tped_tfam <- function(snp_matrix, output_prefix) {
  num_individuals <- nrow(snp_matrix)
  num_snps <- ncol(snp_matrix)
  
  # Ensure column names (SNP IDs) don't have prepended "X"
  rownames(snp_matrix) <- sub("^X", "", rownames(snp_matrix))
  
  # Generate .tfam file
  tfam_file <- paste0(output_prefix, ".tfam")
  tfam_data <- data.frame(
    FID = rep(0, num_individuals),
    IID = rownames(snp_matrix),
    PID = rep(0, num_individuals),
    MID = rep(0, num_individuals),
    Sex = rep(0, num_individuals),
    Phenotype = rep(-9, num_individuals)
  )
  write.table(tfam_data, tfam_file, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # Generate .tped file
  tped_file <- paste0(output_prefix, ".tped")
  tped_data <- data.frame(
    Chromosome = rep(1, num_snps),
    SNP_ID = colnames(snp_matrix),
    Genetic_Distance = rep(0, num_snps),
    Basepair_Position = 1:num_snps
  )
  
  # Add genotype information
  tped_genotypes <- t(snp_matrix)
  tped_genotypes[tped_genotypes == 0] <- "A A"  # ref
  tped_genotypes[tped_genotypes == 1] <- "T T"  # alt
  tped_genotypes[which(!(tped_genotypes == "A A" | tped_genotypes == "T T"))] <- "0 0"
  
  tped_data <- cbind(tped_data, tped_genotypes)
  
  write.table(tped_data, tped_file, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
}

snp_matrix <- read.delim(file = "~/Desktop/simulating_klebsiella_data/pres_abs_matrix.txt",
                         row.names = 1)
snp_matrix <- t(snp_matrix)
output_prefix <- "~/Desktop/simulating_klebsiella_data/simulated_klebsiella_SNPs"
snp_matrix_to_tped_tfam(snp_matrix, output_prefix)

########################################## CHOOSING CAUSAL SNPS ##########################################

# load snp matrix
snp_matrix <- read.delim(file = "~/Desktop/simulating_klebsiella_data/pres_abs_matrix.txt",
                         row.names = 1)
snp_matrix <- t(snp_matrix)

# check which variants are present in multiple lineages
poppunk_clusters <- read.delim("~/Desktop/simulating_klebsiella_data/PP_clusters.txt",
                               sep = "\t", header = F)

# get a random sample of lineages with a population of samples of those lineages
# that is roughly half of the overall population
# these will become the resistant/case phenotype lineages
length(unique(poppunk_clusters[,2]))
set.seed(34)
resistant_lineages <- sample(unique(poppunk_clusters[,2]), size = 170, replace = F)
length(which(poppunk_clusters[,2]%in%resistant_lineages))

# identify the samples belonging to the target lineages
poppunk_clusters$resistance <- "not_resistant"
poppunk_clusters$resistance[which(poppunk_clusters$V2 %in% resistant_lineages)] <- "resistant"

# subset the SNP frame to only resistant samples
rownames(snp_matrix) <- sub("^X", "", rownames(snp_matrix))
snp_matrix_resistant <- snp_matrix[match(poppunk_clusters[which(poppunk_clusters$resistance=="resistant"),1], rownames(snp_matrix)),]
nrow(snp_matrix_resistant)
snp_matrix_nonresistant <- snp_matrix[match(poppunk_clusters[which(poppunk_clusters$resistance=="not_resistant"),1], rownames(snp_matrix)),]
nrow(snp_matrix_nonresistant)

# get MAF of SNPs in resistant and nonresistant samples
maf_resist <- colMeans(snp_matrix_resistant)
maf_nonresist <- colMeans(snp_matrix_nonresistant)

maf_resist_potential <- which(maf_resist > 0.10 & maf_nonresist < 0.05)
maf_resist_potential <- sort(maf_resist[maf_resist_potential]/maf_nonresist[maf_resist_potential], decreasing = T)

top_8 <- maf_resist_potential[1:8]
top_20 <- maf_resist_potential[1:20]

# also want to ensure that these snps are present together in decent LD
# otherwise it will be hard to simulate a phenotype that fits all of them well
potential_causal <- snp_matrix[,which(colnames(snp_matrix)%in%names(maf_resist_potential))]

# reuse some old D'-calculating code
V <- as.numeric(ncol(potential_causal))
N <- as.numeric(nrow(potential_causal))

### Define functions

# Calculate the allele frequencies
find_allele_frequencies <- function(data, labels = TRUE) {
  allele_frequencies <- matrix(data=NA, nrow = V, ncol = 2)
  
  allele_frequencies[,2] <- colSums(data)
  allele_frequencies[,1] <- N - allele_frequencies[,2]
  allele_frequencies <- allele_frequencies/N
  
  if (labels == TRUE) {
    row.names(allele_frequencies) <- colnames(data)
    colnames(allele_frequencies) <- c("0", "1")
    return(allele_frequencies)
  } else {
    return(allele_frequencies)
  }
}

# Calculate haplotype frequencies for a given pair of variants
find_haplotype_frequencies <- function(data, variantA, variantB, compact = TRUE) {
  if (compact == FALSE) {
    haplotype_frequencies <- as.data.frame(table(data[,variantA], data[,variantB]))
    haplotype_frequencies[,3] <- haplotype_frequencies[,3]/N
    return(haplotype_frequencies)
  } else {
    haplotype_frequencies <- table(data[,variantA], data[,variantB])/N
    return(haplotype_frequencies)
  }
}

# Calculate D_prime
#standalone function
#note that variantA and variantB are the index number of the variant (1:V)

find_D_prime <- function(data, variantA, variantB, hap_freq_data, allele_freq_data) {
  if (missing(hap_freq_data)) {
    # find haplotype frequencies for the combination
    #haplotype_frequencies <- table(data[,variantA], data[,variantB])/N
    haplotype_frequencies <- table(factor(data[,variantA], levels = c(0,1)), factor(data[,variantB], levels = c(0,1)))/N
  } else {
    haplotype_frequencies <- hap_freq_data
  }
  if (missing(allele_freq_data)) {
    # find allele frequency for both variants
    allele_frequencies <- find_allele_frequencies(data)
  } else {
    allele_frequencies <- allele_freq_data
  }
  
  #calculate D
  D <- as.numeric(haplotype_frequencies[1,1]*haplotype_frequencies[2,2] - haplotype_frequencies[2,1]*haplotype_frequencies[1,2])
  #calculate D_max
  if (D < 0) {
    set_elements <- as.numeric(c((-1 * allele_frequencies[variantA,1] * allele_frequencies[variantB,1]) , (-1 * allele_frequencies[variantA,2] * allele_frequencies[variantB,2])))
    D_max <- max(set_elements)
  } else {
    if (D>0) {
      set_elements <- as.numeric(c((allele_frequencies[variantA,1] * allele_frequencies[variantB,2]) , (allele_frequencies[variantA,2] * allele_frequencies[variantB,1])))
      D_max <- min(set_elements)
    } else {
      D_max <- 1
    }
  }
  #calculate D_prime
  D_prime <- D/D_max
  return(D_prime)
}


### Generate pairwise matrix of D_prime values

D_prime_matrix <- matrix(data=NA, nrow = V, ncol = V) # set up matrix
pb = txtProgressBar(min = 0, max = V, initial = 0, style = 3) # set up progress bar

for (i in 1:V) {
  setTxtProgressBar(pb,i) # for progress bar
  for (j in i:V) {
    D_prime_matrix[i,j] <- find_D_prime(data = potential_causal, variantA = i, variantB = j)
  }
  close(pb)
}

colnames(D_prime_matrix) <- colnames(potential_causal)
rownames(D_prime_matrix) <- colnames(potential_causal)

# subset to get just top 8 or 20
top_8_d_primes <- D_prime_matrix[names(top_8), names(top_8)]
top_20_d_primes <- D_prime_matrix[names(top_20), names(top_20)]

# top 8 are almost all in perfect LD
# top 20 are in mostly high LD, but with some lower values (eg 0.22)
# try both in simulations

causal_snps_EF_100_ncausal_8 <- data.frame("causal_snps" = names(top_8), "effect_size" = rep(100, times=8))
causal_snps_EF_100_ncausal_20 <- data.frame("causal_snps" = names(top_20), "effect_size" = rep(100, times=20))

write.table(causal_snps_EF_100_ncausal_8, file = "~/Desktop/simulating_klebsiella_data/causal_snps_EF_100_ncausal_8.txt", sep = "\t", row.names = F, quote = F)
write.table(causal_snps_EF_100_ncausal_20, file = "~/Desktop/simulating_klebsiella_data/causal_snps_EF_100_ncausal_20.txt", sep = "\t", row.names = F, quote = F)

causal_snps_EF_30_ncausal_8 <- data.frame("causal_snps" = names(top_8), "effect_size" = rep(30, times=8))
causal_snps_EF_30_ncausal_20 <- data.frame("causal_snps" = names(top_20), "effect_size" = rep(30, times=20))

write.table(causal_snps_EF_30_ncausal_8, file = "~/Desktop/simulating_klebsiella_data/causal_snps_EF_30_ncausal_8.txt", sep = "\t", row.names = F, quote = F)
write.table(causal_snps_EF_30_ncausal_20, file = "~/Desktop/simulating_klebsiella_data/causal_snps_EF_30_ncausal_20.txt", sep = "\t", row.names = F, quote = F)

########################################## CHECK IDS TO TEST (get cppRATE subset for desktop) ##########################################

which(colnames(snp_matrix) %in% names(top_8))

which(colnames(snp_matrix) %in% names(top_20))

