### generate a binary phenotype based on EUCAST guidelines
# not sure why 8 samples are missing in the pyseer resistance dataframe or how it was made, so making my own binary from the sparc metadata

# for benzylpenicillin, 2025 eucast breakpoints are:
#S if <0.06
#I if 0.06-1
#R if >1

# i'm going to go <0.06 and >0.06

phenotype <- read.csv(file = "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/phenotype/sparc_metadata.csv", na.strings = "")

penicillin <- data.frame(samples = phenotype$Lane.name, MIC = phenotype$Benzylpenicillin.MIC..ug.mL.)
penicillin <- penicillin[complete.cases(penicillin),] # remove NAs

# show original values
print("True MIC values:")
print(sort(unique(penicillin$MIC)))

# clean duplicate values and remove symbols
penicillin$MIC <- replace(penicillin$MIC, penicillin$MIC==">4.00", 4)
penicillin$MIC <- replace(penicillin$MIC, penicillin$MIC==">4", 4)
penicillin$MIC <- replace(penicillin$MIC, penicillin$MIC=="<=0.03", 0.03)
penicillin$MIC <- replace(penicillin$MIC, penicillin$MIC=="<0.03", 0.03)
penicillin$MIC <- replace(penicillin$MIC, penicillin$MIC=="<0.016", 0.016)
penicillin$MIC <- replace(penicillin$MIC, penicillin$MIC=="<.016", 0.016)
penicillin$MIC <- replace(penicillin$MIC, penicillin$MIC=="< 0.01", 0.01)

# make numeric
penicillin$MIC <- as.numeric(penicillin$MIC)

# test different bin breaks and make the cuts
# 4 bins (from 26 total)
MICs <- sort(unique(penicillin$MIC))
print(unique(cut(MICs, breaks=c(0, 0.065, 10))))

# cut according to give breakpoints and add as separate column
penicillin_complete <- penicillin %>% mutate(MIC_factor = cut(penicillin$MIC, breaks=c(0, 0.065, 10)))

# replace factors with characters
penicillin_complete$MIC_factor <- print(penicillin_complete$MIC_factor)

# create sorted vector of all MIC bins
MIC_intervals <- as.character(sort(unique(penicillin_complete$MIC_factor)))

# convert MIC to sample,s ordinal category (1:24) in resistances
for (i in 1:nrow(penicillin_complete)) {
  penicillin_complete$MIC_ord_category[i] <- which(MIC_intervals == penicillin_complete$MIC_factor[i])
}

penicillin_complete$MIC_binary_category <- as.character(as.numeric(penicillin_complete$MIC_ord_category)-1)

# switch order of columns and rename
penicillin <- data.frame(ID=penicillin_complete$samples, resistance=penicillin_complete$MIC_binary_category)

write_csv(x = penicillin, file = "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/phenotype/sparc_penicillin_eucast_binary.csv")



eucast <- read.csv(file = "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/phenotype/sparc_penicillin_eucast_binary.csv")
pyseer <- read.delim(file = "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/pyseer_data/pyseer_tutorial/resistances.pheno", na.strings = "")

match(eucast$ID, pyseer$samples)

test_difference <- eucast[match(eucast$ID, pyseer$samples),]
test_difference <- test_difference[complete.cases(test_difference),]
test_difference$eucast_ids <- test_difference$ID
test_difference$resistance_eucast <- test_difference$resistance
test_difference <- test_difference[,-c(1:2)]
test_difference$resistance_pyseer <- pyseer$penicillin
test_difference$pyseer_ids <- pyseer$samples

cor(as.numeric(as.character(test_difference$resistance_eucast)), y = as.numeric(as.character(test_difference$resistance_pyseer)), method = "spearman", use = "complete.obs")

# cursed summary
# shrodinger's resistance table
table(eucast=test_difference$resistance_eucast, pyseer=test_difference$resistance_pyseer)

