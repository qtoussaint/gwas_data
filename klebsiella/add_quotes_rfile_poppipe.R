old_rfile <- read.delim(file = "/nfs/research/jlees/jacqueline/gwas_data/klebsiella/lineages/rfile_poppipe.tsv",
           header = F)

write.table(old_rfile, 
            file = "/nfs/research/jlees/jacqueline/gwas_data/klebsiella/lineages/rfile_with_quotes.tsv", 
            sep="\t",
            row.names = F, 
            col.names = F)
