# activate mash env

# define outdir
outdir=/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype/mash_distances

# sketch genomes
#mash sketch -s 10000 -o $outdir/mash_sketch /nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype/sequences/*.fa

# create distance matrix
#mash dist $outdir/mash_sketch.msh $outdir/mash_sketch.msh > $outdir/mash.tsv

# make square
awk -v od="$outdir/" '{
  gsub(od, "", $1)
  gsub(od, "", $2)
  print $1 "\t" $2 "\t" $3
}' $outdir/mash.tsv > $outdir/mash_GRM.tsv

# strip path prefixes
$outdir/read_in_mash.R

# final file is mash_distance_matrix.csv