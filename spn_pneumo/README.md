# Streptococcus pneumoniae GWAS Data

## Data source
SPARC cohort. Sequence assemblies (FASTA contigs) for each sample.
Phenotype metadata: `phenotype/sparc_metadata.csv`

---

## Pipeline steps

### Step 1: Variant calling
Script: `code/gwas_data/spn_pneumo/01_genotyping/snippy_multi.sh`
- Runs snippy on each sample assembly against `ATCC_700669_genbank.gbff` reference
- Output: `genotype/snippy/results_*/snps.vcf.gz` (per-sample VCFs)
- Note: `genotyping.sh` (uses bactsnp) also exists in this directory but was not used;
  the output VCF files (`snps_snippy_*`) confirm snippy was the approach used.

### Step 2: VCF merging, filtering, and annotation
Script: `code/gwas_data/spn_pneumo/02_filtering_annotation/snippy_bcftools_snpEff.sh`
- Merges per-sample VCFs with bcftools merge (missing genotypes → ref)
- Fixes ploidy (diploid → haploid)
- Annotates with snpEff (Streptococcus_pneumoniae_atcc_700669)
- Extracts field annotations with SnpSift (effect, impact, gene, HGVS)
- Filters MAF > 0.05; removes MODIFIER-only sites using BED intersection
- Note: most commands are commented out (earlier steps); the active (uncommented) section
  generates the `NOMODIFIERS_NOSYNONYMOUS` genotype files
- Key outputs (in `genotype/`):
  - `snps_snippy_maf05_NOMODIFIERS_multiallelic.vcf.gz` — MAF>0.05, positions with ≥1 non-MODIFIER effect
  - `genotype_maf05_NOMODIFIERS_multiallelic.tsv` — presence/absence matrix ← **USED IN BAYESIAN GWAS PAPER**
  - `genotype_maf05_NOMODIFIERS_NOSYNONYMOUS_multiallelic.tsv` — additionally removes synonymous/LOW impact sites
  - `fields_filtered_maf05_multiallelic.txt` — snpEff field annotations

### Step 3: Population structure (subclustering)
Tool: PopPIPE (wraps PopPunk)
Config notes: `code/gwas_data/spn_pneumo/03_lineages/poppipe_sublineages.sh`
- This is a notes file, not a runnable script
- Input: PopPunk database at `lineages/poppunk/`, rfile at `lineages/poppunk/rfile.txt`
- Final run used `min_cluster_size=3` (the notes file records `min_cluster_size=6`, which is outdated)
- Output: `lineages/sublineages_from_poppipe/min_cluster_3/all_clusters.txt` ← **USED IN BAYESIAN GWAS PAPER**
- Also available: coarse clusters at `lineages/database_clusters.csv` (used in older scripts only)

### Step 4: Binary phenotype generation
Script: `code/gwas_data/spn_pneumo/04_phenotype/phenotyping_binary_eucast.R`
- Input: `phenotype/sparc_metadata.csv`
- Applies EUCAST 2025 benzylpenicillin breakpoint: S if MIC < 0.06 μg/mL, R if ≥ 0.06
- Output: `phenotype/sparc_penicillin_eucast_binary.csv` ← **USED IN BAYESIAN GWAS PAPER**

### Step 5: Data bundling (in bayesian_gwas_paper)
Script: `bayesian_gwas_paper/code/bayesian_gwas_paper/00_data_preprocessing/write_inference_jsons.R`

Inputs:
- `genotype/genotype_maf05_NOMODIFIERS_multiallelic.tsv`
- `lineages/sublineages_from_poppipe/min_cluster_3/all_clusters.txt`
- `phenotype/sparc_metadata.csv` (columns: `Benzylpenicillin.MIC..ug.mL.`, `Trimethoprim.MIC..ug.mL.`)
- `phenotype/sparc_penicillin_eucast_binary.csv`

MIC cleaning applied:
- `>4.00` → 4, `<=0.03` → 0.03, `<0.016` → 0.016, `< 0.01` → 0.01

Output JSONs in `bayesian_gwas_paper/00_data/inference/`:
| JSON | Description |
|------|-------------|
| `01_spn_penicillin_subclusters_K5` | 5 MIC bins, breakpoints [0, 0.032, 0.065, 0.2, 2, 6] |
| `02_spn_penicillin_subclusters_binary` | Binary S/R classification |
| `03_spn_penicillin_subclusters_continuous` | Continuous log MIC |
| `04_spn_trimethoprim_subclusters_K9` | 9 MIC bins, breakpoints [0, 0.04, 0.08, 0.12, 0.16, 0.24, 0.32, 0.4, 0.5, 2] |
| `05_spn_trimethoprim_subclusters_binary` | Binary S/R classification |
| `06_spn_trimethoprim_subclusters_continuous` | Continuous log MIC |
| `07_spn_penicillin_clusters_binary` | Binary, coarse PopPunk clusters (not sublineages) |

---

## Pyseer GWAS (also used in paper)
Scripts: `code/gwas_data/spn_pneumo/05_pyseer/`
- `data_prep_pyseer.R` — prepares genotype/phenotype files for pyseer
- `run_pyseer.sh` — runs pyseer LMM with phylogenetic similarity correction
- `mash.sh` — calculates mash distances for phylogeny
- `scripts/phylogeny_distance.py` — helper script
Output: `pyseer_data/`

---

## Outdated scripts
`code/gwas_data/spn_pneumo/outdated/`

| Script | What it does | Why outdated |
|--------|--------------|--------------|
| `data_bundling_spn_sparc.R` | Bundles data for cmdstanr | Reads old `sparc_snps.vcf.gz`, uses coarse `database_clusters.csv` (not sublineages), outputs to `gwas_bundles/`; superseded by `write_inference_jsons.R` |
