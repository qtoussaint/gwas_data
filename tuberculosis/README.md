# Mycobacterium tuberculosis GWAS Data

## Data source
CRyPTIC consortium regenotyped VCFs (release June 2022).
- Raw VCF paths: `/nfs/ftp/public/databases/cryptic/release_june2022/reproducibility/...`
  (accessed via datamover partition; 12,287 per-sample VCFs)
- Full file list: `file.txt`
- Phenotype metadata: `CRyPTIC_reuse_table_20240917.csv` (most recent version; older scripts use `20231208.csv`)

---

## Pipeline steps

### Step 1: Lineage assignment
Script: `code/gwas_data/tuberculosis/01_lineages/fastlin_lineage_clustering.sh`
- Tool: fastlin with MTBC barcodes
- Barcode file: `01_lineages/MTBC_barcodes_for_fastlin.tsv`
- Input: sequence reads (in `/nfs/research/jlees/jacqueline/gwas_datasets/tuberculosis/genotype/sequence_reads`)
  Note: path is `gwas_datasets`, not `gwas_data`
- Note: uses older LSF job scheduler (`#BSUB` headers)
- Output: `01_lineages/output_fastlin.txt` ← **USED IN BAYESIAN GWAS PAPER**
- Samples with multiple lineage designations (~300 contaminated samples) are excluded downstream
  inside `write_inference_jsons.R`

Step 0 (optional/old): `01_lineages/sequence_reads_download_pipeline.sh` — downloads CRyPTIC reads
from public repositories (only needed if reads are not already present)

### Step 2: Merge CRyPTIC VCFs + filter + annotate
Script: `code/gwas_data/tuberculosis/02_vcf_processing/cryptic_chunked_merge.sh`
(SLURM array job, run on datamover partition)
- Merges 12,287 per-sample regenotyped VCFs listed in `file.txt` in 20 genomic chunks
- Renames chromosome (`NC_000962.3`) to match snpEff DB name
- Fixes ploidy (diploid → haploid)
- Adds MAF tag → filters MAF > 0.0004
- Concatenates 20 filtered chunks → `cryptic_regeno_MAF_0.0004.allchunks.vcf.gz`
- Annotates with snpEff (Mycobacterium_tuberculosis_h37rv)
- Converts missing genotypes to 0, any alt call to 1 (binary presence/absence)
- Extracts snpEff field annotations (effect, impact, gene, HGVS)
- Note: most commands in the script are commented out (earlier runs); the active section
  handles the final steps
- Key outputs in `cryptic_regeno_snpeff/`:
  - `cryptic_regeno_MAF_0.0004.allchunks_effannotated_binary.vcf.gz`
  - `fields_all.txt`, `fields_filtered.txt`

### Step 3: Extract presence/absence matrix
Scripts:
- `code/gwas_data/tuberculosis/03_presence_absence/presence_absence_from_VCF.R`
- `code/gwas_data/tuberculosis/03_presence_absence/presence_absence_from_VCF_step2.sh`

- Reads `cryptic_regeno_MAF_0.0004.allchunks_effannotated_binary.vcf.gz`
- Creates variant IDs as `CHROM_POS_REF_ALT`
- Output: `cryptic_regeno_snpeff/cryptic_regeno_MAF_0.0004_effannotated_presabs.tsv` (variants × samples)
- `presence_absence_from_VCF_step2.sh` then:
  - Transposes with datamash → samples × variants
  - Reorders samples to match `sample_index.csv`
  - Output: `cryptic_regeno_snpeff/subset_K8.tsv`
    (samples that have RIF MIC data AND lineage AND genotype data)

### Step 4: Filter to RIF K8 sample subset
Script: `code/gwas_data/tuberculosis/03_presence_absence/filter_samples.sh`
- `sample_index_RIF_K8.txt` — created manually from `CRyPTIC_reuse_table_20240917.csv`;
  contains all samples that have RIF MIC metadata available (no generating script)
- Filters VCF to those samples, extracts presence/absence
- Outputs in `cryptic_regeno_snpeff/subsampled_presabs/`:
  - `presence_absence_final_RIF_K8.tsv` — samples × variants (rows have no header)
    ← **USED IN BAYESIAN GWAS PAPER**
  - `variant_index_presence_absence_RIF_K8.tsv` — variant names only

Note on the two output files:
- `presence_absence_final_RIF_K8.tsv` (step 4) = all samples with RIF MIC metadata
- `subset_K8.tsv` (step 3) = further filtered to also have lineage and genotype data
- The bayesian_gwas_paper uses `presence_absence_final_RIF_K8.tsv`; lineage filtering
  happens inside `write_inference_jsons.R`

### Step 5: Convert TSV to JSON
Scripts: `cryptic_regeno_snpeff/subsampled_presabs/tsv_to_cmdstan_json.py` / `tsv_to_cmdstan_json.sh`
- These scripts live inside the data subdirectory (not in `code/`) because they are run there
- Converts `presence_absence_final_RIF_K8.tsv` → `presence_absence_final_RIF_K8.json`
- JSON contains only `{"variant_matrix": [[...], ...]}` (genotype matrix only; other fields
  are concatenated in step 6)

### Step 6: Data bundling (in bayesian_gwas_paper)
Script: `bayesian_gwas_paper/code/bayesian_gwas_paper/00_data_preprocessing/write_inference_jsons.R`

Inputs:
- `cryptic_regeno_snpeff/subsampled_presabs/presence_absence_final_RIF_K8.tsv` (ordinal/continuous)
  or `presence_absence_final_RIF_K8.json` (binary — genotype matrix only, other fields concatenated)
- `cryptic_regeno_snpeff/subsampled_presabs/variant_index_presence_absence_RIF_K8.tsv`
- `CRyPTIC_reuse_table_20240917.csv` (RIF phenotype)
- `code/gwas_data/tuberculosis/01_lineages/output_fastlin.txt` (top-level lineage: "1", "2", "3", "4")

MIC cleaning applied:
- `<=0.06` → 0.06, `>2` → 2, `>4` → 4, `>8` → 8

Output JSONs in `bayesian_gwas_paper/00_data/inference/`:
| JSON | Description |
|------|-------------|
| `16_tb_rifampicin_subclusters_K8` | 8 unique MIC values mapped to categories 1:8 |
| `17_tb_rifampicin_subclusters_binary` | Binary S/R classification |
| `18_tb_rifampicin_subclusters_continuous` | Continuous log MIC |
| `19_tb_rifampicin_subclusters_K5` | 5 MIC bins |

---

## Outdated scripts
`code/gwas_data/tuberculosis/outdated/`

| Script | What it does | Why outdated |
|--------|--------------|--------------|
| `data_bundling_tb.R` | JSON bundling | Uses older `tb_mhunt_genotypes_cleaned.txt` dataset (not CRyPTIC) and rifabutin (not rifampicin); outputs to `gwas_bundles/`; superseded |
| `data_bundling_tb_sublineage.R` | JSON bundling with sublineages | Same as above |
| `pyseer_data_prep.R` | Pyseer data prep | Uses `tb_mhunt_genotypes_cleaned.txt`; outputs to `pyseer_data/`; superseded |
| `run_pyseer_tb.sh` | Runs pyseer GWAS | Pyseer pipeline using parsnp tree; superseded |
| `run_shovill_for_pyseer.sh` | Assembly from reads | Pyseer-specific prep; superseded |
| `run_fasttree_for_pyseer.sh` | Builds phylogenetic tree | Pyseer-specific prep; superseded |
| `generate_sample_ids.sh` | Generates sample ID list | Associated with pyseer pipeline; superseded |
| `rename_assemblies.sh` | Renames assembly files | Associated with pyseer pipeline; superseded |

Legacy scripts are also in `old/` (kept as-is):
- `old/snippy_pipeline.sh`, `old/bcftools_pipeline.sh` — old variant calling; superseded by CRyPTIC pre-existing VCFs
- `old/create_genotype_matrix.sh`, `old/create_gt_matrix.sh` — old genotype extraction; superseded
