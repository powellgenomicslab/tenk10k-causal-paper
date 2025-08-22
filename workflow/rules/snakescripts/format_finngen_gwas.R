# | Column name             | Description
# | ----------------------- | -------------------
# | #CHR                    | Chromosome
# | POS                     | Position
# | REF                     | Reference base
# | ALT                     | Alternative allele(s)
# | SNP                     | Variant in format #CHR:POS:REF:ALT
# | FINNGEN_beta            | FinnGen summary beta
# | FINNGEN_sebeta          | FinnGen summary standard deviation
# | FINNGEN_pval            | FinnGen summary p-value
# | FINNGEN_af_alt          | FinnGen alternative allele frequency
# | FINNGEN_af_alt_cases    | FinnGen alternative allele frequency in cases
# | FINNGEN_af_alt_controls | FinnGen alternative allele frequency in controls
# | UKBB_beta               | UKBB summary beta
# | UKBB_sebeta             | UKBB summary standard deviation
# | UKBB_pval               | UKBB summary p-value
# | UKBB_af_alt             | UKBB alternative allele frequency
# | ESTBB_beta              | ESTBB summary beta
# | ESTBB_sebeta            | ESTBB summary standard deviation
# | ESTBB_pval              | ESTBB summary p-value
# | ESTBB_af_alt            | ESTBB alternative allele frequency
# | all_meta_N              | Number of studies used for meta-analysis
# | all_inv_var_meta_beta   | Inverse-variance weighted meta-analysis beta
# | all_inv_var_meta_sebeta | Inverse-variance weighted meta-analysis standard deviation
# | all_inv_var_meta_p      | Inverse-variance weighted meta-analysis p-value
# | all_inv_var_meta_mlogp  | Inverse-variance weighted meta-analysis -log10(p-value)
# | all_inv_var_het_p       | Inverse-variance weighted meta-analysis Cochran's Q test p-value
# | leave_[COHORT]_*         | Leave-one-out meta-analysis statistics
# | rsid                    | Reference SNP ID assigned by dbSNP

library(tidyverse)
library(data.table)
library(yaml)

INPUT <- snakemake@input
OUTPUT <- snakemake@output[[1]]
PHENO <- snakemake@wildcards[["pheno"]]
setDTthreads(snakemake@threads)

# for testing / interactive use
# setDTthreads(4)
# INPUT <- list(gwas = "resources/sumstats/finngen_gwas_extract/tenk10k_phase1/cmelanoma.tsv.gz",
#               metadata = "resources/misc/meta_analysis_mvp_ukbb_FinnGen_R12_MVP_UKBB_mapping.tsv",
#               config = "workflow/config/finngen_meta_path.yaml",
#               frq = "resources/genotypes_frq/tenk10k_phase1.frq")
# OUTPUT <- "resources/sumstats/finngen_gwas_extract/tenk10k_phase1/cmelanoma.ma"
# PHENO <- "cmelanoma"
df <- fread(INPUT$gwas)

df_trait_meta <- fread(INPUT$metadata)
pheno_prefix <- read_yaml(INPUT$config)$pheno_prefix[[PHENO]]

df_n <- df_trait_meta %>% 
    filter(fg_phenotype == pheno_prefix) %>% 
    rename_with(~str_replace_all(.x, "_n", ".n")) %>% 
    pivot_longer(contains(".n"),
                 names_to = c("cohort", "n_type"), names_sep = "\\.") %>% 
    filter(!is.na(value))
setDT(df_n)

n_case <- df_n[n_type == "n_cases", sum(value)]
n_control <- df_n[n_type == "n_controls", sum(value)]
n_effective <- floor(4 / (1/n_case + 1/n_control))

# match with tenk10k snp id
df_tenk10k_freq <- fread(INPUT$frq)

df[df_tenk10k_freq,
   `:=`(snp_id = i.SNP, A1 = i.A1, A2 = i.A2),
   on = c("#CHR" = "CHR", "POS")]

df[, `:=`(Amin.x = pmin(REF, ALT),
          Amax.x = pmax(REF, ALT),
          Amin.y = pmin(A1, A2),
          Amax.y = pmax(A1, A2))]

df_out <- df %>% 
    filter(SNP == snp_id | (SNP != snp_id & Amin.x == Amin.y & Amax.x == Amax.y)) %>% 
    # calculate average AF for alt (effect) allele
    mutate(avg_alt_freq = rowMeans(select(., ends_with("af_alt")), na.rm = TRUE),
           b = fifelse(ALT == A1, all_inv_var_meta_beta, -all_inv_var_meta_beta),
           freq = fifelse(ALT == A1, avg_alt_freq, 1-avg_alt_freq),
           n = n_effective) %>% 
    select(SNP = snp_id,  A1, A2, freq, b,
           se = all_inv_var_meta_sebeta, p = all_inv_var_meta_p, n)


# Write results
fwrite(df_out, OUTPUT, sep = "\t", na = "NA", quote = FALSE)