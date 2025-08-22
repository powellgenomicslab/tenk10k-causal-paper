# Purpose: Compute scDRS cell scores for a given phenotype using preprocessed data and gene sets
# Input: Preprocessed AnnData (pickle), gene set file, config YAML
# Output: Cell score and cell type statistics files

import scdrs
import scanpy as sc
import pandas as pd
from pathlib import Path
import numpy as np
import yaml
import pickle

INPUT = snakemake.input
OUTPUT = snakemake.output

REG_PKL = INPUT['reg_pkl']
GS = INPUT['gs']
CONFIG = INPUT['config']
OUT_CELL = OUTPUT['cell']
OUT_CELLTYPE = OUTPUT['celltype']
TRAIT = snakemake.wildcards.phenotype

# Interactive test
# REG_PKL = "resources/scdrs/h5ad/tenk10k_phase1.reg.h5ad.pkl" 
# GS = "resources/scdrs/gs_chunked/tenk10k_phase1/crohns.gs"
# CONFIG = "resources/scdrs/config/tenk10k_phase1.yaml"
# OUT_CELL =  "results/scdrs/cell_score/tenk10k_phase1/crohns.cell_score.tsv.parquet.gz"
# OUT_CELLTYPE =  "results/scdrs/cell_type_stats/tenk10k_phase1/crohns.cell_type_stats.tsv"
# Load config
with open(CONFIG, 'r') as f:
    config = yaml.safe_load(f)

# Load preprocessed pickled data
with open(REG_PKL, 'rb') as f:
    adata = pickle.load(f)

# Load geneset
dict_gs = scdrs.util.load_gs(
    GS,
    src_species=config['gs_src_species'],
    dst_species=config['gs_dst_species'],
    to_intersect=adata.var_names,
)

gene_list, gene_weights = dict_gs[TRAIT]

df_score = scdrs.score_cell(
        data=adata,
        gene_list=gene_list,
        gene_weight=gene_weights,
        ctrl_match_key="mean_var",
        n_ctrl=1000,
        weight_opt="vs",
        return_ctrl_raw_score=False,
        return_ctrl_norm_score=True,
        verbose=False,
    )

# cols_output = ['raw_score', 'norm_score', 'mc_pval', 'pval', 'nlog10_pval', 'zscore']

# write to output
Path(OUT_CELL).parent.mkdir(parents=True, exist_ok=True)
df_score \
    .to_parquet(OUT_CELL, compression = "gzip", index=True)

# compute cell type stats
df_stats = scdrs.method.downstream_group_analysis(
    adata=adata,
    df_full_score=df_score,
    group_cols=[config['group_stats_col']],
)[config['group_stats_col']] 

Path(OUT_CELLTYPE).parent.mkdir(parents=True, exist_ok=True)
df_stats \
    .reset_index() \
    .to_csv(OUT_CELLTYPE, sep='\t', index=False)