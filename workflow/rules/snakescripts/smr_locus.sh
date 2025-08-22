#!/usr/bin/env bash

set -Eeuxo pipefail

exec 2> "${snakemake_log[0]}"  # send all stderr from this script to the log file

PREFIX_BFILE_CHR="${snakemake_params[prefix_bfile_chr]}"
PREFIX_BESD_CHR="${snakemake_params[prefix_besd_chr]}"
PREFIX_OUT="${snakemake_params[prefix_out]}"

# cleanup() {
#     rm -f "${PREFIX_OUT}"*
# }

# trap cleanup ERR

PROBE="${snakemake_wildcards[probe]}"
i=$(awk -v probe="$PROBE" '$4 == probe {print $1}' "${snakemake_input[genelist]}")


smr --bfile "${PREFIX_BFILE_CHR}${i}" \
    --gwas-summary ${snakemake_input[ma]} \
    --beqtl-summary "${PREFIX_BESD_CHR}${i}" \
    --probe ${snakemake_wildcards[probe]} \
    --gene-list ${snakemake_input[genelist]} \
    --probe-wind ${snakemake_params[probe_wind]} \
    --plot \
    --thread-num ${snakemake[threads]} \
    --out "${PREFIX_OUT}"  \
    2>&1 | tee "${snakemake_log[0]}" 

# Merge the results
# awk 'NR == FNR || FNR > 1' "${PREFIX_OUT_CHR}"*.msmr > ${snakemake_output[smr]}
# awk 'NR == FNR || FNR > 1' "${PREFIX_OUT_CHR}"*.snp_failed_freq_ck.list > ${snakemake_output[fail]}
