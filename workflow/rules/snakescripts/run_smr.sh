#!/usr/bin/env bash

set -Eeuxo pipefail

exec 2> "${snakemake_log[0]}"  # send all stderr from this script to the log file

PREFIX_BFILE_CHR="${snakemake_params[prefix_bfile_chr]}"
PREFIX_BESD_CHR="${snakemake_params[prefix_besd_chr]}"
PREFIX_OUT_CHR="${snakemake_params[prefix_out_chr]}"

cleanup() {
    rm -f "${PREFIX_OUT_CHR}"*
}

trap cleanup EXIT ERR

for i in {1..22}; do
    smr --bfile "${PREFIX_BFILE_CHR}${i}" \
        --gwas-summary ${snakemake_input[ma]} \
        --beqtl-summary "${PREFIX_BESD_CHR}${i}" \
        --peqtl-smr $(cat ${snakemake_input[pthresh]}) \
        --extract-probe ${snakemake_input[probe]} \
        --maf ${snakemake_params[maf]} \
        --diff-freq-prop ${snakemake_params[diff_freq_prop]} \
        --smr-multi \
        --thread-num ${snakemake[threads]} \
        --out "${PREFIX_OUT_CHR}${i}" \
        2>&1 | tee "${snakemake_log[0]}" 
done

# Merge the results
awk 'NR == FNR || FNR > 1' "${PREFIX_OUT_CHR}"*.msmr > ${snakemake_output[smr]}
awk 'NR == FNR || FNR > 1' "${PREFIX_OUT_CHR}"*.snp_failed_freq_ck.list > ${snakemake_output[fail]}

