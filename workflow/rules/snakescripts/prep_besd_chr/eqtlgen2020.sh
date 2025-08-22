#!/usr/bin/env bash

set -Eeuxo pipefail

exec 2> "${snakemake_log[0]}" 

PREFIX_BESD="resources/download/eqtl/eqtlgen_vosa2020/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense"
OUTDIR="${snakemake_output[0]}"
THREADS=${snakemake[threads]} 
# OUTDIR="resources/besd/eqtlgen2020"
# THREADS=8 
CHAIN_FILE="resources/misc/hg19ToHg38.over.chain"

mkdir -p "${OUTDIR}/bulk_wb"

# liftover to b38
awk -v OFS='\t' '{print $1, $4, $4, $2}' "${PREFIX_BESD}.esi" > "${PREFIX_BESD}.esi.b37"

# conda init && conda activate pydata
CrossMap bed --chromid s $CHAIN_FILE "${PREFIX_BESD}.esi.b37" "${PREFIX_BESD}.esi.b38"

# extract snp list available in b38 and in tenk10k
awk -v PREFIX="${OUTDIR}/bulk_wb/chr" \
    '$1 <= 22 {print $4 > PREFIX $1 ".snp.b38"}' \
    "${PREFIX_BESD}.esi.b38"

# remove chromosome not in 1-22
rm ${OUTDIR}/bulk_wb/chr*_*.snp.b38

GENO_DIR="resources/genotypes/tenk10k_phase1_common"
# Split besd into chromosomes
for i in {1..22}; do
    # extract snp list available in b38 and in tenk10k
    awk -v OFS='\t' \
        'NR == FNR {a[$4] = $2; next}
         $2 in a {print $4, a[$2], $2}' \
        "${GENO_DIR}/chr${i}.bim" "${PREFIX_BESD}.esi.b38" \
        > "${OUTDIR}/bulk_wb/chr${i}.snp.b38"
    
    smr --beqtl-summary "${PREFIX_BESD}" \
        --probe-chr ${i} \
        --snp-chr ${i} \
        --extract-snp <(cut -f1 "${OUTDIR}/bulk_wb/chr${i}.snp.b38") \
        --thread-num ${THREADS} \
        --make-besd \
        --out "${OUTDIR}/bulk_wb/chr${i}" \
        2>&1 | tee "${snakemake_log[0]}"
    
    # update esi to b38 and tenk10k id
    cp "${OUTDIR}/bulk_wb/chr${i}.esi" "${OUTDIR}/bulk_wb/chr${i}.bak.esi"

    awk -v OFS='\t' \
        'NR==FNR {a[$1]=$2; b[$1]=$3; next}
         $2 in a {x = $2; $2 = a[x]; $4 = b[x];  print}' \
        "${OUTDIR}/bulk_wb/chr${i}.snp.b38" "${OUTDIR}/bulk_wb/chr${i}.bak.esi" \
        > "${OUTDIR}/bulk_wb/chr${i}.esi"

done
