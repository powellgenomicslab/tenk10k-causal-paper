rule prep_smr_input:
    """
    Prepare SMR input files (p-value threshold, probe list)
    Input: R script for SMR input preparation
    Output: SMR input directory
    """
    input: "workflow/rules/snakescripts/prep_smr_input/{study}.R"
    output: directory("resources/smr/{study}/")
    conda: "renv"
    script: "snakescripts/prep_smr_input/{wildcards.study}.R"

rule prep_besd_chr:
    """
    Prepare BESD files for SMR per chromosome
    Input: Shell script for BESD preparation
    Output: BESD directory
    """
    input: "workflow/rules/snakescripts/prep_besd_chr/{study}.sh"
    output: directory("resources/besd/{study}/")
    threads: 8
    conda: "pydata"
    log: "logs/prep_besd_chr/{study}.log"
    script: "snakescripts/prep_besd_chr/{wildcards.study}.sh"

rule run_smr:
    """
    Run SMR per chromosome, cell type, and phenotype.
    Input: BESD, EPI, ESI, p-threshold, probe, genotype files, GWAS summary.
    Output: SMR results and failed SNP list.
    """
    input:
        besd = expand("resources/besd/{{study}}/{{biosample}}/chr{chr}.{ext}", chr=range(1, 23), ext=["besd", "epi", "esi"]),
        pthresh = "resources/smr/{study}/{biosample}/pthresh.txt",
        probe = "resources/smr/{study}/{biosample}/probe.txt",
        bfile = expand("resources/genotypes/{{study}}_common/chr{chr}.{ext}", chr=range(1, 23), ext=["bed", "bim", "fam"]),
        ma = "resources/ma/{pheno}.ma"
    output:
        smr = "results/smr/{study}/{biosample}/{pheno}/all_chr.msmr",
        fail = "results/smr/{study}/{biosample}/{pheno}/all_chr.snp_failed_freq_ck.list"
    threads: 8
    log: "logs/smr/{study}/{biosample}/{pheno}.log"
    resources:
        mem="16G"
    params:
        prefix_bfile_chr = "resources/genotypes/{study}_common/chr",
        prefix_besd_chr = "resources/besd/{study}/{biosample}/chr",
        prefix_out_chr = "results/smr/{study}/{biosample}/{pheno}/chr",
        maf = 0.01,
        diff_freq_prop = 0.2
    script: "snakescripts/run_smr.sh"

def target_smr(x):
    """
    Target function to run SMR for all biosamples and phenotypes
    """
    BIOSAMPLES = [d.name for d in Path(f"resources/smr/{x.study}").iterdir() if d.is_dir()]
    PHENOS = [f.with_suffix('').name \
              for f in list(Path("resources/ma/").glob("*.ma"))]
    TARGET = expand(f"results/smr/{x.study}/{{biosample}}/{{pheno}}/all_chr.msmr", \
                     biosample = BIOSAMPLES, \
                    pheno = PHENOS)
    return TARGET

rule run_smr_all:
    """
    Run SMR for all biosamples and phenotypes
    """
    input: target_smr
    output: touch("results/smr/{study}/.done")
    shell: "touch {output}"


rule concat_smr_all:
    """
    Run SMR for all biosamples and phenotypes
    """
    input: target_smr
    output: "results/aggregate/{study}.msmr.gz.parquet"
    conda: "renv"
    log: "logs/aggregate/{study}.msmr.log"
    resources:
        mem="64G",
        ncpus=16
    script: "snakescripts/aggregate/concat_smr_all.R"


rule smr_to_parquet:
    """
    Convert SMR output to parquet format
    """
    input: "TenK10K_SMR_brenner/results/combined/smr_all_traits_with_egene_specificity.csv"
    output: "results/smr/smr_combined_specificity.gz.parquet"
    conda: "renv"
    script: "snakescripts/smr_to_parquet.R"

rule smr_prepare_genelist:
    """
    Prepare gene list from GTF annotation for SMR analysis.
    Input: GTF file
    Output: Gene list text file
    """
    input: "resources/smr_misc/{study}.gtf.gz"
    output: "resources/smr_misc/{study}.genelist.txt"
    conda: "renv"
    script: "snakescripts/smr/prep_genelist.R"

rule smr_extract_locus:
    """
    extract SMR locus information
    see: https://yanglab.westlake.edu.cn/software/smr/#SMRlocusplot19
    """
    input:
        besd = expand("resources/besd/{{study}}/{{biosample}}/chr{chr}.{ext}",
                      chr = range(1, 23), ext = ["besd", "epi", "esi"]),  
        genelist = "resources/smr_misc/{study}.genelist.txt",
        bfile = expand("resources/genotypes/{{study}}_common/chr{chr}.bed", \
                       chr = range(1, 23), ext = ["bed", "bim", "fam"]),
        ma = "resources/ma/{pheno}.ma"
    output:
        "results/smr_locus/{study}/{biosample}/{pheno}/plot/locus.{probe}.txt"
    threads: 8
    log: "logs/smr_locus/{study}/{biosample}/{pheno}.{probe}.log"
    params:
        prefix_bfile_chr = "resources/genotypes/{study}_common/chr",
        prefix_besd_chr = "resources/besd/{study}/{biosample}/chr",
        prefix_out = "results/smr_locus/{study}/{biosample}/{pheno}/locus",
        probe_wind = 500
    script:
        "snakescripts/smr_locus.sh"

rule get_gene_universe:
    """
    Generate gene universe by intersecting MAGMA and MSMR gene sets.
    Input: MAGMA and MSMR gene-level parquet files
    Output: Gene universe text file
    """
    input:
        msmr = "results/aggregate/{study}.msmr.gz.parquet",
        magma = "results/aggregate/{study}.magma.gz.parquet"
    output:
        gene_universe = "results/enrichment/{study}/gene_universe.txt"
    conda: "renv"
    script: "snakescripts/enrichment/get_gene_universe.R"

rule msmr_sig:
    """
    Get significant SMR results,
    filtered to gene_universe (based on intersection between MAGMA and MR tested genes)
    """
    input:
        msmr = "results/aggregate/{study}.msmr.gz.parquet",
        gene_universe = "results/enrichment/{study}/gene_universe.txt"
    output:
        msmr_sig = "results/aggregate/msmr_sig/{study}~q_{q_thresh}~heidi_{heidi_thresh}.msmr_sig.tsv"
    conda: "renv"
    script: "snakescripts/aggregate/msmr_sig.R"

rule locus_zoom_extract:
    """
    extract SMR locus information
    for manual locus zoom plot
    """
    input:
        saige = "resources/saige_eqtl/{study}/{biosample}/common_raw.tsv",
        gtf = "resources/smr_misc/{study}.gtf.gz",
        ma = "resources/ma/{pheno}.ma",
        dir_geno_chr = "resources/genotypes/{study}",
        dir_besd_chr = "resources/besd/{study}/{biosample}"
    output:
        gwas = "results/smr_locus/{study}/{biosample}/{pheno}/{probe}.gwas.tsv",
        eqtl = "results/smr_locus/{study}/{biosample}/{pheno}/{probe}.eqtl.tsv",
        gtf = "results/smr_locus/{study}/{biosample}/{pheno}/{probe}.gtf",
        ld = "results/smr_locus/{study}/{biosample}/{pheno}/{probe}.ld"
    threads: 8
    log: "logs/smr_locus/{study}/{biosample}/{pheno}.{probe}.locus_zoom.log"
    conda: "renv"
    params:
        probe_flank_kb = 500
    script:
        "snakescripts/smr_locus/extract_region.R"