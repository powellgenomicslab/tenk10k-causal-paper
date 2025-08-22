...existing code...

rule format_gwas:
    """
    Format GWAS summary statistics for a given phenotype.
    Input: Raw GWAS file, trait metadata, chain file, liftover script
    Output: Formatted .ma file
    """
    input:
        gwas = "resources/sumstats/gwas/{pheno}.gwas",
        trait_metadata = "resources/metadata/trait_metadata_curated.xlsx",
        hg19tohg38 = "resources/misc/hg19ToHg38.over.chain",
        liftover_script = "workflow/rules/snakescripts/hg19tohg38.R"
    output: "resources/pipeline_ma/{pheno}.ma"
    threads: 8
    resources:
        mem = "16G",
        jobfs = "8G",
        ncpus = 8
    conda: "renv"
    log: "logs/format_gwas/{pheno}.log"
    script: "snakescripts/format_gwas/{wildcards.pheno}.R"

rule mk_tabix_region:
    """
    Create tabix region file for genotype extraction.
    Input: Genotype BIM files
    Output: Tabix region file
    """
    input: expand("resources/genotypes/{{study}}_common/chr{chr}.bim", chr = range(1,23))
    output: "resources/misc/tabix_region/{study}_common.region"
    threads: 8
    shell: "cut -f1,4 {input} | sort -nk1 -nk2 > {output}"

def get_finngen_gwas(x):
    with open("workflow/config/finngen_meta_path.yaml", 'r') as f:
        config = yaml.full_load(f)
    base_dir = Path(config['base_dir'])
    return {
        "gwas": str(base_dir / f"{config['pheno_prefix'][x.pheno]}_meta_out.tsv.gz"),
        "tbi": str(base_dir / f"{config['pheno_prefix'][x.pheno]}_meta_out.tsv.gz.tbi")
    }

rule extract_finngen_gwas:
    """
    Extract Finngen GWAS data for a given phenotype and region.
    Input: Finngen GWAS file, tabix region file
    Output: Extracted GWAS TSV file
    """
    input:
        unpack(get_finngen_gwas),
        region = "resources/misc/tabix_region/{study}_common.region"
    output: temp("resources/sumstats/finngen_gwas_extract/{study}/{pheno}.tsv.gz")
    threads: 8
    envmodules: "htslib"
    shell:
        """
        tabix -hR {input.region} {input.gwas} | \
            gzip -c > {output}
        """

# rule mk_tabix_region:
#     """
#     Make tabix region file based on genotype file for extraction
#     """
#     input: expand("resources/genotypes/{{study}}_common/chr{chr}.bim", chr = range(1,23))
#     output: "resources/misc/tabix_region/{{study}}_common.region"
#     threads: 8
#     envmodules: "htslib"
#     shell:

rule format_finngen_gwas:
    """
    Format GWAS summary statistics from Finggen meta-analysis (with UKB + MVP EUR)
    for use in SMR
    """
    input:
        unpack(get_finngen_gwas),
        gwas_subset = rules.extract_finngen_gwas.output,
        metadata = "resources/misc/meta_analysis_mvp_ukbb_FinnGen_R12_MVP_UKBB_mapping.tsv",
        frq = "resources/genotypes_frq/{study}.frq",
    output: "resources/sumstats/finngen_gwas_extract/{study}/{pheno}.ma"
    threads: 8
    conda: "renv"
    log: "logs/format_finngen_gwas/{study}/{pheno}.log"
    script: "snakescripts/format_finngen_gwas.R"
