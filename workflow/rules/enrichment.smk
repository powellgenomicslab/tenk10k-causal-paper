...existing code...

rule gget_enrichr:
    """
    Perform gene set enrichment analysis using Enrichr via gget.
    Input: Gene universe, MSMR significant genes, gene set list
    Output: Enrichment results per phenotype and biosample
    Note: Requires internet access to connect to external databases.
    """
    input:
        gene_universe = "results/enrichment/{study}/gene_universe.txt",
        msmr_sig = "results/aggregate/msmr_sig/{study}~q_{q_thresh}~heidi_{heidi_thresh}.msmr_sig.tsv",        
        # full gene_set list: https://maayanlab.cloud/Enrichr/#libraries
        gene_set = "resources/misc/enrichr.gene_set.txt"
    conda:
        "pydata"
    output:
        phenotype = "results/enrichment/{study}~q_{q_thresh}~heidi_{heidi_thresh}.enrichr.phenotype.tsv",
        biosample = "results/enrichment/{study}~q_{q_thresh}~heidi_{heidi_thresh}.enrichr.biosample.tsv"
    script: "snakescripts/enrichment/gget_enrichr.py"

checkpoint gget_enrichr_dir_pheno:
    """
    Prepare input files for phenotype-specific Enrichr enrichment analysis using gget.
    Input: MSMR significant genes
    Output: Directory of gene-celltype files per phenotype
    """
    input: msmr_sig = "results/aggregate/msmr_sig/{study}~q_{q_thresh}~heidi_{heidi_thresh}.msmr_sig.tsv"
    log: "logs/enrichment/gget_enrichr_dir_pheno/{study}~q_{q_thresh}~heidi_{heidi_thresh}.log"
    output:
        dir_pheno = directory("results/enrichment_pheno/{study}~q_{q_thresh}~heidi_{heidi_thresh}/gene_celltype")
    shell:
        """
        mkdir -p {output.dir_pheno}
        cut -f1,2,5 {input.msmr_sig} | \
            sort -u | \
            awk 'NR > 1 {{print $1,$3 > "{output.dir_pheno}/" $2 ".txt"}}'
        """

rule gget_enrichr_pheno:
    """
    Perform phenotype-specific gene set enrichment analysis using Enrichr via gget.
    Input: Gene universe, MSMR significant genes per phenotype, gene set list
    Output: Enrichment results per phenotype
    Note: Requires internet access to connect to external databases.
    """
    input:
        gene_universe = "results/enrichment/{study}/gene_universe.txt",
        msmr_sig = "results/enrichment_pheno/{study}~q_{q_thresh}~heidi_{heidi_thresh}/gene_celltype/{pheno}.txt",        
        # full gene_set list: https://maayanlab.cloud/Enrichr/#libraries
        gene_set = "resources/misc/enrichr.gene_set.txt"
    conda: "pydata"
    log: "logs/enrichment/gget_enrichr_pheno/{study}~q_{q_thresh}~heidi_{heidi_thresh}/{pheno}.log"
    params: min_gene = 5
    resources:
        queue = "copyq",
        ncpus = "1",
        mem = "4GB",
        time = "06:00:00"
    output:
        pheno = "results/enrichment_pheno/{study}~q_{q_thresh}~heidi_{heidi_thresh}/enrichr/{pheno}.tsv"
    script: "snakescripts/enrichment/gget_enrichr_pheno.py"

def gget_enrichr_pheno_aggregate(x):
    DIR = Path(checkpoints.gget_enrichr_dir_pheno.get(**x).output[0])
    pheno = [f.with_suffix('').name for f in DIR.glob("*.txt")]
    return [f"results/enrichment_pheno/{x.study}~q_{x.q_thresh}~heidi_{x.heidi_thresh}/enrichr/{p}.tsv"
            for p in pheno]
                  
rule aggregate_gget_enrichr:
    """Aggregate Enrichr enrichment analysis per phenotype using gget
    note: this requires internet access to connect to external databases
    """
    input: gget_enrichr_pheno_aggregate
    output: "results/aggregate/{study}/enrichr.q_{q_thresh}~heidi_{heidi_thresh}.tsv.gz"
    log: "logs/aggregate/enrichr.{study}.q_{q_thresh}.heidi_{heidi_thresh}.log"
    shell: "awk 'NR == FNR || FNR > 1' {input} | gzip -c > {output}"
