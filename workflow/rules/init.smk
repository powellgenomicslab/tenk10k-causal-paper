rule mk_gene_annot:
  """Make gene type annotation from gencode gtf"""
  input: "resources/misc/gencode.v44.basic.annotation.gtf"
  output: "resources/misc/gencode.{version}.gene_type.tsv"
  conda: "renv"
  script: "snakescripts/mk_gene_annot.R"

