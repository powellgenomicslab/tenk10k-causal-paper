# Prepare Crohn's DEG for downstream analysis 

library(tidyverse)

crohns_dir <- "resources/crohns_case_study"

# Find sig deg
#deg2 <- read.csv("crohns_disease_deg/data/publication_data/Kong_2023_Immunity_Crohnssc/Kongetal_DEGenes.csv")

deg <- readxl::read_xlsx(paste0(crohns_dir, "/deg/Kongetal2023_supplementary/Kongetal2023_Crohns-sc_DEG.xlsx"), sheet = 1)

colnames(deg)
# remove random cols
deg <- deg[,1:10]
# get immune cells only
deg <- deg %>% rename(scRNAseq_cellid = `Cell subset`)
deg <- deg %>% filter(grepl(pattern = "Immune", scRNAseq_cellid))
deg$scRNAseq_cellid <- gsub("Immune.", "", as.character(deg$scRNAseq_cellid))

deg <- deg %>% select(-`Continuous DE coefficients`, -`Continuous DE coefficients p value`, -`Continuous FDR` )

# get the transformed p values
deg$disc_p_transform <- -log10(deg$`Discrete DE coefficients p value`)*sign(deg$`Discrete DE coefficients`) 


# Renaming the contrast 
deg <- deg %>% mutate(Location = ifelse(Location == "CO", "Colon", "Terminal Ileum"), Contrast = ifelse(Contrast == "Infl vs. Heal", "Inflamed Tissue vs Healthy", "Non-Inflamed Tissue vs Healthy"))
deg$Location <- as.factor(deg$Location)
deg$Contrast <- as.factor(deg$Contrast)

#add the "tenk10k harmonised" cell types - basically their cell type groups 

#first export it, then resave it after adding 
#write.csv(sort(unique(deg$scRNAseq_cellid)), 'crohns_disease_deg/metadata/deg_celltype_groups_raw.csv', row.names = F)
# the major cell types in tenk are: 

cell_map <- read.delim("resources/metadata/cell_map.tsv")
major_cell_types_tenk = unique(cell_map$major_cell_type)
# [1] "CD4 T"            "Unconventional T" "CD8 T"            "NK"               "B"                "Monocyte"         "Dendritic"        "HSPC"  

# manually add the major cell types, then read in the updated file

deg_cell_map <- read.csv(paste0(crohns_dir, "/deg/cell_annotation_or_features/deg_celltype_groups.csv"))
deg$major_cell_type <- deg_cell_map$major_cell_type[match(deg$scRNAseq_cellid, deg_cell_map$scRNAseq_cellid)]

deg$major_cell_type <- as.factor(deg$major_cell_type)
deg$scRNAseq_cellid <- as.factor(deg$scRNAseq_cellid)


# create a text file of the common cell types maybe useful later.
# intersect major cell types from tenk and deg
common_cell_types <- intersect(cell_map$major_cell_type, deg_cell_map$major_cell_type)
write.table(common_cell_types, paste0(crohns_dir, "/deg/common_cell_types_with_tenk.txt"), sep = "\t", col.names = F, quote = F, row.names = F)


# filter to significant only and get rid of cell types not in tenk10k and sign of both discrete and continous are matched
deg <- deg %>% filter(`Discrete FDR` < 0.05, major_cell_type %in% common_cell_types, Location == "Colon")

saveRDS(deg, paste0(crohns_dir, "/deg/crohns_deg_pre-processed.RDS"))



