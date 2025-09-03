library(tidyverse)
library(data.table)
library(patchwork)
library(paletteer)

# create heatmap with Helvetica font 
new_annot <- data.table(
  old = c("Canonical or Drug Pathway", "Differentially Expressed in Crohn's Disease"),
  new = c("Canonical or Drug Pathway", "Differentially Expressed in Crohn's Disease")
)
load("resources/crohns_case_study/heatmap_objects.RData")

# plot_data <- readRDS("resources/crohns_case_study/plot_data.RDS")
canon_pathway_deg_genes <- c("ENSG00000100365", "ENSG00000245532", "ENSG00000206503",
                             "ENSG00000188906", "ENSG00000115232", "ENSG00000096968")

mark <- "*"
setDT(plot_data)
plot_data[new_annot, annot := i.new, on = .(annot = old)]
plot_data[probeID %in% canon_pathway_deg_genes, Gene := paste0(mark, Gene)]

plot_data[, annot := factor(annot, levels = c(new_annot$new))]

# plot_data[, .(n_gene = n_distinct(Gene)), by = annot]

source("scripts/preprocess.R")
source("scripts/crohns_publication/locus_zoom.R")

# max_ptransform <- df_msmr[phenotype == "crohns", max(-log10(p_SMR_multi))]

# df_crohns_summary <- df_msmr %>% 
#   filter(phenotype == "crohns") %>% 
#   group_by(probeID) %>% 
#   summarise(n_celltypes = n_distinct(cell_type),
#             celltypes = list(cell_type))

# get celltypes with n_celltype == 1
# df_celltype_specific <- df_msmr %>% 
#   filter(phenotype == "crohns") %>% 
#   group_by(probeID) %>% 
#   filter(n() == 1) %>% 
#   group_by(cell_type) %>% 
#   tally(name = "n_specific_mr_genes")
# 
# df_celltype_genes <- df_msmr %>% 
#   filter(phenotype == "crohns") %>%
#   group_by(cell_type) %>% 
#   tally(name = "n_mr_genes")
# 
# df_celltype_summary <- df_celltype_genes %>% 
#   left_join(df_celltype_specific, by = "cell_type") %>% 
#   mutate(n_specific_mr_genes = replace_na(n_specific_mr_genes, 0),
#            prop = n_specific_mr_genes / n_mr_genes)

max_ptransform <- max(abs(plot_data$p_transform), na.rm = TRUE)

plot_data[df_cell_map, major_cell_type := i.major_cell_type, on = .(cell_type)]
plot_data[, fct_gene := fct_reorder(Gene, sig, ~sum(.x, na.rm = T), .na_rm = TRUE)]

(p_heatmap <- ggplot(plot_data, aes(y = fct_gene, x = cell_type)) +
  facet_grid(rows = vars(annot), 
             scales = "free", space = "free") +
  geom_blank() +
  # fake geom_point to set the legend
  geom_point(aes(shape = "MR evidence"),
             data = ~inner_join(.x, df_targets %>% select(probeID = probe, biosample)),
             size = 1, stroke = 1, color = "black") +
  # geom_tile(fill = "grey95", colour = "white", linewidth = 0.25) +
  geom_tile(aes(fill = p_transform), linewidth = 0.25, color = "gray", data = ~filter(.x, !is.na(p_transform), sig == FALSE)) +
  geom_tile(aes(fill = p_transform), linewidth = 0.5, color = "black",
            data = ~filter(.x, !is.na(p_transform), sig == TRUE)) +
  geom_point(aes(shape = "Example"),
            data = ~inner_join(.x, df_targets %>% select(probeID = probe, biosample)),
            size = 3, color = "#6A1B9A", stroke = 1.5) +
  scale_shape_manual(values = c("MR evidence" = 22, "Example" = 1),
                     guide = guide_legend(direction = "vertical", override.aes = list(size= 2.5))) +
  # geom_tile(data = ~filter(.x, !is.na(p_transform)), color = "black", linewidth = 0.5) +
  # geom_point(aes(shape = "MR Gene"), size = 1, 
  #            data = ~filter(.x, sig))+
  paletteer::scale_fill_paletteer_c(
    "ggthemes::Red-Blue-White Diverging",
    na.value = "grey90", limits = c(-max_ptransform, max_ptransform), direction = -1,
    guide = guide_colorbar(theme = theme(legend.key.width = unit(7.5, "lines"),
                                        legend.key.height = unit(0.75, "lines")))) +
  scale_x_discrete() +
  labs(colour = NULL, x = NULL, fill = bquote(-log[10]~italic(P) %*% "direction of effect"), size = 5, shape = NULL) +
  theme_minimal() +
  coord_fixed(clip = "off") +
  geom_text(aes(label = paste0(fct_gene, " "), x = -Inf),
            data = ~.x %>% filter(probeID %in% unique(df_targets$probe)) %>% 
              distinct(fct_gene, annot),
            size = 8/.pt, color = "red3", fontface = "bold.italic", hjust = 1) +
  geom_text(aes(label = paste0(fct_gene," "), x = -Inf),
            data = ~.x %>% filter(!probeID %in% unique(df_targets$probe) & !is.na(probeID)) %>%
              distinct(fct_gene, annot),
            size = 7/.pt, color = "gray30", fontface = "italic", hjust = 1) +
  theme(axis.ticks = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 7, color = NA, vjust = 0.5, hjust = 1, face = "italic"),
        plot.title = element_text(size = 50, face = "bold"),
        legend.position = "bottom",
        legend.title.position = "top",
        strip.clip = "off",
        # strip.text.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.spacing.x = unit(1, "lines"),
        plot.margin = margin(),
        legend.text = element_text(size = 7, margin = margin(t = 1)),
        legend.title = element_text(size = 7, vjust = 0.75)) 
)


# panel c: DEG umap label
p_umap_annot <- png::readPNG("figures/crohns/deg_umap_annotation_only.png", native =TRUE) %>% 
  wrap_elements(full = .)

# panel d: DEG umap label
p_umap_deg <- png::readPNG("figures/crohns/deg_umap_genes_only_merged_colourbar.png", native =TRUE) %>% 
  wrap_elements(full = .)

panel_c <- (p_umap_annot | p_umap_deg) +
  plot_layout(ncol = 2, guides = "keep",
              widths = c(1.4, 2))

panel_ab <- (wrap_elements(full = p_heatmap) | wrap_elements(full = mr_example_plots)) +
  plot_layout(nrow = 1, width = c(1, 1.8))

(plots <- (panel_ab / panel_c) +
  plot_layout(ncol = 1, height = c(3.3, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 14, face = "bold"),
        plot.tag.position = c(0, 1))
)

# ggsave("figures/crohns/p_heatmap_annot_helvetica.png", p_heatmap, device = ragg::agg_png(),
#        width = 5.0, height = 15, bg = "white", scaling = 1.5, dpi = 300)


ggsave("figures/crohns/crohns_combined_v2.png", plots, device = ragg::agg_png(),
       width = 9, height = 12.5, bg = "white", scaling = 1, dpi = 300)

ggsave("figures/publication_pdf/Fig5 - Crohn's disease.pdf", plots, device = cairo_pdf,
       width = 9, height = 12.5, bg = "white", scale = 1, dpi = 320)
 
# (p_n_celltype <- plot_data[, .(n = sum(sig, na.rm = TRUE)), by = list(fct_gene, annot)] %>% 
#   ggplot(aes(y = fct_gene, x = n)) +
#   geom_col(aes(fill = n), width = 0.5, show.legend = FALSE) +
#   # geom_tile(aes(fill = n), color = "white") +
#   # scale_fill_viridis_c(direction = -1, na.value = "grey90") +
#   scale_fill_stepsn(colours = paletteer_d("MexBrewer::Frida")) +
#   facet_grid(rows = vars(annot), scales = "free", space = "free") +
#   # geom_text(aes(label = n), size = 2.5, color = "black") +
#   theme_classic() +
#   labs(x = "Number of cell types\nwith MR evidence", y = NULL) +
#   layer_scales(p_heatmap)$y +
#   scale_x_continuous(
#     breaks = c(1,5,10,15,20),
#     expand = expansion(mult = c(0, 0.05))) +
#   theme(strip.text.y = element_text(angle = -90),
#         strip.background = element_blank(),
#         axis.line.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.text.y = element_blank(),
#         plot.margin = margin()
#         )
# )

# (p_heatmap + p_n_celltype + 
#   plot_layout(nrow = 1, widths = c(5, 1))
# )
# 
# ggsave("causal_inference_manuscript/figs/main_corrected/heat_facet_annot_2_categories_Roboto.png", FIG1, device = ragg::agg_png(),
#        width = 5.0, height = 16, bg = "white", scaling = 1.2, dpi = 300)

