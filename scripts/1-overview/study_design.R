# get numbers for manuscript writing

source("scripts/0-preprocess/preprocess_results.R")

library(patchwork)
library(paletteer)
library(geomtextpath)

# n genes by cell type
df_gene <- df_msmr_tenk10k %>%
  distinct(cell_type, major_cell_type, probeID, Gene, gene_type) %>%
  mutate(gene_cat = case_when(
    gene_type == "protein_coding" ~ "Protein-coding",
    TRUE ~ "Non-coding"
  ) %>% factor(levels = c("Protein-coding", "Non-coding")))

df_overall <- df_gene %>%
  distinct(probeID, gene_cat) %>%
  mutate(
    cell_type = "Overall",
    major_cell_type = "Overall"
  )

pals <- c("#D6D8D0FF", "#A4ABB0FF", "#4C6C94FF", "#435E7FFF", "#2F415FFF", "#232C43FF", "#0B1829FF")
(p_gene <- bind_rows(df_gene, df_overall) %>%
  mutate(
    cell_type = factor(cell_type, c("Overall", levels(df_gene$cell_type))),
    major_cell_type = factor(major_cell_type, c("Overall", levels(df_gene$major_cell_type)))
  ) %>%
  ggplot(aes(x = cell_type, group = fct_rev(gene_cat))) +
  theme_classic() +
  geom_bar(aes(fill = fct_rev(gene_cat))) +
  facet_grid(cols = vars(major_cell_type), scales = "free_x", space = "free_x") +
  labs(y = "Number of eGenes", x = NULL) +
  scale_fill_manual(values = pals[c(2, 4)], name = "Gene type") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  coord_cartesian(clip = "off") +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_text(size = 8, face = "bold"),
    axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1, size = 7),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 8),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 9),
    strip.clip = "off",
    legend.key.size = unit(0.75, "lines"),
    legend.key.spacing.y = unit(0.25, "lines"),
    panel.grid.major.y = element_line(color = "gray90"),
    legend.position = "right"
  )
)

df_traits <- df_trait_map %>%
  left_join(df_magma[, .(n_magma = .N), by = phenotype],
    by = c(trait_id = "phenotype")
  ) %>%
  mutate(pheno_cat = factor(cat_rev, cat_order)) %>%
  mutate(pheno_label = fct_reorder(label, n_magma), .by = pheno_cat) %>%
  mutate(scaled_n_eff = scales::rescale(n_eff, to = c(7e3, 1e4)))

(p_traits <- df_traits %>%
  ggplot(aes(x = pheno_label, fill = pheno_cat)) +
  theme_classic() +
  # geom_point() +
  geom_segment(aes(xend = pheno_label, y = n_magma, yend = -Inf),
    color = "gray90", linetype = "solid"
  ) +
  geom_textline(aes(y = n_magma, label = pheno_cat, group = pheno_cat),
    hjust = 0.4, vjust = -0.8, size = 9 / .pt, fontface = "bold"
  ) +
  geom_point(aes(y = n_magma, size = n_eff),
    shape = "circle filled"
  ) +
  labs(y = "Number of GWAS genes", x = NULL) +
  coord_cartesian(clip = "off") +
  # geom_point(aes(y = scaled_n_eff)) +
  facet_grid(cols = vars(pheno_cat), space = "free", scales = "free") +
  scale_fill_manual(values = deframe(df_trait_cat[, .(cat_order, color)]), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_x_discrete(expand = expansion(add = 1)) +
  scale_size_continuous(
    labels = scales::number_format(scale = 1e-3, suffix = "k"),
    name = "Effective GWAS sample size",
    guide = guide_legend(nrow = 1, override.aes = list(
      fill = "gray", color = "black"
    ))
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 8),
    strip.text = element_blank(),
    strip.background = element_blank(),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8, margin = margin(l = 0)),
    strip.clip = "off",
    panel.grid.major.y = element_line(color = "gray90"),
    legend.position = "inside",
    plot.margin = margin(t = 2, unit = "lines"),
    legend.position.inside = c(1, 1),
    legend.justification = c(1, 0.5)
  )
)

infographic <- png::readPNG("figures/biorender/study_design.png", native = TRUE)

plots <- (wrap_elements(full = infographic) /
  wrap_elements(full = p_gene) /
  wrap_elements(full = p_traits)) +
  plot_layout(heights = c(0.6, 0.2, 0.45)) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(face = "bold"),
    plot.tag.position = c(0, 0.97),
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave(plots,
  filename = "figures/main/1-study_design.pdf",
  width = 9.5 / 1.05, height = 12 / 1.05, units = "in", dpi = 320,
  device = cairo_pdf
)

# write supplementary table
source("scripts/util/write_table.R")

df_celltype <- df_gene %>%
  group_by(gene_cat, cell_type) %>%
  tally() %>%
  ungroup() %>%
  pivot_wider(names_from = gene_cat, values_from = n, values_fill = 0) %>%
  mutate(n_gene = `Protein-coding` + `Non-coding`) %>%
  left_join(df_cell_map, by = "cell_type")

write_table(df_celltype, "celltype_summary", 1)

df_traits %>%
  mutate(across(starts_with("n_"), as.numeric)) %>%
  write_table("trait_summary", 2)
