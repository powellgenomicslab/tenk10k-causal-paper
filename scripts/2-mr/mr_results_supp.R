# smr overview plot
source("scripts/0-preprocess/preprocess_results.R")

library(scales)
library(patchwork)
library(ragg)
library(paletteer)
library(geomtextpath)

df_intersect <- df_msmr %>%
    rename(gene = probeID) %>%
    group_by(gene, phenotype, pheno_label, pheno_cat) %>%
    tally(name = "n_celltypes") %>%
    full_join(df_magma %>% select(gene = GENE, phenotype, qval)) %>%
    mutate(
        n_celltypes = replace_na(n_celltypes, 0),
        gwas_gene = !is.na(qval)
    )

df_intersect_overall <- df_intersect %>%
    group_by(n_celltypes, gwas_gene) %>%
    tally() %>%
    mutate(prop = n / sum(n))

# quantify number of genes that are not in MAGMA
df_smr_gwas_by_pheno <- df_intersect %>%
    group_by(phenotype, gwas_gene) %>%
    tally() %>%
    arrange(phenotype, gwas_gene) %>%
    mutate(prop = n / sum(n))

# df_smr_gwas_by_pheno %>%
#   group_by(gwas_gene) %>%
#   summarise(median = median(prop)) %>%
#   arrange(gwas_gene, desc(prop)) %>%
#   View()

n_genes <- df_msmr %>%
    group_by(phenotype, pheno_label, pheno_cat, supercategory, cell_type, major_cell_type) %>%
    tally() %>%
    ungroup() %>%
    complete(nesting(phenotype, pheno_label, pheno_cat, supercategory),
        nesting(cell_type, major_cell_type),
        fill = list(n = 0)
    )

# df_msmr[,.N, by = list(phenotype, cell_type)] %>% summary()

setDT(n_genes)

# df_msmr[, .N, by = list(sign(b_SMR))] %>%
#   mutate(prop = N / sum(N))

df_summary_gene_pheno <- df_msmr_tenk10k %>%
    group_by(pheno_label, pheno_cat, probeID) %>%
    summarise(mr = max(sig), gwas = max(magma_gene)) %>%
    summarise(
        n_mr = sum(mr),
        n_gwas = sum(gwas),
        prop_mr_gwas = sum(mr & gwas) / sum(gwas),
        prop_gwas_mr = sum(mr & gwas) / sum(mr),
        mr_gwas_ratio = sum(mr) / sum(gwas),
        intersect_mr_gwas = sum(mr & gwas),
        mr_only = sum(mr & !gwas),
        gwas_only = sum(gwas & !mr),
        union_mr_gwas = sum(mr | gwas),
        jaccard_mr_gwas = sum(mr & gwas) / sum(mr | gwas)
    ) %>%
    ungroup() %>%
    mutate(pheno_label = fct_reorder(pheno_label, n_mr))

setDT(df_summary_gene_pheno)

df_summary_gene_pheno[df_trait_map, n_eff_gwas := i.n_eff,
    on = c(pheno_label = "label")
]

# scatter plot
rho_gwas_mr <- cor.test(df_summary_gene_pheno$n_gwas, df_summary_gene_pheno$n_mr)

(p_cor_gwas_mr <- ggplot(
    df_summary_gene_pheno,
    aes(x = n_gwas, y = n_mr)
) +
    theme_bw() +
    geom_point(aes(color = pheno_cat), size = 3, alpha = 0.8) +
    geom_smooth(color = "gray20", linewidth = 0.5) +
    paletteer::scale_color_paletteer_d("MetBrewer::Hiroshige",
        name = "Phenotype category"
    ) +
    labs(
        x = "Number of GWAS genes",
        y = "Number of MR genes"
    ) +
    scale_x_continuous(expand = 0.01) +
    scale_y_continuous(expand = 0.01) +
    geom_abline(slope = 1, intercept = 0, color = "gray50", linewidth = 0.5, linetype = "dashed") +
    coord_equal() +
    theme(
        panel.grid.minor = element_blank(),
        legend.position = "inside",
        # aspect.ratio = 1,
        legend.position.inside = c(0.01, 0.99),
        legend.justification = c(0, 1)
    )
)

ggsave(p_cor_gwas_mr,
    filename = "figures/smr/smr_gwas_correlation.png",
    width = 6, height = 8, device = agg_png, scaling = 1.2
)

# get numbers
df_summary_gene_pheno[, list(
    n_trait = .N,
    n_mr_trait = sum(n_mr),
    avg_mr = mean(n_mr)
), ,
by = pheno_cat
]

# Comparison between MR genes - GWAS genes and TenK10K MR genes - eQTLgen

df_gene_summary_by_pheno <- df_msmr %>%
    rename(gene = probeID) %>%
    group_by(gene, phenotype, pheno_label, pheno_cat) %>%
    summarise(
        n_celltypes = n(),
        gwas = max(magma_gene, na.rm = TRUE) %>% as.logical(),
        eqtlgen = max(eqtlgen_mr, na.rm = TRUE) %>% as.logical()
    ) %>%
    pivot_longer(c(gwas, eqtlgen),
        names_to = "annotation", values_to = "value"
    ) %>%
    group_by(phenotype, pheno_label, pheno_cat, annotation, value) %>%
    tally(name = "n_genes") %>%
    ungroup() %>%
    complete(nesting(phenotype, pheno_label, pheno_cat, annotation),
        nesting(value),
        fill = list(n_genes = 0)
    ) %>%
    group_by(phenotype, pheno_label, pheno_cat, annotation) %>%
    mutate(
        n_mr_genes = sum(n_genes),
        prop_genes = n_genes / sum(n_genes)
    ) %>%
    mutate(fct_pheno = fct_reorder(pheno_label, n_mr_genes, mean) %>% fct_rev())

pals <- paletteer_d("MexBrewer::Frida")[c(3, 2)] %>% set_names(c("TRUE", "FALSE"))

mkplot_gene_summary <- function(df, x_col, y_col = fct_pheno,
                                xlab = NULL, col_scheme = pals,
                                labs = waiver(), name = NULL) {
    df %>%
        ggplot(aes(y = {{ y_col }}, x = {{ x_col }})) +
        theme_bw() +
        facet_wrap(~pheno_cat, nrow = 1, scale = "free_y", space = "free_y") +
        geom_col(aes(fill = value), width = 1, color = "black", linewidth = 0.4) +
        scale_fill_manual(
            values = pals, breaks = names(pals),
            labels = labs, name = NULL
        ) +
        labs(y = NULL, x = xlab) +
        scale_x_continuous(expand = expansion(0)) +
        theme(
            strip.text = element_text(
                hjust = 0, face = "bold", size = 9,
                margin = margin(b = 0.2, unit = "lines")
            ),
            strip.background = element_blank(),
            legend.key.size = unit(0.5, "lines"),
            legend.text = element_text(size = 9, margin = margin(l = 2)),
            legend.key.spacing.y = unit(0.25, "lines"),
            axis.title.x = element_text(size = 10),
            axis.text.x = element_text(size = 9),
            # panel.grid.major.y = element_blank(),
            panel.grid.minor = element_blank(),
            panel.spacing = unit(0.2, "lines")
        )
}


(p_gwas_n <- df_gene_summary_by_pheno %>%
    filter(annotation == "gwas") %>%
    mkplot_gene_summary(
        x_col = n_genes,
        xlab = "Number of unique MR genes",
        labs = c(
            "FALSE" = "MR only",
            "TRUE" = "MR & GWAS"
        )
    ) +
    theme(
        plot.margin = margin(),
        legend.position = "none"
    )
)

(p_gwas_prop <- df_gene_summary_by_pheno %>%
    filter(annotation == "gwas") %>%
    mkplot_gene_summary(
        x_col = prop_genes,
        xlab = "Proportion of MR genes",
        labs = c(
            "FALSE" = "MR only",
            "TRUE" = "MR & GWAS"
        )
    ) +
    scale_x_continuous(labels = percent, expand = expansion(0)) +
    geom_label(
        aes(
            label = ifelse(value == TRUE,
                paste0(" ", number(n_genes, 1, big.mark = ",")),
                paste0(number(n_genes, 1, big.mark = ","), " ")
            ),
            x = ifelse(value == TRUE, -Inf, Inf),
            hjust = ifelse(value == TRUE, 0, 1)
        ),
        linewidth = 0, label.padding = unit(0, "lines"),
        fill = alpha("white", 0.5),
        vjust = 0.5, size = 7 / .pt,
        position = position_nudge(x = 0.05)
    ) +
    theme(
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        strip.text = element_text(color = NA),
        axis.ticks.y = element_line(color = "gray"),
        legend.direction = "horizontal",
        legend.position = "inside",
        legend.position.inside = c(1, 1),
        legend.justification = c(1, 0),
        plot.margin = margin(r = 0.5, unit = "lines"),
        axis.ticks.length.y = unit(1, "lines")
    )
)

p_gwas_mr <- p_gwas_n + p_gwas_prop +
    plot_layout(nrow = 1, guides = "keep") &
    layer_scales(p_gwas_n)$y

ggsave("figures/mr_overview/gwas_mr_gene_summary.png",
    p_gwas_mr,
    scaling = 1.05,
    width = 8, height = 12, device = agg_png
)

# eQTLgen MR genes
(p_eqtlgen_n <- df_gene_summary_by_pheno %>%
    filter(annotation == "eqtlgen") %>%
    mkplot_gene_summary(
        x_col = n_genes,
        xlab = "Number of unique MR genes",
        labs = c(
            "FALSE" = "TenK10K MR only",
            "TRUE" = "TenK10K & eQTLgen MR"
        )
    ) +
    theme(
        plot.margin = margin(),
        legend.position = "none"
    )
)

(p_eqtlgen_prop <- df_gene_summary_by_pheno %>%
    filter(annotation == "eqtlgen") %>%
    mkplot_gene_summary(
        x_col = prop_genes,
        xlab = "Proportion of MR genes",
        labs = c(
            "FALSE" = "TenK10K MR only",
            "TRUE" = "TenK10K & eQTLgen MR"
        )
    ) +
    scale_x_continuous(labels = percent, expand = expansion(0)) +
    geom_label(
        aes(
            label = ifelse(value == TRUE,
                paste0(" ", number(n_genes, 1, big.mark = ",")),
                paste0(number(n_genes, 1, big.mark = ","), " ")
            ),
            x = ifelse(value == TRUE, -Inf, Inf),
            hjust = ifelse(value == TRUE, 0, 1)
        ),
        linewidth = 0, label.padding = unit(0, "lines"),
        fill = alpha("white", 0.5),
        vjust = 0.5, size = 7 / .pt,
        position = position_nudge(x = 0.05)
    ) +
    theme(
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        strip.text = element_text(color = NA),
        axis.ticks.y = element_line(color = "gray"),
        legend.direction = "horizontal",
        legend.position = "inside",
        legend.position.inside = c(1, 1),
        legend.justification = c(1, 0),
        plot.margin = margin(r = 0.5, unit = "lines"),
        axis.ticks.length.y = unit(1, "lines")
    )
)


(p_eqtlgen <- p_eqtlgen_n + p_eqtlgen_prop +
    plot_layout(nrow = 1, guides = "keep") &
    layer_scales(p_eqtlgen_n)$y
)

ggsave("figures/mr_overview/eqtlgen_mr_gene_summary.png",
    p_eqtlgen,
    scaling = 1.05,
    width = 8, height = 12, device = agg_png
)

# get number:

setDT(df_gene_summary_by_pheno)
tab_gene_summary_by_pheno <- df_gene_summary_by_pheno %>%
    pivot_wider(
        names_from = c(annotation, value), values_from = c(n_genes, prop_genes),
        names_sep = "."
    ) %>%
    setDT()

# write to supplementary table
source("scripts/util/helper.R")
write_gs(tab_gene_summary_by_pheno, "mr_gwas_eqtlgen", 4)


# OLD: Proportion MR - GWAS
# mkplot <- function(df, xlab = NULL) {
#   df %>%
#     ggplot(aes(y = pheno_label, x = x)) +
#     theme_bw() +
#     facet_wrap(~pheno_cat, nrow = 1, scale = "free_y", space = "free_y") +
#     geom_col(aes(fill = name), width = 1, color = "black", linewidth = 0.4) +
#     scale_fill_manual(
#       values = paletteer_d("ggthemes::Tableau_10")[1:3] %>%
#         set_names(c("gwas_only", "intersect_mr_gwas", "mr_only")),
#       breaks = c("mr_only", "intersect_mr_gwas", "gwas_only"),
#       labels = c("gwas_only" = "GWAS only",
#                  "intersect_mr_gwas" = "MR & GWAS",
#                  "mr_only" = "MR only"),
#       name = "Gene set") +
#     labs(y = NULL, x = xlab) +
#     scale_x_continuous(expand = expansion(0)) +
#     theme(strip.text = element_text(hjust = 0, face = "bold", size = 9,
#                                     margin = margin(b = 0.2, unit = "lines")),
#           strip.background = element_blank(),
#           # panel.grid.major.y = element_blank(),
#           panel.grid.minor = element_blank(),
#           panel.spacing = unit(0.2, "lines"))
# }
#
# (p1 <- df_summary_gene_pheno %>%
#   pivot_longer(c(mr_only, intersect_mr_gwas, gwas_only),
#                values_to = "x") %>%
#   mkplot(xlab = "Number of genes")
# )
# (p2 <- df_summary_gene_pheno %>%
#   pivot_longer(c(intersect_mr_gwas, mr_only)) %>%
#   mutate(x = value / sum(value), .by = pheno_label) %>%
#   mkplot(xlab = "Proportion of MR genes") +
#   scale_x_continuous(labels = percent, expand = expansion(0)) +
#   guides(fill = "none") +
#   theme(axis.text.y = element_blank(),
#         strip.text = element_text(color = NA),
#         axis.ticks.y = element_line(color = "gray90"),
#         axis.ticks.length.y = unit(2, "lines"))
# )
#
# (p3 <- df_summary_gene_pheno %>%
#   pivot_longer(c(intersect_mr_gwas, gwas_only)) %>%
#   mutate(x = value / sum(value), .by = pheno_label) %>%
#   mutate(name = fct_rev(name)) %>%
#   mkplot(xlab = "Proportion of GWAS genes") +
#   scale_x_continuous(labels = percent, expand = expansion(0)) +
#   guides(fill = "none") +
#   theme(axis.text.y = element_blank(),
#         strip.text = element_text(color = NA),
#         axis.ticks.y = element_line(color = "gray90"),
#         axis.ticks.length.y = unit(2, "lines"))
# )
#
# plots_mr_gwas <- (p1 + p2 + p3) +
#   plot_layout(nrow = 1, guides = "collect") &
#   layer_scales(p1)$y &
#   theme(plot.margin = margin(),
#         axis.title = element_text(size = 11))
#
# ggsave("figures/smr/smr_gwas_mr_overview.png",
#        plots_mr_gwas, scaling = 0.9,
#        width = 9, height = 12, device = agg_png)
#
# # get numbers
# df_summary_gene_pheno[prop_mr_gwas > prop_gwas_mr]
#
# df_summary_gene_pheno[, sum(mr_only)]
# df_summary_gene_pheno[, sum(gwas_only)]
# df_summary_gene_pheno[, sum(intersect_mr_gwas)]
#
# df_summary_gene_pheno[mr_only >0, c(list(.N), summary(mr_only) %>% as.list)]
# df_summary_gene_pheno[gwas_only>0, c(list(.N),summary(gwas_only))]
# df_summary_gene_pheno[intersect_mr_gwas>0, c(list(.N), summary(intersect_mr_gwas))]


# direction of effect
df_msmr[, .N, by = list(sign(b_SMR), pheno_cat)]

# cell type summary
df_cell_summary <- df_msmr_tenk10k[,
    .(
        n_genes = length(unique(probeID)),
        n_mr_genes = .SD[sig == TRUE, length(unique(probeID))],
        n_sig = sum(sig),
        n_tested = .N,
        prop_sig = sum(sig) / .N
    ),
    by = list(cell_type, major_cell_type)
]

(p_scatter_cell <- df_cell_summary %>%
    ggplot(aes(x = n_sig, y = n_tested)) +
    theme_bw() +
    geom_smooth(color = "gray20", linewidth = 0.5) +
    geom_point(aes(color = cell_type)) +
    scale_color_manual(
        values = deframe(df_cell_map[, .(cell_type, color)]),
        name = NULL
    ) +
    labs(
        x = "Number of significant MR associations",
        y = "Number of MR tests"
    ) +
    scale_x_continuous(labels = label_number(scale = 1e-3, suffix = "k")) +
    scale_y_continuous(labels = label_number(scale = 1e-3, suffix = "k")) +
    # geom_abline(slope = 1, intercept = 0, color = "gray50", linewidth = 0.5, linetype = "dashed") +
    theme(
        aspect.ratio = 1,
        legend.text = element_text(size = 9)
    )
)

ggsave("figures/smr/smr_celltype_scatter.png",
    p_scatter_cell,
    scaling = 1.5,
    width = 10, height = 6, device = agg_png
)

custom_summary <- function(dt) {
    c(
        summary(dt$n) %>% as.list(),
        list(n_trait = length(unique(dt$phenotype)))
    )
}

n_genes[, custom_summary(.SD), by = pheno_cat]

ggplot(n_genes, aes(
    y = fct_reorder(pheno_label, n, sum),
    x = n
)) +
    geom_col() +
    theme_bw() +
    facet_grid(
        rows = vars(pheno_cat),
        scale = "free",
        space = "free"
    )
# +
#   geom_text(aes(label = n), data = ~.x[, .(n = sum(n)),by = pheno_label],
#             hjust = 0)


n_genes[, summary(n) %>% as.list(), by = supercategory]
# decile
n_genes[, quantile(n, probs = seq(0, 1, by = 0.1)) %>% as.list(), by = supercategory]

n_genes[, decile := ntile(n, 10), by = supercategory]

n_genes[supercategory == "disease"] %>%
    ggplot() +
    geom_histogram(aes(x = n), fill = "gray", color = "black", bins = 50)

# by celltype
n_genes %>%
    group_by(cell_type, supercategory) %>%
    summarise(median_n = median(n)) %>%
    # pivot_wider(names_from = "supercategory", values_from = "median_n") %>%
    arrange(desc(median_n))

n_genes %>%
    group_by(phenotype, supercategory) %>%
    summarise(median_n = median(n)) %>%
    arrange(desc(median_n)) %>%
    print(n = 100)


df_tally_pheno <- df_msmr %>%
    group_by(cell_type, major_cell_type, pheno_label, pheno_cat) %>%
    tally() %>%
    group_by(pheno_label, pheno_cat) %>%
    mutate(
        prop = n / sum(n),
        total_n = sum(n)
    )

(p_tally_pheno <- df_tally_pheno %>%
    ggplot(aes(x = pheno_label, y = prop)) +
    theme_bw() +
    geom_col(aes(fill = fct_rev(cell_type))) +
    scale_fill_manual(
        values = deframe(df_cell_map[, .(cell_type, color)]),
        name = NULL,
        guide = guide_legend(reverse = TRUE, nrow = 2),
    ) +
    facet_grid(cols = vars(pheno_cat), space = "free", scale = "free") +
    scale_y_continuous(
        expand = expansion(0),
        labels = percent,
        guide = guide_axis(cap = "both")
    ) +
    labs(y = "Proportion of cell types", x = NULL) +
    theme(
        legend.position = "bottom",
        strip.background = element_blank(),
        axis.line.y = element_line(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1),
        strip.text = element_blank(),
        panel.border = element_blank(),
        plot.margin = margin(t = 0),
        strip.clip = "off"
    )
)

df_tally_total <- df_tally_pheno %>%
    group_by(pheno_label, pheno_cat) %>%
    summarise(total_n = sum(n))

# collapse celltypes
df_tally_total_gene <- df_msmr %>%
    group_by(pheno_label, pheno_cat) %>%
    # group_by(pheno_label, pheno_cat, n_celltypes) %>%
    # tally(n = "n_genes") %>%
    summarise(total_n = length(unique(probeID)))

(p_tally_total <- df_tally_total_gene %>%
    ggplot(aes(x = pheno_label, y = log10(total_n))) +
    theme_classic() +
    geom_col(aes(fill = pheno_cat), show.legend = FALSE) +
    scale_fill_paletteer_d("ggthemes::Tableau_10") +
    geom_text(aes(label = number(total_n, 1, big.mark = ",")),
        hjust = 0,
        angle = 90,
        position = position_nudge(y = 0.02), size = 3
    ) +
    facet_grid(cols = vars(pheno_cat), space = "free", scale = "free", switch = "x") +
    scale_y_continuous(
        guide = guide_axis(cap = "both"), limits = c(0, 4),
        expand = expansion(add = c(0.1, 0.1)), breaks = 0:5
    ) +
    labs(y = bquote(log[10] ~ "N unique eGenes"), x = NULL) +
    coord_cartesian(clip = "off") +
    theme(
        strip.background = element_blank(),
        strip.clip = "off",
        strip.text = element_text(face = "bold"),
        panel.grid.major.x = element_blank(),
        plot.margin = margin(b = 0, t = 10),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.minor.x = element_blank()
    )
)

plots <- (
    (p_tally_total / p_tally_pheno) &
        scale_x_discrete(expand = expansion(add = 1)) &
        theme(
            panel.spacing.x = unit(0.5, "lines"),
            panel.grid.major.y = element_line(),
            legend.position = "bottom",
            legend.key.size = unit(0.8, "lines"),
            legend.text = element_text(size = 8)
        )
) +
    plot_layout(guides = "collect")

ggsave("figures/smr_tally_celltype_heidi.png",
    plots,
    scaling = 1.2,
    width = 16, height = 9, device = agg_png
)

# check concordance of effect direction
df_concordance <- df_msmr %>%
    select(phenotype,
        gene = probeID, cell_type, major_cell_type, b_SMR, se_SMR,
        p_SMR, p_SMR_multi
    ) %>%
    mutate(
        p_smr_calc = 2 * pnorm(-abs(b_SMR / se_SMR)),
        beta_sign = sign(b_SMR)
    ) %>%
    group_by(phenotype, gene) %>%
    summarise(concordant_beta = var(beta_sign) == 0)

df_concordance %>% filter(!concordant_beta)

# df_tally_celltype %>%
#   ggplot(aes(x = log10(n_genes), y = pheno_label)) +
#   theme_bw() +
#   facet_grid(rows = vars(pheno_cat), space = "free", scale = "free") +
#   geom_col(aes(fill = n_celltypes),
#            position = "stack") +
#   scale_fill_viridis_b(n.breaks = 7) +
#   theme(legend.position.inside = c(0.9,0.9),
#         legend.position = "inside",
#         panel.grid.major.y = element_line(),
#         legend.justification = c(1,1))
#
# ungroup(df_sum_pheno) %>%
#   ggplot(aes(y = pheno_label, x=log10(total_n))) +
#   theme_bw() +
#   geom_col(aes(fill = pheno_cat)) +
#   scale_fill_paletteer_d("ggthemes::Tableau_10") +
#   facet_grid(rows = vars(pheno_cat), space = "free", scale = "free") +
#   theme(axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.title.y = element_blank())
#
