# ============================================================
# Codon frequency analysis: TE fold change by codon usage
# Input:  04.TE_contrasts_summary.tsv from 5.TE_contrasts.R
#         ./reference/GRCh38.cds_longest_transcript.codon.frq.txt
# Output: 05TE_analysis/06.codon_freq_analysis.pdf
# ============================================================
suppressMessages({ library(data.table); library(dplyr); library(ggplot2); library(tidyr) })

OUT <- "./05TE_analysis"
contrasts_summary <- fread(file.path(OUT, "04.TE_contrasts_summary.tsv"))
codon.frq <- fread("./reference/GRCh38.cds_longest_transcript.codon.frq.txt")

# ---- Select codons to analyze ----
# Modify to analyze codons of interest (e.g. corresponding to perturbed tRNAs)
codons_of_interest <- c("AGG", "CAG")  # Example: Arg-AGG and Gln-CAG

# ---- Prepare codon frequency groups (3 tiles) ----
codon_groups <- lapply(codons_of_interest, function(cd) {
  codon.frq %>%
    filter(gene_id %in% contrasts_summary$gene) %>%
    dplyr::select(gene_id, all_of(cd)) %>%
    dplyr::mutate(frq_group = as.character(dplyr::ntile(.data[[cd]], 3))) %>%
    dplyr::rename(frequency = all_of(cd)) %>%
    mutate(codon = cd)
})
codon_groups <- bind_rows(codon_groups)

# ---- Join TE data with codon frequency ----
te_codon <- contrasts_summary %>%
  inner_join(codon_groups, by = c("gene" = "gene_id"))

# ---- Wilcoxon test: TE_FC across frequency groups (3 tiles) ----
test_codon <- function(df) {
  df %>%
    group_by(contrast, codon) %>%
    summarise(
      p_1_2 = wilcox.test(TE_FC[frq_group == "1"], TE_FC[frq_group == "2"])$p.value,
      p_2_3 = wilcox.test(TE_FC[frq_group == "2"], TE_FC[frq_group == "3"])$p.value,
      p_1_3 = wilcox.test(TE_FC[frq_group == "1"], TE_FC[frq_group == "3"])$p.value,
      .groups = "drop"
    ) %>%
    pivot_longer(cols = starts_with("p_"), names_to = "comparison", values_to = "p") %>%
    mutate(
      label = case_when(
        p < 0.001 ~ "***",
        p < 0.01  ~ "**",
        p < 0.05  ~ "*",
        TRUE ~ "ns"
      ),
      x = case_when(
        comparison == "p_1_2" ~ 1.5,
        comparison == "p_2_3" ~ 2.5,
        comparison == "p_1_3" ~ 2
      )
    )
}

p_df <- test_codon(te_codon)

# ---- Plot: TE_FC by codon frequency tercile ----
p1 <- ggplot(te_codon, aes(x = frq_group, y = log2(TE_FC), color = frq_group)) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.15) +
  scale_color_manual(values = c('#fdbe85','#fd8d3c','#d94701')) +
  facet_wrap(contrast ~ codon, ncol = 2, scales = 'free_y') +
  theme_classic2() +
  geom_text(
    data = p_df,
    aes(x = x, y = max(log2(te_codon$TE_FC), na.rm = TRUE) * 0.9, label = label),
    inherit.aes = FALSE, color = "black", size = 4
  ) +
  labs(x = 'Codon frequency tercile', y = 'log2(TE fold change)', title = 'TE by codon usage') +
  theme(legend.position = 'none')

# ---- Output ----
pdf(file.path(OUT, "06.codon_freq_analysis.pdf"), width = 12, height = 6)
print(p1)
dev.off()

cat("Saved:", file.path(OUT, "06.codon_freq_analysis.pdf"), "\n")
