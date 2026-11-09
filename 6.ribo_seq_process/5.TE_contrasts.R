# ============================================================
# TE contrasts: treatment vs control (NTC) mean
# Input: 02.gene_TE_per_rep.tsv from 4.merge_tpm_TEcalc.R
# Output: per-rep contrasts and summary statistics
# ============================================================
suppressMessages({ library(data.table); library(ggplot2) })

OUT <- "./05TE_analysis"
te_all <- fread(file.path(OUT, "02.gene_TE_per_rep.tsv"))

# ---- Contrast: treatment rep vs NTC mean ----
make_contrast <- function(ctrl, trt, name) {
  # NTC mean across replicates
  ntc <- te_all[
    cond == ctrl,
    .(
      symbol   = dplyr::first(symbol),
      RNA_c    = mean(RNA_tpm[is.finite(RNA_tpm) & RNA_tpm > 0], na.rm = TRUE),
      Ribo_c   = mean(Ribo_tpm[is.finite(Ribo_tpm) & Ribo_tpm > 0], na.rm = TRUE)
    ),
    by = gene
  ][, TE_ctrl := Ribo_c / RNA_c
  ][, .(gene, symbol, RNA_c, Ribo_c, TE_ctrl)]

  # Treatment replicates
  trt_data <- te_all[
    cond == trt,
    .(gene, symbol, rep, TE_trt = TE, RNA_t = RNA_tpm, Ribo_t = Ribo_tpm)
  ]

  # Merge + filter + fold changes (raw FC, no log2)
  m <- inner_join(trt_data, ntc, by = c("gene", "symbol"))
  m <- m[is.finite(TE_ctrl) & is.finite(TE_trt) & TE_ctrl > 0 & TE_trt > 0]
  m[, TE_FC   := TE_trt / TE_ctrl]
  m[, RNA_FC  := RNA_t / RNA_c]
  m[, Ribo_FC := Ribo_t / Ribo_c]
  m[, contrast := name]
  m
}

# ---- Define contrasts ----
# Example: compare each treatment condition to the control
# Modify ctrl/trt names to match your experimental design
contrasts_detail <- rbind(
  make_contrast("NTC", "treatment1", "treatment1_vs_NTC"),
  make_contrast("NTC", "treatment2", "treatment2_vs_NTC")
  # Add more contrasts as needed
)
contrasts_detail <- contrasts_detail[order(contrast, -abs(TE_FC))]
fwrite(contrasts_detail, file.path(OUT, "03.TE_contrasts_per_rep.tsv"), sep = "\t")

# Summary: mean FC across reps (raw scale, then log2 for reporting)
contrasts_summary <- contrasts_detail[, .(
  TE_FC = mean(TE_FC), RNA_FC = mean(RNA_FC), Ribo_FC = mean(Ribo_FC),
  TE_FC_sd = sd(TE_FC), n_reps = .N,
  rep_values = paste(sprintf("%.3f", TE_FC), collapse = ",")
), by = .(gene, symbol, contrast)]
contrasts_summary[, log2FC := log2(TE_FC)]
contrasts_summary <- contrasts_summary[order(contrast, -abs(TE_FC))]
fwrite(contrasts_summary, file.path(OUT, "04.TE_contrasts_summary.tsv"), sep = "\t")

cat("\nContrasts (per rep):", nrow(contrasts_detail),
    "  Summary (gene level):", nrow(contrasts_summary), "\n")


