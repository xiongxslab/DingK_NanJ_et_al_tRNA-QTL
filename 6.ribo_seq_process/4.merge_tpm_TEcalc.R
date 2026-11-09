# ============================================================
# Merge TPM from kallisto abundance.tsv (RNA-seq + Ribo-seq)
# Output: merged TPM matrices and per-replicate TE table
# ============================================================
suppressMessages({ library(data.table); library(dplyr); library(tidyr) })

# MANE tx2gene annotation
tx2g <- fread("./reference/tx2gene_v.tsv")
setnames(tx2g, c("tx", "gene", "symbol"))
cat("tx2gene:", nrow(tx2g), "\n")

# ---- Load TPM ----
load_tpm <- function(dir, files) {
  lst <- lapply(files, function(f) {
    d <- fread(file.path(dir, f, "abundance.tsv"), select = c("target_id", "tpm"))
    setnames(d, c("tx", f))
    d
  })
  m <- Reduce(function(a, b) merge(a, b, by = "tx", all = TRUE), lst)
  m[is.na(m)] <- 0
  m
}

# Input directories
path.rna <- "./04rna_kallisto_MANE"
path.ribo <- "./04ribo_kallisto_MANE"
rna_samples <- sort(basename(list.dirs(path.rna, recursive = FALSE)))
ribo_samples <- sort(basename(list.dirs(path.ribo, recursive = FALSE)))
cat("RNA samples:", length(rna_samples), "  Ribo samples:", length(ribo_samples), "\n")

tpm_rna <- load_tpm(path.rna, rna_samples)
tpm_ribo <- load_tpm(path.ribo, ribo_samples)

# Merge RNA + Ribo TPM by transcript
tpm_all <- merge(tpm_rna, tpm_ribo, by = "tx", all = TRUE)
tpm_all <- merge(tpm_all, tx2g, by = "tx")

# ---- TE per replicate ----
# Match Ribo samples to RNA samples by name (e.g. "cond-rep1" in both)
base_of <- function(x) sub("-rep[0-9]+$", "", x)

te_list <- list()
for (ribo_rep in ribo_samples) {
  if (!ribo_rep %in% colnames(tpm_all)) {
    cat("  [SKIP]", ribo_rep, "- no RNA match\n")
    next
  }
  rna_tpm <- tpm_all[[ribo_rep]]
  ribo_tpm <- tpm_all[[ribo_rep]]
  # Use .x/.y columns if RNA and Ribo sample names are identical
  if (paste0(ribo_rep, ".x") %in% colnames(tpm_all)) {
    rna_tpm <- tpm_all[[paste0(ribo_rep, ".x")]]
    ribo_tpm <- tpm_all[[paste0(ribo_rep, ".y")]]
  }
  te <- ifelse(rna_tpm > 0, ribo_tpm / rna_tpm, NA)
  te_list[[ribo_rep]] <- data.table(
    tx = tpm_all$tx, gene = tpm_all$gene, symbol = tpm_all$symbol,
    cond = base_of(ribo_rep), rep = sub(".*-rep", "", ribo_rep), sample = ribo_rep,
    RNA_tpm = rna_tpm, Ribo_tpm = ribo_tpm, TE = te
  )
}
te_all <- rbindlist(te_list)
cat("TE rows (per rep, tx-level):", nrow(te_all), "\n")

# ---- Aggregate to gene level ----
te_all <- te_all[, .(
  symbol   = dplyr::first(symbol),
  RNA_tpm  = sum(RNA_tpm),
  Ribo_tpm = sum(Ribo_tpm)
), by = .(gene, cond, rep, sample)][, TE := ifelse(RNA_tpm > 0, Ribo_tpm / RNA_tpm, NA)]
cat("TE rows (per rep, gene-level):", nrow(te_all), "\n")

# Save
OUT <- "./05TE_analysis"
dir.create(OUT, showWarnings = FALSE)
saveRDS(tpm_rna, file.path(OUT, "01.rna_tpm.rds"))
saveRDS(tpm_ribo, file.path(OUT, "01.ribo_tpm.rds"))
fwrite(te_all, file.path(OUT, "02.gene_TE_per_rep.tsv"), sep = "\t")
cat("Saved:", OUT, "\n")
