# ============================================================
# Ribo-seq QC: RPF length, gene region distribution, P-site QC
# Input:  04.QC/1.rpf_length_distribution.txt
#         04.QC/2.rpf_region_counts.txt
#         06.RiboseQC/26-32nt/P_site_distribution.tsv
# Output: 04.QC/3.RiboSeq_QC.pdf, 3.RiboSeq_QC.png
# ============================================================
suppressMessages({
  library(data.table); library(ggplot2); library(ggpubr)
})

setwd("./04.QC")

# Input files
len   <- fread("1.rpf_length_distribution.txt")
reg   <- fread("2.rpf_region_counts.txt")
psite <- fread("./06.RiboseQC/26-32nt/P_site_distribution.tsv")

# Sample grouping: derive from sample names (e.g. "cond-rep1" → cond="cond", rep="rep1")
grp.lvl <- unique(sub('-rep[0-9]$', '', len$sample))
mkgrp <- function(x) factor(sub('-rep[0-9]$', '', x), levels = grp.lvl)
mkrep <- function(x) sub('.*-(rep[0-9])$', '\\1', x)
len[,   group := mkgrp(sample)]
reg[,   group := mkgrp(sample)]
psite[, group := mkgrp(sample)]
len[,   rep := mkrep(sample)]
psite[, rep := mkrep(sample)]
ngrp <- length(grp.lvl)

# Colors
reg.col  <- c(CDS = '#2ECC71', UTR5 = '#3498DB', UTR3 = '#E74C3C', other = '#95A5A6')
reg.col2 <- c(CDS = '#F4A460', UTR5 = '#8B7E66', UTR3 = '#CDBA96')   # mRNA-only colors

############################################################
## p1: RPF length distribution (0-50nt, grey band = 26-32nt window)
############################################################
p1 <- ggplot(len, aes(x = length, y = pct * 100, color = rep)) +
  annotate('rect', xmin = 26, xmax = 32, ymin = -Inf, ymax = Inf,
           fill = 'grey80', alpha = 0.4) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.1) +
  theme_classic() +
  facet_wrap(vars(group), ncol = ngrp) +
  scale_color_manual(values = c(rep1 = '#1f78b4', rep2 = '#e31a1c')) +
  scale_x_continuous(breaks = seq(0, 50, 10), limits = c(0, 50)) +
  labs(x = 'RPF length (nt)', y = '% trimmed reads', color = NULL,
       title = 'RPF length distribution (pre 26-32nt selection)') +
  theme(plot.title = element_text(size = 10.5, hjust = 0.5),
        legend.position = 'bottom')

############################################################
## p2: RPF 5' end region distribution (pie charts)
############################################################
regl <- melt(reg, id.vars = c('sample','group'),
             measure.vars = c('UTR5', 'CDS', 'UTR3'), variable.name = 'region')
regl[, pct := value / sum(value), by = sample]
p2.data <- regl[, .(pct = mean(pct)), by = .(group, region)]
p2.data[, lbl := sprintf('%.1f%%', pct * 100)]
p2.data$region <- factor(p2.data$region, levels = c('UTR5', 'CDS', 'UTR3'))

p2 <- ggplot(p2.data, aes(x = '', y = pct, fill = region)) +
  geom_col(width = 1, color = 'white') +
  coord_polar(theta = 'y', start = 0) +
  scale_y_continuous(expand = expansion(0, 0)) +
  geom_text(aes(label = lbl), position = position_stack(vjust = 0.5), size = 2.6) +
  facet_wrap(vars(group), ncol = ngrp) +
  scale_fill_manual(values = reg.col) +
  theme_void() +
  labs(fill = NULL, title = "RPF 5' end distribution (5'UTR/CDS/3'UTR)") +
  theme(plot.title = element_text(size = 11, hjust = 0.5),
        legend.position = 'bottom')

############################################################
## p3: P-site genomic distribution (all regions, RiboseQC)
############################################################
psite.l <- melt(psite, id.vars = c('sample','group','rep','total_p_sites'),
                measure.vars = c('UTR5', 'CDS', 'UTR3', 'other'),
                variable.name = 'region', value.name = 'pct')
psite.l$region <- factor(psite.l$region, levels = c('UTR5', 'CDS', 'UTR3', 'other'))

p3 <- ggplot(psite.l, aes(x = rep, y = pct, fill = region)) +
  geom_col(width = 0.7, color = 'white') +
  theme_classic() +
  facet_wrap(vars(group), ncol = ngrp) +
  scale_fill_manual(values = reg.col) +
  labs(x = NULL, y = '% P-sites', fill = NULL,
       title = 'P-site genomic distribution (RiboseQC)') +
  theme(plot.title = element_text(size = 11, hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'bottom')

############################################################
## p4: P-site distribution within mRNA (exclude 'other', renormalize)
############################################################
psite.m <- melt(psite[, .(sample, group, rep, CDS, UTR5, UTR3)],
                id.vars = c('sample','group','rep'),
                measure.vars = c('CDS', 'UTR5', 'UTR3'),
                variable.name = 'region', value.name = 'pct')
psite.m[, pct := pct / sum(pct) * 100, by = sample]
psite.m$region <- factor(psite.m$region, levels = c('CDS', 'UTR5', 'UTR3'))

p4 <- ggplot(psite.m, aes(x = rep, y = pct, fill = region)) +
  geom_col(width = 0.7, color = 'white') +
  geom_text(aes(label = sprintf('%.1f', pct)),
            position = position_stack(vjust = 0.5), size = 2.6) +
  theme_classic() +
  facet_wrap(vars(group), ncol = ngrp) +
  scale_fill_manual(values = reg.col2) +
  labs(x = NULL, y = '% P-sites in mRNA', fill = NULL,
       title = "P-site distribution within mRNA ('other' excluded)") +
  theme(plot.title = element_text(size = 11, hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'bottom')

# ---- Output ----
ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2,
          labels = c('p1', 'p2', 'p3', 'p4'), font.label = list(size = 12))
ggsave("3.RiboSeq_QC.pdf", width = 10, height = 8)
ggsave("3.RiboSeq_QC.png", width = 10, height = 8, dpi = 150)

cat("mRNA CDS proportion:\n")
print(psite[, .(sample, cds_in_mrna = round(CDS / (CDS + UTR5 + UTR3) * 100, 1))])
cat("\nSaved: 3.RiboSeq_QC.pdf, 3.RiboSeq_QC.png\n")
