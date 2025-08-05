#### codon frequency ####
reverse_string <- function(input_string) {
  char_vector <- strsplit(input_string, "")[[1]]
  reversed_vector <- rev(char_vector)
  complement_vector <- sapply(reversed_vector, function(base) {
    if (base == "A") {
      return("T")
    } else if (base == "T") {
      return("A")
    } else if (base == "C") {
      return("G")
    } else if (base == "G") {
      return("C")
    } else {
      return(base)  # 非碱基字符保持不变
    }
  })
  complement_string <- paste(complement_vector, collapse = "")
  return(complement_string)
}
pair <- fread("02.EUR.tp.mr.info.res") %>%
  na.omit() %>%
  mutate(codon = apply(as.matrix(exposure), 1, function(x){reverse_string(str_split_i(x,'[-]',3))}),
         AA = str_split_i(exposure,'[-]',2)) %>%
  dplyr::select(gene_id,codon) %>%
  distinct() %>%
  mutate(Is_sign = 'Yes')
frq <- fread("02.all_gene.mrna.and.codon_frq.txt") %>%
  left_join(pair,by=c('gene_id','variable'='codon')) %>%
  replace_na(list(Is_sign = "No")) %>%
  na.omit() %>%
  mutate(mrna.tpm_quantile = ntile(mrna.tpm, 3)) %>%
  group_by(gene_id,Is_sign,mrna.tpm_quantile) %>%
  mutate(median.frq = median(frq))

frq$Is_sign <- factor(frq$Is_sign, levels = c('Yes','No'))

p1 <- ggplot(frq, aes(y = median.frq, x = Is_sign,fill=Is_sign)) +
  geom_violin(trim=FALSE,color="white") +
  facet_grid(~mrna.tpm_quantile) +
  geom_boxplot(varwidth = FALSE, outlier.shape = NA, width=0.2,fill="white") +
  scale_fill_manual(values = c('Yes'="#EFD496",'No'="#9BB89C")) +
  theme_bw() +
  xlab(NULL) +
  ylab('Codon frequency') +
  theme(legend.position = "none",
    panel.grid.major=element_line(colour=NA),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    panel.grid.minor = element_blank()) +
  stat_summary(
    fun = median,
    geom = "text",  # 将均值显示为文本
    aes(label = round(..y.., 5)),  
    vjust = -0.5
  ) +
  geom_signif(comparisons = list(c("Yes", "No")),
  map_signif_level = FALSE)

ggplot2::ggsave(filename = paste0('02.compare.tp.MR.preference(mrna.3tile).pdf'),
  width=7,height = 3,
  plot = p1,
  device = "pdf")
