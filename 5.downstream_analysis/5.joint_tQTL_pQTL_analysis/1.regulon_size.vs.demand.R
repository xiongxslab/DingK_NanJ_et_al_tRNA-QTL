# transform anticodon to codon
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
      return(base)
    }
  })
  complement_string <- paste(complement_vector, collapse = "")
  return(complement_string)
}

pQTL <- fread("00.file.proteinID.trait.pQTL.SOMAscan.txt") %>% filter(pQTL == 'pQTL')
tp.pair <- fread("01.lcl.eur.tp.MR.pairs.txt") %>%
  na.omit() %>%
  mutate(anticodon = paste0(str_split_i(exposure,'-',2),'-',str_split_i(exposure,'-',3))) %>%
  select(anticodon,gene_id) %>%
  distinct() %>%
  group_by(anticodon) %>%
  summarise(N.protein = n()) %>%
  mutate(codon = apply(as.matrix(anticodon),1,function(x){reverse_string(str_split_i(x,'-',2))}))
exp <- readRDS("01.lcl.eur.mRNA.tpm.filter_0.8.rds") %>% rowMeans() %>% as.data.frame()
colnames(exp) <- 'mrna.tpm'
exp$gene_id <- rownames(exp)

cds <- fread("02.GRCh38.cds_longest_transcript.codon.count.txt") %>% 
  inner_join(exp,by=c('gene_id')) %>%
  dplyr::select(-transcript_id,-width)

cds.m <- cds %>% mutate(mrna.tpm = as.numeric(mrna.tpm)) %>% melt(id.vars = c('gene_id','mrna.tpm')) %>%
  mutate(value = as.numeric(value),
         command = mrna.tpm*value) %>%
  group_by(variable) %>%
  reframe(command.median = median(command), 
          command.mean = mean(command),
          command.sum = sum(command)) %>%
  inner_join(tp.pair, by = c('variable'='codon'))
p1 <- ggplot(cds.m,aes(y = log10(N.protein), x = log10(command.sum))) +
  geom_point(color = '#407BD0') +
  stat_cor(method = 'spearman',label.y=-1) +
  geom_text_repel(aes(label = variable),
                  segment.color = "grey50", 
                  segment.size = 0.5) +
  geom_smooth(method='lm', color = "#8B7355", fill = "#EED8AE") +
  theme_classic2()

ggplot2::ggsave(filename = paste0('./02.tp.pleiotropic.demand.comparison.pdf'),
                width=3,height = 3,
                plot = p1,
                device = "pdf")
