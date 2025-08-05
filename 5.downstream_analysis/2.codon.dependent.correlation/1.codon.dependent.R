library(data.table)
library(dplyr)
library(ggplot2)
library(ggsignif)
library(gridExtra)
library(ggpubr)

# reverse anticodon to codon
reverse_string <- function(tRNA_name) {
  parts <- strsplit(tRNA_name, "-")[[1]]
  amino_acid <- parts[2]
  codon <- parts[3]
  
  char_vector <- strsplit(codon, "")[[1]]
  reversed_vector <- rev(char_vector)
  complement_vector <- sapply(reversed_vector, function(base) {
    switch(base,
           "A" = "T",
           "T" = "A",
           "C" = "G",
           "G" = "C",
           base)
  })
  complement_codon <- paste(complement_vector, collapse = "")
  
  if(amino_acid == "iMet"){
    return(paste0("iMet_", complement_codon))
  } else {
    return(complement_codon)
  }
}

############################## mRNA*count & tRNA ############################## 
res <- data.frame()
for(para in c('lcl.eur','lcl.yri','brain.eur')){
  trna.info <- fread(paste0("./01.",para,".tRNA.classification.bed")) %>%
    filter(classification %in% c(10,30))
  trna.exp <- readRDS(paste0("./02.",para,".tRNA.tpm.filter_0.8.inter.rds")) %>%
    as.data.frame() %>%
    filter(rownames(.) %in% trna.info$ID)
  trna.exp$mean <- rowMeans(trna.exp)
  trna.exp$codon <- apply(as.matrix(rownames(trna.exp)),1, function(x){reverse_string(x)})
  trna.exp.mean <- trna.exp %>% select(mean,codon) %>%
    group_by(codon) %>%
    reframe(sum_isodecoder.tpm = sum(mean))
  
  gene.exp <- readRDS(paste0("./02.",para,".mRNA.tpm.filter_0.8.inter.rds")) %>% 
    as.data.frame() 
  gene.exp$mean <- rowMeans(gene.exp)
  gene.exp$gene_id <- rownames(gene.exp)
  
  gene.cds.count <- fread("./02.GRCh38.cds_longest_transcript.codon.count.txt") %>% as.data.frame()
  gene.cds.count.exp <- gene.exp %>% select(gene_id, mean) %>% inner_join(y=gene.cds.count,by=c('gene_id'))

  cutoff = 1
  for(codons in unique(trna.exp.mean$codon)){
    trna_tmp <- subset(trna.exp.mean,codon == codons)
    gene_tmp <- gene.cds.count.exp %>% 
      select(gene_id,mean,codon=all_of(codons)) %>% 
      mutate(demand = mean * codon)
    trna_tmp$codon_demand <- sum(gene_tmp$demand)
    trna_tmp$codon_count = sum(gene_tmp$codon)
    trna_tmp$cutoff = cutoff
    trna_tmp$batch = para
    res <- rbind(res,trna_tmp)
  }
}
saveRDS(res,paste0("01.supply.vs.demand.hc.rds"))


demand.cor <- readRDS('01.supply.vs.demand.hc.rds') %>%
  filter(codon != 'iMet_ATG' & cutoff == 1)
ggplot(demand.cor,aes(x = log10(sum_isodecoder.tpm), y = log10(codon_demand), color = batch)) +
  facet_wrap(vars(batch),scales = 'free') +
  geom_point() +
  stat_cor(method="spearman",number.digits = 3)+
  geom_smooth(method='lm', color = "#8B7355", fill = "#EED8AE")+
  scale_color_manual(values = c('lcl.eur' = '#315c82','lcl.yri' = '#dca44c','brain.eur' = '#af8eb0')) + 
  theme_classic2()+
  theme(
    panel.grid.major = element_blank(),panel.grid.minor = element_blank()
  )

