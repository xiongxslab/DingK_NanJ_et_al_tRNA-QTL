# consistency #
library(ggpubr)
eur.cluster <- fread("01.lcl.eur.tRNA.classification.bed") %>%
  select(ID,classification) %>%
  filter(classification %in% c(10,30))
yri.cluster <- fread("01.lcl.yri.tRNA.classification.bed") %>%
  select(ID,classification) %>%
  filter(classification %in% c(10,30))

eur.lead <- fread("EUR.permute.out.qval_beta.0.1.txt") %>%
  filter(phenotype_id %in% eur.cluster$ID) # 55
yri.lead <- fread("YRI.permute.out.qval_beta.0.1.txt") %>%
  filter(phenotype_id %in% yri.cluster$ID) # 4

gtrna.union <- union(eur.lead$phenotype_id,yri.lead$phenotype_id)
gtrna.inter <- intersect(eur.lead$phenotype_id,yri.lead$phenotype_id)

eur <- fread("EUR.tQTL.sumstat.txt.gz",
               col.names = c('phenotype_id','variant_id','distance','ma_samples','ma_count','maf','pval_nominal','slope','slope_se')) %>%
  dplyr::select(-distance,-ma_samples,-ma_count)
yri <- fread("YRI.tQTL.sumstat.txt.gz",
               col.names = c('phenotype_id','variant_id','distance','ma_samples','ma_count','maf','pval_nominal','slope','slope_se')) %>%
  dplyr::select(-distance,-ma_samples,-ma_count)
df <- inner_join(eur,yri,by=c('phenotype_id','variant_id'),suffix=c('.eur','.yri')) %>%
  filter(phenotype_id %in% gtrna.union) %>%
  mutate(label = ifelse(phenotype_id %in% gtrna.inter, 'Shared',
                        ifelse(phenotype_id %in% eur.lead$phenotype_id, 'EUR specific', 'YRI specific')))

result <- data.frame()
for(para1 in c('EUR specific','YRI specific','Shared')){
  if(para1 == 'YRI specific'){
    ref <- yri.lead
    df.s <- df %>% filter(label == para1) %>%
      inner_join(ref %>% select(phenotype_id,npval.cut),by=c('phenotype_id')) %>%
      filter(pval_nominal.yri <= npval.cut) %>%
      select(-npval.cut)
  }else if(para1 == 'EUR specific'){
    ref <- eur.lead
    df.s <- df %>% filter(label == para1) %>%
      inner_join(ref %>% select(phenotype_id,npval.cut),by=c('phenotype_id')) %>%
      filter(pval_nominal.eur <= npval.cut) %>%
      select(-npval.cut)
  }else{
    eur.tqtl <- fread("EUR.tQTL.qval_beta.0.1.txt")
    yri.tqtl <- fread("YRI.tQTL.qval_beta.0.1.txt")
    tqtl <- inner_join(yri.tqtl,eur.tqtl,by=c('phenotype_id','variant_id'),suffix=c('.yri','.eur')) %>%
      mutate(name = paste(phenotype_id,variant_id,sep="|"))
    df.s <- df %>% filter(label == para1) %>%
      mutate(name = paste(phenotype_id,variant_id,sep="|")) %>%
      filter(name %in% tqtl$name) %>%
      select(-name)
  }
  result <- rbind(result,df.s)
}
result <- result %>% mutate(pval = ifelse(label == 'EUR specific', pval_nominal.yri, pval_nominal.eur))
write.table(result,"03.tQTL.compare.between.yri.and.eur.txt",col.names = T,row.names = F,sep = "\t",quote=F)

result <- fread("03.tQTL.compare.between.yri.and.eur.txt")
result.eur <- result %>%
  filter(label %in% c('EUR specific','Shared'))
result.eur$label <- factor(result.eur$label,levels = c('EUR specific','Shared'))
x_lim <- round(max(abs(-log10(result.eur$pval_nominal.yri)*sign(result.eur$slope.eur))),digits = 0)
y_lim <- round(max(abs(result.eur$slope.yri)),digits = 1)
p1 <- ggplot(result.eur,aes(x=-log10(pval_nominal.yri)*sign(slope.eur), y=slope.yri, col=-log10(pval_nominal.yri),
                              shape=label)) +
  # facet_wrap(vars(label),ncol = 3)+
  geom_point(size=2)+
  scale_colour_gradient2(low="#8D99AE",mid="#EDF2F4",high="#D80032",midpoint = 2)+
  geom_hline(aes(yintercept=0), colour="black", linetype="dashed") +
  geom_vline(aes(xintercept=0), colour="black", linetype="dashed") +
  xlim(c(-x_lim,x_lim)) + ylim(c(-y_lim,y_lim)) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),panel.grid.minor = element_blank()
  )

ggplot2::ggsave(filename = paste0('./03.consistency between ancestries.pdf'),
                width=6,height = 3,
                plot = p1,
                device = "pdf")



