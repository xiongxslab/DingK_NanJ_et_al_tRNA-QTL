res <- data.frame()
for(para1 in c('hc.cis&trans','hc.cis')){
  if(para1 == 'hc.cis&trans'){
    hc.list = c(10,20,30)
  }else{
    hc.list = c(10,30)
  }
  
  for(para in c('lcl.eur','lcl.yri','brain.eur')){
    file = paste0('01.',para,'.tRNA.classification.bed')
    rs.qtl <- fread(file)
    
    sign.count <- nrow(rs.qtl %>% filter(classification %in% hc.list))
    all.count <- nrow(rs.qtl)
    sign.count.df <- rs.qtl %>% filter(classification %in% hc.list) %>% group_by(chr) %>% summarise(N.sign = n())
    all.count.df <- rs.qtl %>% group_by(chr) %>% summarise(N.all = n())
    
    df <- full_join(sign.count.df,all.count.df,by='chr') %>% na.omit()
    
    
    
    df$pval <- NA
    df$OR <- NA
    df$conf.int.l <- NA
    df$conf.int.u <- NA
    for(i in 1:nrow(df)){
      t <- matrix(c(df$N.sign[i],
                    sign.count-df$N.sign[i], 
                    df$N.all[i]-df$N.sign[i], 
                    all.count - sign.count - (df$N.all[i]-df$N.sign[i])),nrow = 2) 
      fisher.t <- fisher.test(t)
      df$pval[i] <- fisher.t$p.value
      df$OR[i] <- fisher.t$estimate
      df$conf.int.l[i] <- fisher.t$conf.int[1]
      df$conf.int.u[i] <- fisher.t$conf.int[2]
    }
    res <- rbind(res,df %>% mutate(label = para, anno = para1))
  }
}
res$se <- (as.numeric(res$conf.int.u) - as.numeric(res$conf.int.l))/(2*1.96)
res$pval.fdr <- p.adjust(res$pval, method = "fdr")
write.table(res,"01.Fisher.bg-all_trna.txt",col.names = T,row.names = F,sep="\t",quote=F)


res <- fread("01.Fisher.bg-all_trna.txt") %>% 
  mutate(sign = case_when(pval.fdr < 0.01 ~ "***",
                          pval.fdr < 0.05 ~ "**",
                          pval < 0.05 ~ "*",
                          TRUE ~ ' '),
         chr = paste0('chr',chr))

library(ggplot2)
res$chr <- factor(res$chr,levels = c(paste0('chr',seq(1:22))))
# 普通柱状图
p <- ggplot(res,aes(x=chr,y=log2(OR),fill=sign))+
  geom_bar(stat = "identity",alpha=1) +
  # geom_jitter(aes(x= V2,y= OR,color= color),width = 0.2, height = 0,size=2)+
  scale_fill_manual(values = c("***"="#8B3626","**"="#DFB6BC","*"="#F4E7E8"," "="grey"))+
  facet_wrap(~label+anno,scales = "free",ncol = 2)+
  theme_classic()+
  # theme(legend.position = 'none') +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  xlab(NULL)


ggplot2::ggsave(filename = paste0('./01.Fisher.bg-all_trna.pdf'),
                width=8,height = 6,
                plot = p,
                device = "pdf")
