library(ggplot2)

start.list = c('pf1','pf2')
pc=5
sum.df = data.frame()
for(start in start.list){
  if(start == 'pf1'){
    end.list <- c(1:20)
  }else if(start == 'pf2'){
    end.list <- c(2:20)
  }else{
    end.list <- c(1:2,4:20)
  }
  for(end in end.list){
    df.merge = data.frame()
    for(chr in 1:22){
      out = read.table(paste0('./output.1M/chr',chr,'.',start,'-',end,'.pc1-',pc,'.out.txt.gz'),head=F)
      if(ncol(out)==1){next}
      df.merge = rbind(df.merge,out)
      print(chr)
    }
    sum.df = rbind(sum.df,data.frame(pf_start=start,pf=end,p_cut=0.0001,nqtl=nrow(df.merge),n_gTRNA=nrow(df.merge[!duplicated(df.merge$V1),]),min_P=min(df.merge$V7)))
    df.merge = df.merge[df.merge$V7<0.00001,]
    sum.df = rbind(sum.df,data.frame(pf_start=start,pf=end,p_cut=0.00001,nqtl=nrow(df.merge),n_gTRNA=nrow(df.merge[!duplicated(df.merge$V1),]),min_P=min(df.merge$V7)))
  }
}
saveRDS(sum.df,'01.tQTL.EdgeR_norm.fastqtl_test.rds')


p1 <- ggplot(sum.df,aes(x=pf,y=nqtl,group=pf_start,color=as.character(pf_start))) + geom_point() + geom_line() + theme_bw() +
  facet_grid(vars(p_cut),scales="free_y") +
  ggtitle('# QTL') +
  labs(color = "pf start")

p2 <- ggplot(sum.df,aes(x=pf,y=n_gTRNA,group=pf_start,color=as.character(pf_start))) + geom_point() + geom_line() + theme_bw() +
  facet_grid(vars(p_cut),scales="free_y") + 
  ggtitle('# g-tRNA')+
  labs(color = "pf start")

p3 <- ggplot(sum.df,aes(x=pf,y=-log10(min_P),group=pf_start,color=as.character(pf_start))) + geom_point() + geom_line() + theme_bw() +
  facet_grid(vars(p_cut),scales="free_y") + ggtitle('Min p-value')+
  labs(color = "pf start")
ggplot2::ggsave(filename = '01.tQTL.EdgeR_norm.fastqtl_test.pdf',
                width=12,height = 4,
                plot = gridExtra::marrangeGrob(list(p1,p2,p3),nrow =1,ncol=3),
                device = "pdf")





