library(ggplot2)
library(ggpubr)
library(dplyr)
library(stringr)

# 1.correlation with tRNA sequencing
res <- data.frame()
for(para1 in c('ARM-seq','DM-tRNA-seq','Brain')){
  if(para1 == 'Brain'){
    val <- readRDS("01.brain.eur.tRNA.tpm.filter_0.8.rds") %>%
      rowMeans() %>% as.data.frame()  %>%
      rownames_to_column('ID') %>%
      mutate(ID = gsub('Homo_sapiens_','',ID),
             family = paste(str_split_i(ID,'-',1),str_split_i(ID,'-',2),str_split_i(ID,'-',3),str_split_i(ID,'-',4), sep = '-')) %>%
      group_by(family) %>%
      summarise(tRNA = sum(`.`)) %>%
      column_to_rownames('family') %>%
      as.data.frame()
  }else{
    val <- readRDS(paste0("00.",para1,".tpm(family).rds")) %>%
      rowMeans() %>% as.data.frame() 
  }
  
  if(para1 == 'Brain'){
    para2s <- 'LCL'
  }else{
    para2s <- c('LCL','Brain','DM-tRNA-seq')
  }
  for(para2 in para2s){
    if(para2 == 'LCL'){
      our <- readRDS("01.lcl.eur.tRNA.tpm.filter_0.8.rds")  %>%
        as.data.frame() %>%
        rowMeans() %>%
        as.data.frame() %>%
        rownames_to_column('ID') %>%
        mutate(ID = gsub('Homo_sapiens_','',ID),
               family = paste(str_split_i(ID,'-',1),str_split_i(ID,'-',2),str_split_i(ID,'-',3),str_split_i(ID,'-',4), sep = '-')) %>%
        group_by(family) %>%
        summarise(tRNA = sum(`.`)) %>%
        column_to_rownames('family')
    }else if(para2 == 'Brain'){
      our <- readRDS("01.brain.eur.tRNA.tpm.filter_0.8.rds")  %>%
        as.data.frame() %>%
        rowMeans() %>%
        as.data.frame() %>%
        rownames_to_column('ID') %>%
        mutate(ID = gsub('Homo_sapiens_','',ID),
               family = paste(str_split_i(ID,'-',1),str_split_i(ID,'-',2),str_split_i(ID,'-',3),str_split_i(ID,'-',4), sep = '-')) %>%
        group_by(family) %>%
        summarise(tRNA = sum(`.`)) %>%
        column_to_rownames('family')
    }else{
      our <- readRDS(paste0("00.",para2,".tpm(family).rds")) %>%
        as.data.frame() %>%
        rowMeans() %>%
        as.data.frame()
      colnames(our) <- 'tRNA'
    }
    if(para1 == para2) next;
    
    if(para1 %in% c('ARM-seq','DM-tRNA-seq')){
      trna <- intersect(rownames(our),rownames(val))
      df <- cbind(our[match(trna,rownames(our)),],
                  val[match(trna,rownames(val)),]) %>%
        as.data.frame()
      colnames(df) <- c('para2.tpm','para1.tpm')
      rownames(df) <- trna
      # df <- df %>% filter('ARM_seq' > 0)
    }else{
      trna <- intersect(rownames(our),rownames(val))
      df <- cbind(our[match(trna,rownames(our)),],
                  val[match(trna,rownames(val)),]) %>%
        as.data.frame()
      colnames(df) <- c('para2.tpm','para1.tpm')
      rownames(df) <- trna
    }
    res <- rbind(res,df %>% mutate(x = para2, y = para1, i = paste0(para2,'_',para1)))
  }
}

p <- ggplot(res,aes(x=log10(para2.tpm) ,y=log10(para1.tpm)))+geom_point(color='#CD853F')+
  facet_wrap(vars(i), ncol = 3,scale = 'free') +
  stat_cor(number.digits = 3,label.y=0.2)+
  stat_cor(number.digits = 3,label.y=0.6,method = 'spearman')+
  geom_smooth(method='lm', color = "#8B7355", fill = "#EED8AE")+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),panel.grid.minor = element_blank()
  ) +
  xlab('log10(para1)') + ylab('log10(para2)')

ggplot2::ggsave(filename = "01.correlation.of.with.validation.in.treatSample.family.pdf",
                width=3*3.5,height=2*3.5,
                plot = gridExtra::marrangeGrob(list(p),nrow =1,ncol=1),
                device = "pdf")
