library(data.table)
library(clusterProfiler)
library("org.Hs.eg.db")
library(ggplot2)
library(dplyr)
###de-duplication for GO results
gene_compare<- function(a,b, cutoff){
  gen_split<- function(x){str_split(x, '\\/') %>% unlist}
  gen.dup<- intersect(gen_split(a), gen_split(b))
  if(length(gen.dup)/min(length(gen_split(a)), length(gen_split(b))) < cutoff){
    g.com<- TRUE
  }else{
    g.com<- FALSE
  }
}

go_select<- function(df){
  df.filter<- df %>%
    mutate(rank=frankv(pvalue, ties.method="min"))
  if(nrow(df.filter)==1){
    df.use<- df.filter
  }else{
    df.filter<- df.filter[order(df.filter$rank),] ## sort by rank
    df.use <- df.filter[1,]
    for(index in 2:nrow(df.filter)){
      df.use1<- df.filter[index,]
      num_rows <- nrow(df.use)
      go.compare<- list()
      for(i in 1:num_rows){
        go.compare1<- gene_compare(df.use1[[1, 9]],df.use[[i, 9]], 0.5)
        go.compare<- c(go.compare, go.compare1)
      }
      if(all(go.compare)){
        df.use<- rbind(df.use, df.use1)
      }else{
        df.use<- df.use
      }
    }
  }
  return(df.use)
}

##### GO #####
pair <- fread("02.EUR.tp.mr.info.res") %>% 
  dplyr::select(exposure,gene_id) %>% distinct()
pair.table <- pair %>% group_by(exposure) %>% summarise(N = n()) %>% filter(N > 5)
pair <- pair %>% filter(exposure %in% pair.table$exposure)
list <- fread("00.file.proteinID.trait.pQTL.SOMAscan.txt") %>%
  na.omit() %>%
  filter(pQTL == 'pQTL')
for(para1 in unique(pair$exposure)){
  trna <- subset(pair,exposure == para1)
  enrich.go <- enrichGO(gene = trna$gene_id,  
                        OrgDb = 'org.Hs.eg.db', 
                        # keyType = 'ENSEMBL',
                        keyType = 'ENSEMBL',
                        ont = 'ALL',  
                        pAdjustMethod = 'fdr',  
                        universe = list$gene_id,
                        pvalueCutoff = 1, 
                        qvalueCutoff = 1,
                        readable = FALSE)
  enrich.go.cutoff <- enrich.go %>% as.data.frame() %>% filter(pvalue < 0.1) %>% 
    mutate(tRNA = para1)
  write.table(enrich.go.cutoff,paste0("1.tp.pQTL.GO.bg-gprotein.sign.res"),sep="\t",quote=F,col.names = T,row.names = F,append = T)
}

trna.info <- fread("hg38-mature-tRNAs-DNA.bed") %>% 
  dplyr::select(V1,V4) %>% 
  dplyr::rename(chr = V1, tRNA=V4) %>%
  mutate(ID.chr = paste0(gsub('tRNA-','',tRNA),'|',chr))
res <- fread("1.tp.pQTL.GO.bg-gprotein.sign.res") %>%
  filter(Description != 'Description' & pvalue < 0.01) %>%
  inner_join(trna.info,by='tRNA') %>%
  go_select()
res.table <- res %>% group_by(Description) %>% 
  summarise(N = n())  %>% filter(N > 1)
res.filter <- res %>% distinct() %>% filter(Description %in% res.table$Description) %>%
  mutate(Count = as.numeric(Count),
         pvalue = as.numeric(pvalue))


pdf("1.tp.pQTL.GO.plot.pdf",width = 7,height = 7)
ggplot(res.filter,aes(x = reorder(ID.chr,Count,decreasing = T), y = reorder(Description,Count),col = ONTOLOGY,size = Count))+
  geom_point() +
  scale_color_manual(values = c('#da5f50','#f7c96e','#37638c'))+
  # scale_color_gradient(low="#DFB6BC", high="#990033")+
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1,size=10),
        axis.text.y = element_text(size = 10)) +
  xlab(NULL) + ylab(NULL)
dev.off()

