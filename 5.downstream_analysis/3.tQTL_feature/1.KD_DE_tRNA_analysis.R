library(DESeq2)

count_matrix <- readRDS(paste0('00.KD.tRNA.count.mat.rds')) %>%
  as.data.frame()
result <- data.frame()
for(sh in shs){
  count_matrix_sub <- count_matrix %>% dplyr::select(starts_with(sh),starts_with('NTC'))
  sample_info <- data.frame(
    row.names = colnames(count_matrix_sub),
    group = c(rep(sh,3), "Control", "Control", "Control")
  )
  dds <- DESeqDataSetFromMatrix(countData = count_matrix_sub,
                                colData = sample_info,
                                design = ~ group)
  
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("group", sh, "Control")) %>% as.data.frame() %>%
    na.omit() %>%
    mutate(label = case_when(padj < 0.01 & log2FoldChange > 1 ~ 'Up',
                             padj < 0.01 & log2FoldChange < -1 ~ 'Down',
                             TRUE ~ 'None'),
           sh_batch = sh,
           gene = gsub('sh','',str_split_i(sh,'-',1)),
           rep = str_split_i(sh,'-',2),
           tRNA = rownames(.)) %>%
    select(tRNA,everything(.))
  result <- rbind(result,res)
}
write.table(result,"05.DEseq2.DE.txt",col.names = T,row.names = F,sep="\t",quote=F)

