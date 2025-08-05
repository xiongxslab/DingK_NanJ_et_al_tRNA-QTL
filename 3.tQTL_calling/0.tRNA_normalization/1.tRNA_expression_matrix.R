countToTPM <- function(count_matrix, gene_lengths) {
  total_counts <- colSums(count_matrix)
  norm_factors <- total_counts / 1e6
  geo_mean <- exp(mean(log(norm_factors)))
  tpm_matrix <- count_matrix / gene_lengths * geo_mean
  tpm_matrix <- tpm_matrix / rowSums(tpm_matrix) * 1e6
  return(tpm_matrix)
}

count <- readRDS("./06exp_matrix/00.tRNA.count.mat.rds") %>%
  as.data.frame()

for(cutoff in c(0.2,0.8)){
  info <- fread("./00reference/hg38-mature-tRNAs-DNA.bed") %>%
    select(ID,length)

  # filter
  row0count.index <- which(as.data.frame(rowSums(count == 0))[,1] > ncol(count)*cutoff)
  if(length(row0count.index) == 0){
    count.sub <- count
  }else{
    count.sub <- count[-row0count.index,]
  }
  
  # 计算tpm
  Len <- as.numeric(info$length[match(rownames(count.sub),info$ID)])
  tpm <- countToTPM(count.sub, Len)
  saveRDS(tpm,paste0("01.tRNA.tpm.scale.filter_",cutoff,".rds"))
  
  #################### log #####################
  info <- fread("./00reference/hg38-mature-tRNAs-DNA.bed") %>% 
    select(V1,V2,V3,V4) %>% mutate(V1 = gsub("chr","",V1))
  colnames(info) <- c('Chr','start','end','ID')

  mat <- log10(tpm + 1)
  mat$ID <- rownames(mat)
  mat.dat <- inner_join(info,mat,by='ID')
  colnames(mat.dat)[[1]] = paste0('#',colnames(mat.dat)[[1]])
  write.table(mat.dat,paste0("./01.tRNA.tpm.scale.log10.filter_",cutoff,".mtx"),row.names=F,sep="\t",quote=F)
}


