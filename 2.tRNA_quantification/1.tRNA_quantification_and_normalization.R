library(data.table)
library(dplyr)
library(reshape2)

countToTPM <- function(counts, effLen)
{
  rate <- counts / effLen
  denom <- sum(rate)
  rate/denom*1e6
}

# convert the matrix
files<-list.files('./05quantification')
df.merge = data.frame()
for(file in files){
 df = read.table(paste0('./05quantification/',file),head=T)
 if(nrow(df)==0) next;
 df$samp = unlist(strsplit(file,'[.]'))[1]
 df.merge = rbind(df.merge,df)
}

mat = acast(df.merge %>% mutate(tRNA = gsub('Homo_sapiens_','',tRNA)), tRNA~samp, value.var="count")
mat[is.na(mat)] = 0
saveRDS(mat,'./06exp_matrix/00.tRNA.count.mat.rds')


info <- fread("./00reference/hg38-mature-tRNAs-DNA.bed") %>%
  dplyr::select(ID,length)

# tRNAs that were not expressed in over 80% of the samples were excluded.
row0count.index <- which(as.data.frame(rowSums(mat == 0))[,1] > ncol(mat)*0.8)
if(length(row0count.index) == 0){
  count <- mat
}else{
  count <- mat[-row0count.index,]
}  
saveRDS(count,"./06exp_matrix/01.tRNA.count.filter_0.8.rds")

# transform TPM
Len <- as.numeric(info$length[match(rownames(count),info$ID)])
tpm <- apply(count, 2, function(x) {
  countToTPM(x, length)
})
saveRDS(tpm,"./06exp_matrix/02.tRNA.tpm.filter_0.8.rds")

