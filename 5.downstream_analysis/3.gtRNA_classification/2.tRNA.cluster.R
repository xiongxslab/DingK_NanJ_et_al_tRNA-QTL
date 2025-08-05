# pipeline
# 1. integrated and compared the significance signals between cis and trans tQTL
# 2. checked whether there was an isodecoder-tRNAs nearby trans tQTL region and colocalization analysis
# 3. merge the information
# 4. classification
library(data.table)
library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)

############# cluster1 ###############
df.cis.sign<-fread(paste0("../3.fastqtl_nocut/tQTL.qval_beta.0.1.txt"))
df.trans.sign<-fread(paste0("../4.tensorqtl_trans/01.sign.trans_qtl_pairs.txt"))
df.trna<-read.table(paste0("../../06exp_matrix/02.match_sample/tpm_scale.log10/01.tRNA.tpm.scale.log10.filter_0.8.mtx"))[,1:4]
colnames(df.trna) <- c('chr','start','end','ID')
df.cis.sign.stat<-unique(df.cis.sign$phenotype_id)
df.trans.sign.stat<-unique(df.trans.sign$phenotype_id) 

# SN/NN
tRNA.cisS.transN<-setdiff(df.cis.sign.stat,df.trans.sign.stat) 
tRNA.cisN.transN<-setdiff(df.trna$ID,union(df.cis.sign.stat,df.trans.sign.stat)) 
tRNA.cisN.transS<-setdiff(df.trans.sign.stat,df.cis.sign.stat) 
tRNA.cisS.transS<-intersect(df.trans.sign.stat,df.cis.sign.stat)

df.trna$label<-'NA'
df.trna$label[match(tRNA.cisS.transN,df.trna$ID)]<-'cisS_transN'
df.trna$label[match(tRNA.cisS.transS,df.trna$ID)]<-'cisS_transS'
df.trna$label[match(tRNA.cisN.transS,df.trna$ID)]<-'cisN_transS'
df.trna$label[match(tRNA.cisN.transN,df.trna$ID)]<-'cisN_transN'

cis.sign <- subset(df.cis.sign,phenotype_id %in% tRNA.cisS.transS)
cis.sign.p <- cis.sign %>% group_by(phenotype_id) %>% arrange(-log10(pval_nominal)) %>% slice_tail(n=1)
trans.sign <- subset(df.trans.sign,phenotype_id %in% tRNA.cisS.transS)
trans.sign.p<- trans.sign %>% group_by(phenotype_id) %>% arrange(-log10(pval)) %>% slice_tail(n=1)
tmp<-cis.sign.p[which(cis.sign.p$pval_nominal>trans.sign.p$pval),1]
df.trna$label[match(tmp$phenotype_id,df.trna$ID)]<-'ciss_transS'
tmp<-cis.sign.p[which(cis.sign.p$pval_nominal<trans.sign.p$pval),1]
df.trna$label[match(tmp$phenotype_id,df.trna$ID)]<-'cisS_transs'
write.table(df.trna,paste0('01.tRNA.class1.bed'),col.names = T,row.names = F,sep = '\t',quote=F)

################### trans loci higher and has isodecoder --> coloc ###################
library("coloc")
trna <- fread(paste0("../../06exp_matrix/02.match_sample/tpm_scale.log10/01.tRNA.tpm.scale.log10.filter_0.8.mtx"))[,1:4] %>%
  mutate(aa = apply(as.matrix(ID),1,function(x){unlist(strsplit(x,'-'))[2]}),
         anticodon = apply(as.matrix(ID),1,function(x){unlist(strsplit(x,'-'))[3]}),
         iso = paste(aa,anticodon,sep='-'))

trna.cluster <- fread(paste0('01.tRNA.class1.bed'))
trans.has.iso.trna <- fread(paste0("../4.tensorqtl_trans/01.sign.tRNA.trans_qtl.has.isodecoder.mtx"))$ID
trna.cluster$trans.has.iso <- ifelse(trna.cluster$ID %in% trans.has.iso.trna, 1, 0)
write.table(trna.cluster,paste0('01.tRNA.class2.bed'),col.names = T,row.names = F,sep="\t",quote=F)

############# cluster2 ############### 
library(GenomicRanges)
# merge information of isodecoder
# Determine whether there is an isodecoder at the trans position. 
# If so, what is the pp4?
rm(list=ls())
trna <- fread(paste0("../../06exp_matrix/02.match_sample/tpm_scale.log10/01.tRNA.tpm.scale.log10.filter_0.8.mtx"))[,1:4] %>%
  mutate(aa = apply(as.matrix(ID),1,function(x){unlist(strsplit(x,'-'))[2]}),
         anticodon = apply(as.matrix(ID),1,function(x){unlist(strsplit(x,'-'))[3]}),
         iso = paste(aa,anticodon,sep='-'))

trna.cluster <- fread(paste0('01.tRNA.class2.bed'))
trna.cluster$iso_num = 0
trna.cluster$iso_name = 0
trna.cluster$iso_pp4 = 0

# tRNA cis-tQTL
tqtl.cis <- fread(paste0("../3.fastqtl_nocut/tQTL.sumstat.txt.gz"),
                  col.names=c('phenotype_id','variant_id','distance','ma_samples','ma_count',
                              'maf','pval','b','b_se'))
# tRNA trans-sign-tQTL
tqtl.trans <- fread(paste0("../4.tensorqtl_trans/01.sign.trans_qtl_pairs.txt"))

# coloc result
coloc <- fread(paste0('02.trna.iso.coloc.res'),
               col.names = c('tRNA.trans','tRNA.iso.cis', 'nsnps', 'PP.H0.abf',
                             'PP.H1.abf','PP.H2.abf','PP.H3.abf', 'PP.H4.abf'))

index <- which(trna.cluster$trans.has.iso == 1) # filter trans-significant tRNA
for(i in index){
  trna.i <- trna.cluster[i,ID]
  trna.i.iso <- subset(trna,
                       iso==paste0(unlist(strsplit(trna.i,'-'))[2],'-',unlist(strsplit(trna.i,'-'))[3])&
                         ID!=trna.i) 
  if(nrow(trna.i.iso)==0) next;
  
  # extract one tRNA trans-QTL
  tqtl.trans.i <- subset(tqtl.trans,phenotype_id == trna.i)
  
  # determine whether there is an isodecoder in the trans loci
  site.gr <- GRanges(seqnames = tqtl.trans.i$chr,ranges = IRanges(start = tqtl.trans.i$pos, end = tqtl.trans.i$pos))
  region.gr <- GRanges(seqnames = trna.i.iso$`#Chr`,ranges = IRanges(start = trna.i.iso$start-1e6, end = trna.i.iso$end+1e6))
  
  if(length(findOverlaps(site.gr,region.gr)) != 0) {
    res <- as.data.frame(findOverlaps(site.gr,region.gr))
    trna.cluster$iso_num[i] = length(unique(res$subjectHits))
    coloc.iso = subset(coloc,tRNA.trans == trna.i & tRNA.iso.cis %in% trna.i.iso[unique(res$subjectHits),ID])
    if(nrow(coloc.iso) > 1){
      trna.cluster$iso_name[i] = paste(coloc.iso$tRNA.iso.cis,collapse =',')
      trna.cluster$iso_pp4[i] = paste(coloc.iso$PP.H4.abf,collapse =',')
    } else {
      trna.cluster$iso_name[i] = coloc.iso$tRNA.iso.cis
      trna.cluster$iso_pp4[i] = coloc.iso$PP.H4.abf
    }
  } else next;
}
write.table(trna.cluster, paste0("02.tRNA.class.bed"),col.names = T,row.names = F,sep = '\t',quote=F)


############# cluster3 ############### 
rm(list=ls())
# compare the significant of cis and trans
trna.cluster <- fread(paste0("02.tRNA.class.bed"))%>% as.data.frame()
trna.cluster$`cis.pval` = NA
trna.cluster$`trans.pval` = NA
trna.cluster$`tqtl.higher`='NA'

# top cis-tQTL
cis <- fread(paste0("../3.fastqtl_nocut/tQTL.sumstat.txt.gz"),
             col.names = c('phenotype_id',	'variant_id',	'distance',	'ma_samples',
                           'ma_count',	'maf',	'pval_nominal',	'slope',	'slope_se')) %>%
       group_by(phenotype_id) %>% arrange(pval_nominal) %>% slice_head(n=1)
for(i in 1:nrow(cis)){
  trna <- cis$phenotype_id[i]
  index.in.cluster.df <- match(trna,trna.cluster$ID)
  trans <- fread(paste0("../4.tensorqtl_trans/split.tRNA/",trna,'.trans.txt.gz'),
                 col.names = c("variant_id","phenotype_id","pval","b","b_se","af")) %>%
           arrange(pval) %>% slice_head(n = 1)
  trna.cluster$cis.pval[index.in.cluster.df] = cis$pval_nominal[i]
  trna.cluster$trans.pval[index.in.cluster.df] = trans$pval
  trna.cluster$tqtl.higher[index.in.cluster.df] <- ifelse(cis$pval_nominal[i]<=trans$pval,'cis','trans')
  print(i)
}

# colnames(trna.cluster) <- c("Chr" ,"Start","End","ID","cis.trans","iso_num","iso_name","iso_pp4","transSS","tqtl.higher")
write.table(trna.cluster,paste0("03.tRNA.class.bed"),col.names = T,row.names = F,quote=F,sep="\t")


rm(list=ls())
############# final classification ############### 
# First position: high confidence (1:cis,2:trans,3:cis+trans) #
# Second position: low confidence (1:cis,2:trans,3:misalign) #
# 0: not genetic regulated #
trna.cluster <- fread(paste0("03.tRNA.class.bed"))%>% as.data.frame() %>%
  mutate(hc=0,
         lc=0,
         classification=NA)

## 1.1.high confidence (cis significant: Only cisS --> hc_cis; cisS.transs --> hc_cis + hc_trans) ##
trna.cluster$hc[which(trna.cluster$label == 'cisS_transN')] <- 1
trna.cluster$hc[which(trna.cluster$label == 'cisS_transs')] <- 3


## 1.2.high confidence (Only trans significant & no isodecoder --> hc_trans) ##
trna.cluster$hc[which(trna.cluster$label == 'cisN_transS' & trna.cluster$trans.has.iso == 0)] <- 2


## 1.3.high confidence (cis<trans: ciss.transS & no isodecoder --> hc_cis + hc_trans) ##
trna.cluster$hc[which(trna.cluster$label == 'ciss_transS' & trna.cluster$trans.has.iso == 0)] <- 3


## 2.1.low confidence (cis<trans: has isodecoder & PP4 > 0.5 --> likely misalign; has isodecoder & PP4 < 0.5 --> lc_cis) ##
index = which(trna.cluster$label %in% c('ciss_transS') & trna.cluster$iso_num > 0)
for(i in index){
  trna.cluster$lc[i] <- ifelse(any(as.numeric(unlist(strsplit(trna.cluster$iso_pp4[i],',')))>0.5),3,1)
}

  
## 2.2.low confidence (Only trans siginificant) ##
# has isodecoder, isodecoder's cisN --> likely misalign #
# has isodecoder, isodecoder's cisS & coloc > 0.5 --> likely misalign #
# has isodecoder, isodecoder's cis significant & coloc < 0.5 --> lc_trans #
index = which(trna.cluster$label %in% c('cisN_transS') & trna.cluster$iso_num > 0)
for(i in index){
  iso.name = unlist(strsplit(trna.cluster$iso_name[i],','))
  if(any(subset(trna.cluster,ID %in% iso.name)$label %in% c('cisN_transN','cisN_transS'))){
    trna.cluster$lc[i] = 3
  }else{
    trna.cluster$lc[i] <- ifelse(any(as.numeric(unlist(strsplit(trna.cluster$iso_pp4[i],','))) > 0.5),3,2)
  }
}
trna.cluster$classification <- paste0(trna.cluster$hc,trna.cluster$lc)
print(table(trna.cluster$classification))

write.table(trna.cluster,paste0("04.tRNA.classification.bed"),col.names = T,row.names = F,quote=F,sep="\t")

