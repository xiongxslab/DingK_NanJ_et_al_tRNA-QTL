################### trans loci has isodecoder --> coloc ###################
library(coloc)
# tRNA cis-tQTL
tqtl.cis <- fread(paste0("../3.fastqtl_nocut/tQTL.sumstat.txt.gz"),
                  col.names=c('phenotype_id','variant_id','distance','ma_samples','ma_count',
                              'maf','pval','b','b_se'))
# tRNA trans-tQTL
tqtl.trans <- fread(paste0("../4.tensorqtl_trans/02.sign.trans-qtl.has.iso.all.out.trans_qtl_pairs.csv")) %>%
  select(-'__index_level_0__')
trna.size = 535
for(trna.i in trans.has.iso.trna){
  trna.i.iso <- subset(trna,
                       iso==paste0(unlist(strsplit(trna.i,'-'))[2],'-',unlist(strsplit(trna.i,'-'))[3])&
                         ID!=trna.i)
  if(nrow(trna.i.iso)==0) next;
  
  # extract one tRNA trans-QTL
  tqtl.trans.i <- subset(tqtl.trans,phenotype_id == trna.i) %>%
    select(phenotype_id,variant_id,pval,b,b_se,af)
  
  # colocalization analysis
  for(j in 1:nrow(trna.i.iso)){
    # extract tRNA-isodecoder cis-QTL
    tqtl.cis.i.j <- subset(tqtl.cis,phenotype_id == trna.i.iso[j,ID]) %>%
      select(phenotype_id,variant_id,pval,b,b_se,maf) %>% na.omit()
    
    meta <- inner_join(tqtl.trans.i,tqtl.cis.i.j,by=c('variant_id'),suffix=c('.trans','.cis')) %>%
      group_by(variant_id) %>% arrange(pval.trans) %>% slice_head(n=1)
    if(nrow(meta)==0){
      next;
    }else{
      coloc_res <- coloc.abf(dataset1 = list(snp = meta$variant_id, beta = meta$b.trans, varbeta = meta$b_se.trans^2,
                                             MAF = meta$af, N = trna.size, type = "quant"),
                             dataset2 = list(snp = meta$variant_id, beta = meta$b.cis, varbeta = meta$b_se.cis^2,
                                             MAF = meta$maf, N = trna.size, type = "quant"))
      res <- data.frame(tRNA.trans=trna.i,
                        tRNA.iso.cis=trna.i.iso[j,ID],
                        t(as.data.frame(coloc_res$summary)))
      write.table(res, file=paste0('./02.trna.iso.coloc.res'), sep='\t', append=T, col.names = F, row.names = F, quote=F)
    }
  }
}

