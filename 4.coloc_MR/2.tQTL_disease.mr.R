.libPaths("/data/dingk/R")
options(ieugwasr.api = NULL)
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(TwoSampleMR))
suppressMessages(library(ieugwasr))
suppressMessages(library(plinkbinr))

#### batch ####
# setwd("/data/dingk/tRNA-QTL/10td/2.mr")
# disease.info <- fread("/data/dingk/tRNA-QTL/00reference/GWAS_data/1.filtered_diseases.info.txt") %>%
#   filter(!(classification %in% 'Other')) %>%
#   select('filename','n_cases_EUR','n_controls_EUR','n','type','batch','classification')
# trna.info <- fread("/data/dingk/tRNA-QTL/08tQTL/6.tRNA_cluster/04.tRNA.classification.bed") %>%
#   filter(classification %in% c(10,30))
# for(i in 1:nrow(trna.info)){
#   for(d in disease.info$filename){
#     code <- paste0("Rscript 2.tQTL_disease.mr.R ",d," ",trna.info$ID[i],' ',trna.info$classification[i])
#     write.table(code,'0.submit.sh',col.names = F,row.names = F,sep="\t",quote=F,append = T)
#   }
# }

df.a_make<- function(trna,disease){
  print(disease)
  snp.pos <- fread("/data/dingk/tRNA-QTL/07genotype/2.match_sample/05.WGS.AllChr.match_filterSample_filterSNP_match38.bim",
                   col.names = c('chr', 'SNP', 'family', 'pos', 'A1', 'A2')) %>%
    dplyr::select(SNP,A1,A2)
  cis <- fread(paste0("/data/dingk/tRNA-QTL/08tQTL/3.fastqtl_nocut/tQTL.sumstat.txt.gz"),
               col.names = c('phenotype_id','variant_id','distance','ma_samples','ma_count','maf','pval_nominal','slope','slope_se')) %>%
    filter(phenotype_id == trna & pval_nominal < 1e-3) %>%
    dplyr::select(phenotype_id,variant_id,maf,pval_nominal,slope,slope_se)
  trans <- fread(paste0("/data/dingk/tRNA-QTL/08tQTL/4.tensorqtl_trans/split.tRNA/",trna,".trans.txt.gz"),
                 col.names = c("variant_id","phenotype_id","pval","b","b_se","af")) %>%
    filter(pval < 5e-8/421) %>%
    dplyr::select(phenotype_id,variant_id,af,pval,b,b_se)
  df.t <- rbind(cis,trans,use.names = FALSE) %>%
    inner_join(snp.pos,by=c('variant_id'='SNP'))
  colnames(df.t) <- c('Phenotype','SNP','maf','pval','beta','se','A1','A2')

  info.g<- fread('/data/dingk/tRNA-QTL/00reference/GWAS_data/1.filtered_diseases.info.txt', header=T) %>%
    select('filename','n_cases_EUR','n_controls_EUR','n','type','batch','classification') %>%
    filter(filename == disease)
  if(info.g$batch=='gwas'){
    df.d<- fread(paste0('/data/dingk/tRNA-QTL/00reference/GWAS_data/all.disease/',disease,'.txt.gz'), header=T) %>%
      filter(MAF > 0 & MAF < 1) %>%
      mutate(A1.n = ifelse(MAF <= 0.5, toupper(A1), toupper(A2)),
             A2.n = ifelse(MAF <= 0.5, toupper(A2), toupper(A1)),
             MAF.n = ifelse(MAF <= 0.5, MAF, 1-MAF)) %>%
      select(-A1,-A2,-MAF) %>%
      rename(A1 = A1.n,
             A2 = A2.n,
             MAF = MAF.n) %>%
      select(SNP,MAF,pvalues,beta,varbeta,A1,A2)
  }else{
    df.d<- fread(paste0('/data/dingk/tRNA-QTL/00reference/GWAS_data/ukb/1.summary_data/',disease,'.tsv.gz'), header=T) %>%
      filter(maf > 0 & maf < 1) %>%
      mutate(A1.n = ifelse(maf <= 0.5, toupper(A1), toupper(A2)),
             A2.n = ifelse(maf <= 0.5, toupper(A2), toupper(A1)),
             maf.n = ifelse(maf <= 0.5, maf, 1-maf)) %>%
      select(-A1,-A2,-maf) %>%
      rename(A1 = A1.n,
             A2 = A2.n,
             maf = maf.n) %>%
      select(SNP,maf,pval,beta,se,A1,A2)
  }
  colnames(df.d) <- c('SNP','maf','pval','beta','se','A1','A2')
  
  input<- merge(df.t, df.d, by='SNP', suffixes = c('.tQTL','.GWAS')) %>%
    filter((A1.tQTL==A1.GWAS & A2.tQTL==A2.GWAS)|
             (A2.tQTL==A1.GWAS & A1.tQTL==A2.GWAS)) %>%
    mutate(A1 = A1.GWAS,
           A2 = A2.GWAS,
           beta.tQTL =ifelse(A1.tQTL == A1.GWAS, beta.tQTL, -beta.tQTL)) %>%
    dplyr::select(-A1.tQTL,-A2.tQTL,-A1.GWAS,-A2.GWAS)
  
  write.table(input,paste0('./1.input_files/1.',trna,'.',disease,'.cis&trans.input.res'),col.names = T,row.names = F,quote=F,sep="\t")
}

cis_and_trans_mr <- function(trna,disease){
  input <- fread(paste0('./1.input_files/1.',trna,'.',disease,'.cis&trans.input.res'))
  print("Start to proccess……")
  
  info.g<- fread('/data/slurm/dingk/tRNA-QTL/00reference/GWAS_data/1.filtered_diseases.info.txt', header=T) %>%
    select('filename','n_cases_EUR','n_controls_EUR','n','type','batch','classification') %>%
    filter(filename == disease)
  
  # 挑选工具变量
  input.iv <- input  %>% dplyr::select(SNP,pval.tQTL)
  colnames(input.iv) <- c('rsid','pval')
  b_file = paste0("/data/slurm/dingk/tRNA-QTL/07genotype/2.match_sample/05.WGS.AllChr.match_filterSample_filterSNP_match38")
  iv <- ld_clump_local(input.iv,
                       clump_kb = 100,
                       clump_r2 = 0.001, clump_p = 1,
                       bfile = b_file, plink_bin = get_plink_exe())
  
  # 输出工具变量
  dat <- input %>% filter(SNP %in% iv$rsid) %>%
    mutate(Phenotype.GWAS = disease)
  fout.iv <- paste0("./2.iv/01.",disease,"mr.iv.res")
  if(!file.exists(fout.iv)){
    write.table(t(c("SNP","Phenotype","maf.tQTL","pval.tQTL","beta.tQTL","se.tQTL",
                    "maf.GWAS","pval.GWAS","beta.GWAS","se.GWAS","A1","A2","Phenotype.GWAS")),
                fout.iv,
                sep="\t",col.names = F,row.names = F,quote=F)
  }
  write.table(dat,fout.iv,append = T,col.names = F,row.names = F,sep="\t",quote = F)
  
  if(info.g$batch[1] == 'gwas'){
    if(info.g$type[1]=='binary'){
      tqtl_exp_dat <- format_data(dat, type = "exposure", phenotype_col = "Phenotype", snp_col = "SNP",
                                  beta_col = "beta.tQTL", se_col = "se.tQTL",effect_allele_col = "A1",
                                  other_allele_col = "A2", eaf_col = "maf.tQTL", pval_col = "pval.tQTL")
      gwas_out_dat <- format_data(dat, type = "outcome", phenotype_col = "Phenotype.GWAS", snp_col = "SNP",
                                  beta_col = "beta.GWAS", se_col = "se.GWAS",effect_allele_col = "A1",
                                  other_allele_col = "A2", eaf_col = "maf.GWAS", pval_col = "pval.GWAS",
                                  samplesize = info.g$n[1], ncase = info.g$n_cases_EUR[1],
                                  ncontrol = info.g$n_controls_EUR[1])
    }else{
      tqtl_exp_dat <- format_data(dat, type = "exposure", phenotype_col = "Phenotype", snp_col = "SNP",
                                  beta_col = "beta.tQTL", se_col = "se.tQTL",effect_allele_col = "A1",
                                  other_allele_col = "A2", eaf_col = "maf.tQTL", pval_col = "pval.tQTL")
      gwas_out_dat <- format_data(dat, type = "outcome", phenotype_col = "Phenotype.GWAS", snp_col = "SNP",
                                  beta_col = "beta.GWAS", se_col = "se.GWAS",effect_allele_col = "A1",
                                  other_allele_col = "A2", eaf_col = "maf.GWAS", pval_col = "pval.GWAS")
    }
  }else{
    tqtl_exp_dat <- format_data(dat, type = "exposure", phenotype_col = "Phenotype", snp_col = "SNP",
                                beta_col = "beta.tQTL", se_col = "se.tQTL",effect_allele_col = "A1",
                                other_allele_col = "A2", eaf_col = "maf.tQTL", pval_col = "pval.tQTL")
    gwas_out_dat <- format_data(dat, type = "outcome", phenotype_col = "Phenotype.GWAS", snp_col = "SNP",
                                beta_col = "beta.GWAS", se_col = "se.GWAS",effect_allele_col = "A1",
                                other_allele_col = "A2", eaf_col = "maf.GWAS", pval_col = "pval.GWAS",
                                samplesize = info.g$n[1], ncase = info.g$n_cases_EUR[1],
                                ncontrol = info.g$n_controls_EUR[1])
  }
  
  
  mr_dat <- harmonise_data(exposure_dat = tqtl_exp_dat,
                           outcome_dat = gwas_out_dat) %>% filter(mr_keep == TRUE)
  if(length(mr_dat$SNP) >= 1){
    mr_res <- mr(mr_dat, method_list = c("mr_wald_ratio", "mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
    mr_res <- mr_res %>% dplyr::select(-id.exposure, -id.outcome) %>% dplyr::select(exposure, everything())
    fout <- paste0("./01.tRNA_disease.mr.res")
    if(!file.exists(fout)){
      write.table(t(c("exposure","outcome","method","nsnp","b","se","pval")),
                  fout,
                  quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    }
    write.table(mr_res, fout, append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  }                        
}  

cis_mr <- function(trna,disease){
  input <- fread(paste0("/data/slurm/dingk/tRNA-QTL/10td/1.coloc/1.input_files/1.",disease,".hc.cis.tRNA.tQTL.input.res")) %>%
    filter(Phenotype == trna & pval.tQTL < 1e-3) %>%
    na.omit()
      
  if(nrow(input) > 0){
    info.g<- fread('/data/slurm/dingk/tRNA-QTL/00reference/GWAS_data/1.filtered_diseases.info.txt', header=T) %>%
      select('filename','n_cases_EUR','n_controls_EUR','n','type','batch','classification') %>%
      filter(filename == disease)
    input.iv <- input  %>% dplyr::select(SNP,pval.tQTL)
    colnames(input.iv) <- c('rsid','pval')
    b_file = paste0("/data/slurm/dingk/tRNA-QTL/07genotype/2.match_sample/05.WGS.AllChr.match_filterSample_filterSNP_match38")
    iv <- ld_clump_local(input.iv,
                         clump_kb = 100,
                         clump_r2 = 0.001, clump_p = 1,
                         bfile = b_file, plink_bin = get_plink_exe())
    
    # 输出工具变量
    dat <- input %>% filter(SNP %in% iv$rsid) %>%
      mutate(Phenotype.GWAS = disease)
    fout.iv <- paste0("./2.iv/01.",disease,"mr.iv.res")
    if(!file.exists(fout.iv)){
      write.table(t(c("SNP","Phenotype","maf.tQTL","pval.tQTL","beta.tQTL","se.tQTL",
                      "maf.GWAS","pval.GWAS","beta.GWAS","se.GWAS","A1","A2","Phenotype.GWAS")),
                  fout.iv,
                  sep="\t",col.names = F,row.names = F,quote=F)
    }
    write.table(dat,fout.iv,append = T,col.names = F,row.names = F,sep="\t",quote = F)
    
    if(info.g$batch[1] == 'gwas'){
      if(info.g$type[1]=='binary'){
        tqtl_exp_dat <- format_data(dat, type = "exposure", phenotype_col = "Phenotype", snp_col = "SNP",
                                    beta_col = "beta.tQTL", se_col = "se.tQTL",effect_allele_col = "A1",
                                    other_allele_col = "A2", eaf_col = "maf.tQTL", pval_col = "pval.tQTL")
        gwas_out_dat <- format_data(dat, type = "outcome", phenotype_col = "Phenotype.GWAS", snp_col = "SNP",
                                    beta_col = "beta.GWAS", se_col = "se.GWAS",effect_allele_col = "A1",
                                    other_allele_col = "A2", eaf_col = "maf.GWAS", pval_col = "pval.GWAS",
                                    samplesize = info.g$n[1], ncase = info.g$n_cases_EUR[1],
                                    ncontrol = info.g$n_controls_EUR[1])
      }else{
        tqtl_exp_dat <- format_data(dat, type = "exposure", phenotype_col = "Phenotype", snp_col = "SNP",
                                    beta_col = "beta.tQTL", se_col = "se.tQTL",effect_allele_col = "A1",
                                    other_allele_col = "A2", eaf_col = "maf.tQTL", pval_col = "pval.tQTL")
        gwas_out_dat <- format_data(dat, type = "outcome", phenotype_col = "Phenotype.GWAS", snp_col = "SNP",
                                    beta_col = "beta.GWAS", se_col = "se.GWAS",effect_allele_col = "A1",
                                    other_allele_col = "A2", eaf_col = "maf.GWAS", pval_col = "pval.GWAS")
      }
    }else{
      tqtl_exp_dat <- format_data(dat, type = "exposure", phenotype_col = "Phenotype", snp_col = "SNP",
                                  beta_col = "beta.tQTL", se_col = "se.tQTL",effect_allele_col = "A1",
                                  other_allele_col = "A2", eaf_col = "maf.tQTL", pval_col = "pval.tQTL")
      gwas_out_dat <- format_data(dat, type = "outcome", phenotype_col = "Phenotype.GWAS", snp_col = "SNP",
                                  beta_col = "beta.GWAS", se_col = "se.GWAS",effect_allele_col = "A1",
                                  other_allele_col = "A2", eaf_col = "maf.GWAS", pval_col = "pval.GWAS",
                                  samplesize = info.g$n[1], ncase = info.g$n_cases_EUR[1],
                                  ncontrol = info.g$n_controls_EUR[1])
    }
    
    mr_dat <- harmonise_data(exposure_dat = tqtl_exp_dat,
                             outcome_dat = gwas_out_dat) %>% filter(mr_keep == TRUE)
    if(length(mr_dat$SNP) >= 1){
      mr_res <- mr(mr_dat, method_list = c("mr_wald_ratio", "mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
      mr_res <- mr_res %>% dplyr::select(-id.exposure, -id.outcome) %>% dplyr::select(exposure, everything())
      fout <- paste0("./01.tRNA_disease.mr.res")
      if(!file.exists(fout)){
        write.table(t(c("exposure","outcome","method","nsnp","b","se","pval")),
                    fout,
                    quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
      }
      write.table(mr_res, fout, append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    }     
  }else{
    write.table(t(c(disease,trna,"has no significant iv for MR")),
                './01.mr.tQTL_GWAS.noSig.res',
                sep="\t",quote=F,col.names = F,row.names = F,append = T)
  }
}


# Alzheimer_disease tRNA-Glu-TTC-4-1 10
args<- commandArgs(T)
classification=args[3]
trna=args[2]
disease=args[1]

setwd("/data/slurm/dingk/tRNA-QTL/10td/2.mr")
if(classification == 10){
  df.a_make(trna,disease)
  cis_mr(trna,disease)
}else{
    df.a_make(trna,disease)
    cis_and_trans_mr(trna,disease)
}
system(paste0("echo ",disease,' ',trna,' ',classification, ' >> 0.submit.done'))


