.libPaths("/data/dingk/R")
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(coloc))
suppressMessages(library(stringr))

# #### batch ####
# disease.info <- fread("/data/dingk/tRNA-QTL/00reference/GWAS_data/1.filtered_diseases.info.txt") %>%
#   select('filename','n_cases_EUR','n_controls_EUR','n','type','batch','classification')
# for(d in disease.info$filename){
#   code <- paste0("Rscript 1.tQTL_disease.coloc.R ",d)
#   write.table(code,'0.submit.sh',col.names = F,row.names = F,sep="\t",quote=F,append = T)
# }

### 1. make molQTL files
setwd("/data/dingk/tRNA-QTL/10td/1.coloc")
df.a_make<- function(disease){
  print(disease)
  df.t <- fread(paste0("00.hc.cis.tRNA.tQTL.sumstat.txt")) %>%
    select(-distance,-ma_samples,-ma_count,-chr,-pos)
  colnames(df.t) <- c('Phenotype','SNP','maf','pval','beta','se','A1','A2')
  
  info.g<- fread('/data/dingk/tRNA-QTL/00reference/GWAS_data/1.filtered_diseases.info.txt', header=T) %>%
    select('filename','n_cases_EUR','n_controls_EUR','n','type','batch','classification') %>%
    filter(filename == disease)
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
  colnames(df.d) <- c('SNP','maf','pval','beta','se','A1','A2')
  
  input<- merge(df.t, df.d, by='SNP', suffixes = c('.tQTL','.GWAS')) %>%
    filter((A1.tQTL==A1.GWAS & A2.tQTL==A2.GWAS)|
             (A2.tQTL==A1.GWAS & A1.tQTL==A2.GWAS)) %>%
    mutate(A1 = A1.GWAS,
           A2 = A2.GWAS,
           beta.tQTL =ifelse(A1.tQTL == A1.GWAS, beta.tQTL, -beta.tQTL)) %>%
    dplyr::select(-A1.tQTL,-A2.tQTL,-A1.GWAS,-A2.GWAS)
  
  write.table(input,paste0('./1.input_files/1.',disease,'.hc.cis.tRNA.tQTL.input.res'),col.names = T,row.names = F,quote=F,sep="\t")
}

coloc<- function(disease){
  input<- fread(paste0('./1.input_files/1.',disease,'.hc.cis.tRNA.tQTL.input.res'))
  if(nrow(input)==0){
    write.table(t(c(disease,'all genes and GWAS', 'have no overlap loci')), 
                file='./01.coloc.tQTL_GWAS.noSig.res',
                append=T, row.names = F,col.names = F)
  }else{
    info.g<- fread('/data/dingk/tRNA-QTL/00reference/GWAS_data/1.filtered_diseases.info.txt', header=T) %>%
      select('filename','n_cases_EUR','n_controls_EUR','n','type','batch','classification') %>%
      filter(filename == disease)
    trna.size = 535
    type.p=info.g$type[1]
    Ncase=info.g$n_cases_EUR[1]
    Ncontrol=info.g$n_controls_EUR[1]
    N.GWAS=info.g$n[1]

    if(any(input$pval.GWAS<1e-4)>0){
      print(paste(disease,'colocalization analysis is starting...'))
      ##select the significant genes for coloc
      for (phen in unique(subset(input,pval.GWAS<1e-4)$Phenotype)){
        input1<- input %>%
          filter(Phenotype %in% phen) %>%
          na.omit()
        
        if(any(input1$pval.GWAS<1e-4)){
          print(phen)
            if(type.p=='quantitative'){
            result<- coloc.abf(dataset1=list(snp=input1$SNP, beta=input1$beta.tQTL, varbeta=input1$se.tQTL^2,
                                             MAF=input1$maf.tQTL, N=trna.size, type='quant'),
                               dataset2=list(snp=input1$SNP, beta=input1$beta.GWAS, varbeta=input1$se.GWAS^2,
                                             MAF=input1$maf.GWAS, N=N.GWAS, type='quant'))
            }else{
            result<- coloc.abf(dataset1=list(snp=input1$SNP, beta=input1$beta.tQTL, varbeta=input1$se.tQTL^2, 
                                             MAF=input1$maf.tQTL, N=trna.size, type='quant'),
                               dataset2=list(snp=input1$SNP, s=Ncase/N.GWAS, 
                                             beta=input1$beta.GWAS, varbeta=input1$se.GWAS^2,
                                             MAF=input1$maf.GWAS, N=N.GWAS, type='cc'))
            }          
          leadvariance = result$results[which.max(result$results$SNP.PP.H4),1]
          write.table(cbind(disease, phen, t(result$summary),leadvariance ,
                            input1$pval.tQTL[match(leadvariance,input1$SNP)],input1$pval.GWAS[match(leadvariance,input1$SNP)]), 
                      file=paste0('./01.hc.cis.tRNA.tQTL.coloc.res'),sep='\t',
                      quote=F, col.names=F, row.names=F,append=T)
          print(paste(disease,'colocalization analysis has done...'))

      }else{
          # print(paste(pop,disease,'has no significant pairs for colocalization analysis...'))
          write.table(t(c(disease,phen,"has no significant pairs for colocalization")),
                      './01.coloc.tQTL_GWAS.noSig.res',
                      sep="\t",quote=F,col.names = F,row.names = F,append = T)
        }
      }
    }else{
        write.table(t(c(disease,'all genes and GWAS', 'have no significant overlap loci')), 
                    file='./01.coloc.tQTL_GWAS.noSig.res',
                    append=T, row.names = F,col.names = F)
      }
  }
}
args<- commandArgs(T)
disease=args[1]
if(!(file.exists(paste0('./01.hc.cis.tRNA.tQTL.coloc.res')))){
  write.table(t(c("disease","tRNA","nsnps","PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf",
                  "leadvariance","pval.tQTL","pval.GWAS")),
              paste0('./01.hc.cis.tRNA.tQTL.coloc.res'),
              sep="\t",quote=F,col.names = F,row.names = F)
}
if(!(file.exists(paste0('./1.input_files/1.',disease,'.hc.cis.tRNA.tQTL.input.res')))){
  df.a_make(args[1])
  coloc(args[1])
}else{
  coloc(args[1])
}
