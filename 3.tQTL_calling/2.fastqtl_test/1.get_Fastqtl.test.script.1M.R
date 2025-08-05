setwd("/data/dingk/tRNA-QTL/08tQTL/1.fastqtl_test")
starts = c('pf1','pf2')
pc=5
for(start in starts){
  if(start == 'pf1'){
    end.list <- c(1:20)
  }else if(start == 'pf2'){
    end.list <- c(2:20)
  }else{
    end.list <- c(1:2,4:20)
  }

  for(end in end.list){
    for(chr in 1:22){
      code <- paste0("/data/dingk/software/fastqtl-master/bin/fastQTL.static --vcf /data/dingk/tRNA-QTL/07genotype/3.split_chr/WGS.chr",
                     chr,".vcf.gz --bed /data/dingk/tRNA-QTL/06exp_matrix/02.match_sample/tpm_scale.log10/01.tRNA.tpm.scale.log10.filter_0.8.mtx.sort.bed.gz --region ",
                     chr," --out output.1M/chr",chr,".",start,
                     "-",end,".pc1-",pc,".out.txt.gz --cov /data/dingk/tRNA-QTL/06exp_matrix/03.PEER/01.covariates.",
                     start,"-",end,".pc1-",pc,".txt --window 1e6 --threshold 0.0001")
      write.table(code,"0.submit.sh",col.names = F,row.names = F,sep="\t",quote=F,append = T)
    }
  }
}



