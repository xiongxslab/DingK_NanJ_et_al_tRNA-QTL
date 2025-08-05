start = c('pf2')
end=10
pc=3
for(chr in 1:22){
  code <- paste0("/data/dingk/software/fastqtl-master/bin/fastQTL.static --vcf /data/dingk/tRNA-QTL/07genotype/3.split_chr/WGS.chr",
                 chr,".vcf.gz --bed /data/dingk/tRNA-QTL/06exp_matrix/02.match_sample/tpm_scale.log10/01.tRNA.tpm.scale.log10.filter_0.8.mtx.sort.bed.gz --region ",
                 chr," --out output.1M/chr",chr,".out.txt.gz --cov /data/dingk/tRNA-QTL/06exp_matrix/03.PEER/01.covariates.",
                 start,"-",end,".pc1-",pc,".txt --window 1e6 --permute 1000 1>chr",chr,
                 ".out 2>chr",chr,".err")
  write.table(code,"0.submit.sh",col.names = F,row.names = F,sep="\t",quote=F,append = T)
}




