df = read.table("PEER_covariates.txt",head=T)
  rownames(df) = df$ID
  df = df[,-1]
# start from 1
  for(i in 1:50){
   df.cut = df[1:i,]
   write.table(data.frame(ID=rownames(df.cut),df.cut),paste0("PEER_covariates.pf1-",i,".txt"),sep="\t",quote=F,row.names=F)
  }
 if(0){
  # start from 2
  for(i in 2:20){
   df.cut = df[2:i,]
   write.table(data.frame(ID=rownames(df.cut),df.cut),paste0("PEER_covariates.pf2-",i,".txt"),sep="\t",quote=F,row.names=F)
  }
 }


