# 1.make bed #
library(GenomicRanges)
library(rtracklayer)
library(Rsamtools)
trna.info <- fread("hg38-mature-tRNAs-DNA.bed") %>%
  select(V4,V6) %>%
  dplyr::rename(ID = V4, strand = V6)
for(para in c('lcl.eur','brain.eur','lcl.yri')){
  for(win in c(200)){
    qtl <- fread(paste0("04.",para,".cis-independent-tQTL.hc.txt")) %>%
      left_join(trna.info,by=c('phenotype_id'='ID')) %>%
      select(chr,pos,strand) %>%
      distinct() %>%
      mutate(start = pos - win,
             end = pos + win) %>%
      select(-pos) %>%
      arrange(chr,start)
    gr <- GRanges(seqnames = paste0("chr",qtl$chr),
                  ranges = IRanges(start = qtl$start, end = qtl$end),
                  strand = qtl$strand)
    qtl.merged <- data.frame(chr = as.character(seqnames(gr)),
                             start = start(gr),
                             end = end(gr),
                             strand = as.character(strand(gr)))
    write.table(gr,paste0("01.",para,".win",win,"bp.bed"),
                col.names = F,row.names = F,sep="\t",quote=F)
    
    genome_fa <- FaFile("./00reference/hg38.fa")
    extracted_seqs <- getSeq(genome_fa, gr)
    names(extracted_seqs) <- paste0(qtl.merged$chr, ":", qtl.merged$start, "-", qtl.merged$end, "(", qtl.merged$strand, ")") # 为序列添加名称（格式：chr:start-end(strand)）
    writeXStringSet(extracted_seqs, filepath = paste0("02.",para,".win",win,"bp.sequences.fasta"), format = "fasta")
  }
}

for(para in c('lcl.eur','brain.eur','lcl.yri')){
  for(win in c(200)){
    system(paste0("~/software/meme-5.5.7/bin/sea -verbosity 1 --text ",
                  "--m ./00reference/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt ",
                  "--p 02.",para,".win",win,"bp.sequences.fasta > ",
                  "03.",para,".win",win,"bp.SEA.output"))
  }
}

for(para in c('lcl.eur','brain.eur','lcl.yri')){
  for(win in c(200)){
    system(paste0("~/software/meme-5.5.7/bin/sea -verbosity 1 -o ./0.output_dir/",
                  para,".win",win,"bp.SEA ",
                  "--m ./tRNA-QTL/00reference/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt ",
                  "--p 02.",para,".win",win,"bp.sequences.fasta"))
  }
}

# 2.merge #
res <- data.frame()
for(para in c('lcl.eur','brain.eur','lcl.yri')){
  for(win in c(200)){
    tmp <- fread(paste0("03.",para,".win",win,"bp.SEA.output")) %>%
      filter(QVALUE < 0.1 & log2(ENR_RATIO) > 1) %>%
      mutate(label = para,
             windows = win)
    res <- rbind(res,tmp)
  }
}
write.table(res,"04.All.merge.sign.enrich.res",col.names = T,row.names = F,sep = "\t",quote=F)
