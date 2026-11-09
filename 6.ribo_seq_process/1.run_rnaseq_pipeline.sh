# RNA-seq PE: trim_galore → hisat2 → featureCounts
for R1 in ./raw_data/*_R1_001.fastq.gz; do
  sample=$(basename "$R1" _R1_001.fastq.gz)
  R2="${R1/_R1_001.fastq.gz/_R2_001.fastq.gz}"
  [ ! -f "$R2" ] && continue

  # trim_galore
  trim_galore --paired --phred33 -o ./results/"$sample" "$R1" "$R2"

  # hisat2
  hisat2 -p 8 -x ./reference/hg38 \
    -1 ./results/"$sample"/*_val_1.fq.gz -2 ./results/"$sample"/*_val_2.fq.gz \
    2> ./results/"$sample"/hisat2.log | \
    samtools sort -@ 4 -o ./results/"$sample"/"$sample".bam -
  samtools index ./results/"$sample"/"$sample".bam
done

# featureCounts
featureCounts -T 8 -p -t exon -g gene_id \
  -a ./reference/Homo_sapiens.GRCh38.110.gtf \
  -o gene_counts.txt ./results/*/*.bam
