# Ribo-seq SE: trim_galore → tRNA/rRNA filter → hisat2

# Build contaminant filter index (tRNA + rRNA)
cat ./reference/hg38-mature-tRNAs-DNA.fa ./reference/ribo_seq_refpanel/human_rRNA.fa > ./reference/contaminant_filter.fa
bowtie2-build ./reference/contaminant_filter.fa ./reference/contaminant_filter

for FQ in ./raw_data/*.fq.gz; do
  sample=$(basename "$FQ" .fq.gz)

  # trim_galore
  trim_galore -o ./01trim_galore "$FQ"

  # filter tRNA/rRNA
  bowtie2 -p 8 --very-sensitive -x ./reference/contaminant_filter \
    -U ./01trim_galore/"$sample"_trimmed.fq.gz 2>/dev/null | \
    samtools view -@ 4 -f 4 -bS | \
    samtools fastq -@ 4 -0 ./02ribo_map_filter/"$sample".clean.fq.gz -s /dev/null

  # hisat2
  hisat2 -p 8 -x ./reference/hg38 \
    -U ./02ribo_map_filter/"$sample".clean.fq.gz \
    2> ./03ribo_map_hg38/"$sample".log | \
    samtools sort -@ 4 -o ./03ribo_map_hg38/"$sample".bam -
  samtools index ./03ribo_map_hg38/"$sample".bam
done
