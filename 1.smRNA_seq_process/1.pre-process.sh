for file in ./00data/*.fastq.gz;do
file=${file##*/}
  # trim-galore
  nohup trim_galore -q 20 --small_rna --phred33 --stringency 3 --length 15 -e 0.1 -o ./01trimmed_data ./00data/${file} >> ./01trimmed_data/trim.out
  
  # 1st alignment
  echo ${file} >> ./02bowtie_arti/bowtie.out
  nohup bowtie -x ./00reference/Hsapi38_arti -v 3 --best ./01trimmed_data/${file%.fastq.gz}_trimmed.fq.gz -S -p 8 2> ./02bowtie_arti/${file%.fastq.gz}.out | samtools view -bS | samtools sort - > ./02bowtie_arti/${file%.fastq.gz}.sort.bam
  
  # filter
  Rscript 1.1.tReasure-Bfilter.r ${file%.fastq.gz}.sort.bam
  
  # 2nd alignment
  nohup bowtie -x ./00reference/hg38-tRNA-DNA -n 1 -m10 --best --strata ./03filter_fastq/${file%.fastq.gz}_trimmed_mature.fastq.gz -S -p 8 2> ./04bowtie_trna/${file%.fastq.gz}.out | samtools view -bS | samtools sort - > ./04bowtie_trna/${file%.fastq.gz}.sort.bam
  
  # quantification
  perl 1.2.calculate.tRNA.count.pl ./04bowtie_trna/${file%.fastq.gz}.sort.bam ./05quantification/${file%.fastq.gz}.count.xls
done
