# Kallisto MANE quantification: Ribo (SE, 26-32nt) + RNA (PE)
# Ribo-seq: --single -l 29 -s 2
# RNA-seq: PE, auto fragment length

for FQ in ./02ribo_map_filter/*.clean.fq.gz; do
  sample=$(basename "$FQ" .clean.fq.gz)
  kallisto quant -i ./reference/MANE_transcripts.idx \
    -o ./04ribo_kallisto_MANE/"$sample" \
    -b 100 --single -l 29 -s 2 "$FQ"
done

for R1 in ./results/*/*_val_1.fq.gz; do
  sample=$(basename "$R1" _val_1.fq.gz)
  R2="${R1/_val_1.fq.gz/_val_2.fq.gz}"
  kallisto quant -i ./reference/MANE_transcripts.idx \
    -o ./04rna_kallisto_MANE/"$sample" \
    -b 100 "$R1" "$R2"
done
