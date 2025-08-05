# format transform #
cutoffs=("0.2 0.8")
for cutoff in ${cutoffs};do
awk 'NR == 1; NR > 1 {print $0 | "sort -k1,1 -k2n,2"}' 01.tRNA.tpm.scale.log10.filter_${cutoff}.mtx > 01.tRNA.tpm.scale.log10.filter_${cutoff}.mtx.sort.bed
bgzip 01.tRNA.tpm.scale.log10.filter_${cutoff}.mtx.sort.bed  && tabix -p bed 01.tRNA.tpm.scale.log10.filter_${cutoff}.mtx.sort.bed.gz
done