for gff in *.gff.gz; do
    zcat $gff | awk -F "\t" -v OFS="\t" '$3=="protein_coding_gene"{split($9, a, /[=;]/); print$1,$4,$5, a[2]}' | sort -k 1,1 -k2,2n > ${gff%.gff.gz}.bed
done
