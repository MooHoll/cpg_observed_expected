#!/bin/bash
# -------------------------------------------------------------------
# Making a fasta files of certain annotations (exons, introns, genes)
# -------------------------------------------------------------------

# To run this script: bash this_script.sh species.gff3 species.fa

# Depends on bedtools

GFF=$1
FASTA=$2
BASE_GFF=$(basename ${GFF} ".gff3")

start=`date +%s`
# -------------------------------------------------------------------

echo "making bed files for genes, exons and introns"
# Make a bed of exons
cat ${GFF} | \
awk 'BEGIN{OFS="\t";} $3=="exon" {print $1,$4-1,$5}' | \
sortBed | \
mergeBed -i - > ${BASE_GFF}_exon.bed
# Make a bed of introns
cat ${GFF} | \
awk 'BEGIN{OFS="\t";} $3=="gene" {print $1,$4-1,$5}' | \
sortBed | \
subtractBed -a stdin -b ${BASE_GFF}_exon.bed > ${BASE_GFF}_intron.bed
# Make a bed of genes
cat ${GFF} | \
awk 'BEGIN{OFS="\t";} $3=="gene" {print $1,$4-1,$5}' | \
sortBed | \
mergeBed -i - > ${BASE_GFF}_genes.bed

# -------------------------------------------------------------------

# Take all of the bed files you've just made and get the corresponsing fasta files
echo "now making fasta files from the bed files"

for FILE in $(ls ${BASE_GFF}*.bed)
do
    BASE_BED=$(basename ${FILE} ".bed")
    bedtools getfasta -fi ${FASTA} -bed ${FILE} > ${BASE_BED}.fa
done

# -------------------------------------------------------------------
echo "remove files not needed later"
rm ${BASE_GFF}*.bed
rm ${FASTA}.fai


end=`date +%s`
runtime=$((((end-start)/60)))
echo "runtime:"$runtime"mins"