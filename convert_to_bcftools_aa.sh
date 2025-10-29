#!/bin/bash

# Usage: ./convert_to_bcftools_aa.sh input.tsv [output_prefix]

input=$1
output_prefix=${2:-ancestral_annotation}

# Convert halAncestralAllele output to bcftools annotate format
# Columns: chrom start(0-based) end refAllele usedAncestor ancestralAllele evidence
awk 'BEGIN {OFS="\t"} 
{
    # Skip header if present
    if (NR==1 && $1 ~ /^#/) next
    
    # Convert 0-based to 1-based position
    pos = $2 + 1
    
    # Use ancestral allele from column 6
    # Filter out "N" (unknown) ancestral alleles
    if ($6 != "N") {
        print $1, pos, $6
    }
}' $input > ${output_prefix}.tsv

# Sort the file (required for tabix)
sort -k1,1 -k2,2n ${output_prefix}.tsv > ${output_prefix}_sorted.tsv
mv ${output_prefix}_sorted.tsv ${output_prefix}.tsv

# Compress with bgzip
bgzip -f ${output_prefix}.tsv

# Index with tabix (chromosome in column 1, position in column 2)
tabix -f -s 1 -b 2 -e 2 ${output_prefix}.tsv.gz

echo "Created bgzip compressed and tabix indexed file: ${output_prefix}.tsv.gz"
echo "Use with:"
echo "bcftools annotate -a ${output_prefix}.tsv.gz -c CHROM,POS,AA -h <(echo '##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral allele\">') -Ob -o output.bcf input.vcf.gz"