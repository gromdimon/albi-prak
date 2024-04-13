#!/bin/bash

# Specify the path to the VCF file containing SNP data for the X chromosome.
vcf_path="/group/albi-praktikum2023/analysis/gruppe_3/aufgabe03/chrX_filtered.vcf"

# Count the number of heterozygous SNPs on the X chromosome.
# This command filters for entries where the genotype is either 0/1 or 1/0.
heterozygous_snps=$(/group/albi-praktikum2023/software/bcftools view -H $vcf_path | awk '{if ($10 ~ "0/1" || $10 ~ "1/0") print $0}' | wc -l)
echo "Number of heterozygous SNPs on the X chromosome: $heterozygous_snps"

# Count the total number of SNPs on the X chromosome.
total_snps=$(/group/albi-praktikum2023/software/bcftools view -H $vcf_path | wc -l)
echo "Total number of SNPs on the X chromosome: $total_snps"

# Calculate the ratio of heterozygous SNPs to total SNPs on the X chromosome.
ratio=$(echo "scale=2; $heterozygous_snps / $total_snps" | bc)
echo "Ratio of heterozygous to total SNPs: $ratio"

# The script can also check if the ratio is higher than a cutoff value, for example 0.3, to determine gender or other attributes.
cutoff=0.3
if (( $(echo "$ratio > $cutoff" | bc -l) )); then
    echo "The ratio of heterozygous to total SNPs is higher than the cutoff ($cutoff). The individual is likely female."
else
    echo "The ratio of heterozygous to total SNPs is not higher than the cutoff ($cutoff). The individual is likely male."
fi
