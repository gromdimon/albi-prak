#!/bin/bash

input_dir="/group/albi-praktikum2023/SNPs"  # Directory containing all compressed VCF files
output_dir="/group/albi-praktikum2023/analysis/gruppe_3/aufgabe04"

# Assuming Individual 3's VCF is named individual_3.vcf.gz and others follow a similar naming pattern
vcf_files=($(ls ${input_dir}/gruppe*.vcf.gz | grep -v "gruppe3.vcf.gz"))
bcftools isec -c none -C  "${input_dir}/gruppe3.vcf.gz" "${vcf_files[@]}" -Oz -p "${output_dir}"
