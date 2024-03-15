#!/bin/bash

# Directories
input_vcf_dir="/group/albi-praktikum2023/SNPs"
output_dir_base="/group/albi-praktikum2023/analysis/gruppe_3/aufgabe04"

# Ensure the base output directory exists
mkdir -p "$output_dir_base"

# Initial count of unique SNPs in Individual 3
initial_unique=$(bcftools view -H "${input_vcf_dir}/gruppe3.vcf" | wc -l)
echo "Initial unique SNPs in Individual 3: $initial_unique"

# Define the list of groups to compare against (excluding group 3 itself if listed)
groups=(1 2 4 5 6 7 8 9 10 11 12)

for i in "${groups[@]}"; do
  # Specific output directory for this comparison
  output_dir="${output_dir_base}/unique_to_3_vs_${i}"
  mkdir -p "$output_dir"

  # Find unique SNPs for Individual 3 compared to Individual i
  bcftools isec -c none -C "${input_vcf_dir}/gruppe3.vcf.gz" "${input_vcf_dir}/gruppe${i}.vcf.gz" -Oz -p "$output_dir"

  # Count the number of unique SNPs for this comparison
  # Assuming the unique SNPs are in the 0000.vcf.gz file created by bcftools isec
  if [[ -f "${output_dir}/0000.vcf.gz" ]]; then
    unique_snps=$(zcat "${output_dir}/0000.vcf.gz" | grep -v "^#" | wc -l)
  else
    unique_snps=0
  fi

  echo "Unique SNPs after including gruppe $i: $unique_snps"
  decrease=$(echo "scale=2; 100 * ($initial_unique - $unique_snps) / $initial_unique" | bc)
  echo "Decrease in unique SNPs: ${decrease}%"
done
