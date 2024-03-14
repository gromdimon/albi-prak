# Algorithmische Bioinformatik: Praktikum

The practical part of the course "Algorithmische Bioinformatik" (Algorithmic Bioinformatics) at the Freie Universität Berlin.

![header](assets/header.png)

## Table of Contents

- [Introduction](#introduction)
- [Tools we used](#tools-we-used)
- [Day1: Quality Control with FastQC](#day1-quality-control-with-fastqc)
  - [Step 1: Run FastQC](#step-1-run-fastqc)
- [Day2: Mapping Workflow](#day2-mapping-workflow)
    - [Step 1: Align the Reads](#step-1-align-the-reads)
    - [Step 2: Process Alignment Output](#step-2-process-alignment-output)
    - [Step 3: Evaluate Mapping Quality](#step-3-evaluate-mapping-quality)
    - [Step 4: Visualize Alignments in IGV](#step-4-visualize-alignments-in-igv)
    - [Step 5: Calculate Coverage for Each Chromosome](#step-5-calculate-coverage-for-each-chromosome)
- [Day3: Variant Calling](#day3-variant-calling)
    - [Step 1: Add Read Groups](#step-1-add-read-groups)
    - [Step 2: Splitting bam file by chromosome](#step-2-splitting-bam-file-by-chromosome)
    - [Step 3: Variant Calling](#step-3-variant-calling)
    - [Step 4: Merge VCFs](#step-4-merge-vcfs)
- [Additional Information](#additional-information)
- [Supervisors](#supervisors)
- [Contributors](#contributors)

## Introduction

In this practice we've learned how to ...

## Tools we used

- BWA (Burrows-Wheeler Aligner)
- SAMtools
- BEDTools
- IGV (Integrative Genomics Viewer)

## Day1: Quality Control with FastQC

### Step 1: Run FastQC

Since the FastQC was already installed on the server, we could run it using the following command, applying the fastq function on 2 input files, in particularly the parts of them with fast.gz extension:

```sh
/group/albi-praktikum2023/software/fastqc -o /group/albi-praktikum2023/analysis/gruppe_3/aufgabe01 --noextract -f fastq /group/albi-praktikum2023/data/1000-Genome-Project/gruppe3/* .fastq.gz
```
![resultsFastQC]()

## Day2: Mapping Workflow

### Step 1: Align the Reads

Since the reference genome was already indexed, we could start align the sequencing reads to the reference genome using BWA MEM. The following command was used to align small subset of reads:

```sh
/group/albi-praktikum2023/software/bwa-mem2 mem -t 50 /group/albi-praktikum2023/data/NCBI-reference/GRCh38.fna /group/albi-praktikum2023/analysis/gruppe_3/aufgabe02/test_gruppe3_1.fastq.gz > /group/albi-praktikum2023/analysis/gruppe_3/aufgabe02/test_output.sam
```

### Step 2: Process Alignment Output

Convert the SAM file to BAM, sort, and index the aligned reads:

```sh
samtools view -S -b aligned_reads.sam > aligned_reads.bam
samtools sort aligned_reads.bam -o aligned_reads_sorted.bam
samtools index aligned_reads_sorted.bam
```

Thus, the filestructure of the folder should look like this:

```sh
hrad04@compute04:/group/albi-praktikum2023/analysis/gruppe_3/aufgabe02$ ls -a
.  .. test_gruppe3_1.fastq.gz  test_output.bam  test_output.sam  test_output_sorted.bam  test_output_sorted.bam.bai
```

### Step 3: Evaluate Mapping Quality

Assess the quality of the mapping using SAMtools:

```sh
samtools flagstat -o json aligned_reads_sorted.bam
```

This was the output of the command:

```json
{
"QC-passed reads": {
  "total": 733319951,
  "primary": 727401006,
  "secondary": 0,
  "supplementary": 5918945,
  "duplicates": 0,
  "primary duplicates": 0,
  "mapped": 731726121,
  "mapped %": 99.78,
  "primary mapped": 725807176,
  "primary mapped %": 99.78,
  "paired in sequencing": 727401006,
  "read1": 363700503,
  "read2": 363700503,
  "properly paired": 709068384,
  "properly paired %": 97.48,
  "with itself and mate mapped": 724443610,
  "singletons": 1363566,
  "singletons %": 0.19,
  "with mate mapped to a different chr": 11264042,
  "with mate mapped to a different chr (mapQ >= 5)": 5927328
 },
 "QC-failed reads": {
  "total": 0,
  "primary": 0,
  "secondary": 0,
  "supplementary": 0,
  "duplicates": 0,
  "primary duplicates": 0,
  "mapped": 0,
  "mapped %": null,
  "primary mapped": 0,
  "primary mapped %": null,
  "paired in sequencing": 0,
  "read1": 0,
  "read2": 0,
  "properly paired": 0,
  "properly paired %": null,
  "with itself and mate mapped": 0,
  "singletons": 0,
  "singletons %": null,
  "with mate mapped to a different chr": 0,
  "with mate mapped to a different chr (mapQ >= 5)": 0
 }
}
```

## Step 4: Visualize Alignments in IGV

1. Launch IGV.
2. Import your sorted and indexed BAM file (`aligned_reads_sorted.bam`).

What we've seen in IGV:

![IGV]()

## Step 5: Calculate Coverage for Each Chromosome

We solved this task in two different ways:

### Using SAMtools:

```sh
samtools depth -a aligned_reads_sorted.bam | awk '{cov[$1]+=$3; count[$1]++} END {for (chr in cov) print chr, cov[chr]/count[chr]}' > coverage_per_chromosome.txt
```

### Using BEDTools:

```sh
bedtools genomecov -ibam aligned_reads_sorted.bam -g reference_genome.fasta.fai > genome_coverage.txt
```

The results can be seen in the files `aufgabe02/samtools_cov_per_chr.txt` and `aufgabe02/bedtools_cov_per_chr.txt`.

## Day3: Variant Calling

The task was to call variants from the aligned reads (the `.bam` file).

### Step 1: Add Read Groups

Before calling variants, we had to add read groups (RG) infromation to the `.bam` file. We used the AddOrReplaceReadGroups from the GATK Picard tools:

```sh
/group/albi-praktikum2023/software/picard.jar AddOrReplaceReadGroups \
      I=/group/albi-praktikum2023/analysis/gruppe_3/gruppe3.bam \
      O=/group/albi-praktikum2023/analysis/gruppe_3/aufgabe03/pickard_output.bam \
      RGID=4 \
      RGLB=lib1 \
      RGPL=illumina \
      RGPU=unit1 \
      RGSM=20
``` 

### Step 2: Splitting bam file by chromosome

To harness the power of parallel computing, we split the `.bam` file into smaller files, each containing reads from a single chromosome. 
For this purpose we wrote following c skript, which was based on the followig command:

```sh
samtools view -b input.bam chrN > chrN.bam
```

```c
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Adjust these values based on your actual requirements
#define NUM_THREADS 24 // Example: Human chromosomes 1-22, X, and Y
char *chromosomes[] = {"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                       "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
                       "chr20", "chr21", "chr22", "chrX", "chrY"}; // List of chromosomes

typedef struct {
    int thread_id;
    char chromosome[10];
} thread_data;

void *splitBamFile(void *threadarg) {
    thread_data *my_data;
    my_data = (thread_data *) threadarg;
    int tid = my_data->thread_id;
    char *chr = my_data->chromosome;

    printf("Thread %d splitting BAM for %s\n", tid, chr);

    char command[256];
    sprintf(command, "/group/albi-praktikum2023/software/samtools view -b /group/albi-praktikum2023/analysis/gruppe_3/aufgabe03/pickard_output.bam %s > %s.bam", chr, chr);
    system(command);

    printf("Thread %d finished splitting BAM for %s\n", tid, chr);
    pthread_exit(NULL);
}

int main () {
    pthread_t threads[NUM_THREADS];
    thread_data thread_data_array[NUM_THREADS];
    int rc;
    long t;

    for(t = 0; t < NUM_THREADS; t++) {
        printf("In main: creating thread %ld for %s\n", t, chromosomes[t]);
        thread_data_array[t].thread_id = t;
        strcpy(thread_data_array[t].chromosome, chromosomes[t]);
        rc = pthread_create(&threads[t], NULL, splitBamFile, (void *)&thread_data_array[t]);
        if (rc) {
            printf("ERROR; return code from pthread_create() is %d\n", rc);
            exit(-1);
        }
    }

    // Wait for all threads to complete
    for(t = 0; t < NUM_THREADS; t++) {
        pthread_join(threads[t], NULL);
    }

    printf("Main: BAM splitting completed. Exiting.\n");
    pthread_exit(NULL);
}
```

### Step 3: Variant Calling

We used the GATK to perform Variant Calling. The following command were crusial for this task:

- `MarkDuplicates` to mark duplicates in the `.bam` file.
```sh
gatk MarkDuplicates \
    -I chrN.bam \
    -O chrN_marked_duplicates.bam \
    -M chrN_marked_dup_metrics.txt
```

- `HaplotypeCaller` to call variants.
```sh
gatk HaplotypeCaller \
    -R reference.fasta \
    -I chrN_marked_duplicates.bam \
    -O chrN_raw_variants.vcf
```

- `VariantFiltration` to filter the variants.
```sh
gatk VariantFiltration \
    -R reference.fasta \
    -V chrN_raw_variants.vcf \
    -O chrN_filtered_snps.vcf \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \
    --filter-name "SNP_FILTER"
```

Also, we had to index each `.bam` file, which was done using the following command:

```sh
samtools index %s_marked.bam
```

Similarly to previous step, we wrote a c skript to run these commands in parallel:

```c
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NUM_THREADS 24 // This should match the number of chromosomes or tasks

typedef struct {
    int thread_id;
    char chromosome[24]; // Adjust size as needed
} thread_data;

void *processChromosome(void *threadarg) {
    thread_data *my_data;
    my_data = (thread_data *) threadarg;
    int tid = my_data->thread_id;
    char *chr = my_data->chromosome;

    printf("Thread %d processing chromosome %s\n", tid, chr);

    // Construct command strings (simplified here)
    char markDupCmd[256];
    sprintf(markDupCmd, "/group/albi-praktikum2023/software/gatk MarkDuplicates -I %s.bam -O %s_marked.bam -M %s_metrics.txt", chr, chr, chr);

    char indexBamFileCmd[256];
    sprintf(indexBamFileCmd, "samtools index %s_marked.bam", chr);

    char callVarCmd[256];
    sprintf(callVarCmd, "/group/albi-praktikum2023/software/gatk HaplotypeCaller -R /group/albi-praktikum2023/data/NCBI-reference/GRCh38.fna -I %s_marked.bam -O %s_raw.vcf", chr, chr);

    char filterSNPsCmd[256];
    sprintf(filterSNPsCmd, "/group/albi-praktikum2023/software/gatk VariantFiltration -R /group/albi-praktikum2023/data/NCBI-reference/GRCh38.fna -V %s_raw.vcf -O %s_filtered.vcf --filter-expression \"QD < 2.0 || FS > 60.0 || MQ < 40.0\" --filter-name \"SNP_FILTER\"", chr, chr);

    // Execute commands
    system(markDupCmd);
    system(indexBamFileCmd);
    system(callVarCmd);
    system(filterSNPsCmd);

    printf("Thread %d finished processing chromosome %s\n", tid, chr);
    pthread_exit(NULL);
}

int main () {
    pthread_t threads[NUM_THREADS];
    thread_data thread_data_array[NUM_THREADS];
    int rc;
    long t;

    char *chromosomes[] = {"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"};

    for(t = 0; t < NUM_THREADS; t++) {
        printf("In main: creating thread %ld\n", t);
        thread_data_array[t].thread_id = t;
        strcpy(thread_data_array[t].chromosome, chromosomes[t]);
        rc = pthread_create(&threads[t], NULL, processChromosome, (void *)&thread_data_array[t]);
        if (rc) {
            printf("ERROR; return code from pthread_create() is %d\n", rc);
            exit(-1);
        }
    }

    // Wait for all threads to complete
    for(t = 0; t < NUM_THREADS; t++) {
        pthread_join(threads[t], NULL);
    }

    printf("Main: program completed. Exiting.\n");
    pthread_exit(NULL);
}
```

### Step 4: Merge VCFs

After calling variants for each chromosome, we merged the resulting VCF files into a single VCF file using the `MergeVcfs` tool from the GATK:

```sh
gatk MergeVcfs \
    -I chr1_filtered_snps.vcf \
    -I chr2_filtered_snps.vcf \
    ... \
    -O merged_filtered_snps.vcf
```

Finally, we received the merged VCF file `merged_filtered_snps.vcf`, which was used for further analysis.

### Step 5: Analyse and filtering of given VCF file
We had to perform in 5 steps on this stage.

1. Checking what our file contains:
- Chromosome, f.e. chr1
- Position, f.e. 21896830
- ID, f.e. rs140261216
- Reference allele, f.e. CTG 
- Alternate allele(s) f.e. CACAC,TATACACAC,TCACACA,TCACACACA,TCACACACACA  
- Quality score (QUAL) f.e. . 
- Filter status (FILTER) f.e. .
- Info fields (INFO) f.e. RS=139323070;RSPOS=6831944;dbSNPBuildID=134;SSR=0;SAO=0;VP=0x050000080005040036000100;GENEINFO=CAMTA1:23261;WGT=1;VC=SNV;INT;ASP;VLD;KGPhase1;KGPhase3;CAF=0.999,0.0009984;COMMON=1;TOPMED=0.99681447502548419,0.00318552497451580

2. Interpret INFO field
- RS (Reference SNP ID) - a unique identifier in the dbSNP database
- RSPOS (Reference SNP Position) - the genomic position of the variant
- RV (Reference allele(s) and Variation allele(s)) - this field is unspecified and can vary in this context. It could specify the reference allele and the variant allele or other relevant information
- dbSNPBuildID - the version number of the dbSNP database in which this variant is included
- SSR (SubSNP rsID) - an internal reference number for the variant in the dbSNP database
- SAO (Variant Allele Origin) - a value indicating whether the variant allele was observed in the reference sequence (SAO=0) or not (SAO=1)
- VP (Variant properties) - a bitflag encoding various information about the variant
- GENEINFO - information about the gene affected by the variant 
- WGT (Weight) - a weighting factor for the variant, usually 1
- VC (Variant Class) - the class of the variant
- ASP (Allele-specific properties) - information about the allele of the variant
- VLD (Variant Locus-specific database) - information about whether the variant is present in a specific locus-specific database.
- G5 (Allele frequency group) - information about the frequency of the variant allele in different population groups.
- KGPhase1 and KGPhase3 - information about the frequency of the variant allele in the 1000 Genome Phase 1 and Phase 3 databases, respectively.
- CAF (Combined Allele Frequency) - the combined allele frequency of the variant allele in different populations.
- COMMON - a flag indicating whether the variant is common in the population (COMMON=1) or not (COMMON=0).
- TOPMED - the allele frequency of the variant allele in the TOPMed database.

3. Testing if the file contains only SNP (SNP)
```sh
zcat /group/albi-praktikum2023/data/dbSNP/dbSNP.vcf.gz | grep -v '^#' | awk '{if ($8 !~ /VT=SNP/ && $8 !~ /VC=SNV/) print}' | cut -f8
```
zcat - command to decompress the file with .gz format on-the-fly
grep -v '^#' - filters out header lines (lines starting with #)
awk - program used to process and format text lines
if ($8 !~ /VT=SNP/ && $8 !~ /VC=SNV/) - if-condition to find variants types and classes, which are not SNP or SNV
print - if statemenet is true, makes output 
cut -f8 - extracts the eighth field (INFO) from each line of input

In output, f.e, were:
RS=200339231;RSPOS=4325724;dbSNPBuildID=137;SSR=0;SAO=0;VP=0x05000000000504003e000200;WGT=1;VC=DIV;ASP;VLD;KGPhase1;KGPhase3;CAF=0.9846,0.003594;COMMON=1;TOPMED=0.99357320336391437,0.00642679663608562

So VC=DIV, it means file contains not only SNP.

4. FIltering file dbSNP.vcf.gz and creating new file dbSNPfiltered.vcf, which contains only SNP(SNV)
```sh
zcat /group/albi-praktikum2023/data/dbSNP/dbSNP.vcf.gz | grep -v '^#' | awk -F '\t' '$8 ~ /VT=SNP/ || $8 ~ /VC=SNV/' > /group/albi-praktikum2023/analysis/gruppe_3/aufgabe03/dbSNPfiltered.vcf
```
-F '\t' - this option sets the field separator to a tab character (\t). This is used to split each line into fields based on tabs
'$8 ~ /VT=SNP/ || $8 ~ /VC=SNV/' - This condition checks if the eighth field ($8) contains the string "VT=SNP" OR "VC=SNV"

5. Testing new created file, if it contains only SNP(SNV)
```sh
cat /group/albi-praktikum2023/analysis/gruppe_3/aufgabe03/dbSNPfiltered.vcf | grep -v '^#' | awk '{if ($8 !~ /VT=SNP/ && $8 !~ /VC=SNV/) print}' | cut -f8
```
Ouput contained no strings, which means new file consist only of SNP

### Step 6 Determination of the Individual's Sex Based on X Chromosome SNP Analysis
We tried to determine the sex of an individual by analyzing the ratio of heterozygous single nucleotide polymorphisms (SNPs) to total SNPs on the X chromosome.

1. Setting the Path to the VCF File
First, specify the path to the VCF (Variant Call Format) file containing SNP data for the X chromosome.
```bash
vcf_path="/group/albi-praktikum2023/analysis/gruppe_3/aufgabe03/chrX_filtered.vcf"
```

2 The following command counts the number of heterozygous SNPs present on the X chromosome.
```bash
heterozygous_snps=$(/group/albi-praktikum2023/software/bcftools view -H $vcf_path | awk '{if($10 ~ "0/1" || $10 ~ "1/0") print $0}' | wc -l)
```

3. Then we count the total number of SNPs present on the X chromosome.
```bash
total_snps=$(/group/albi-praktikum2023/software/bcftools view -H $vcf_path | wc -l)
```

4. The ratio is calculated as follows:
```bash
ratio=$(echo "scale=2; $heterozygous_snps / $total_snps" | bc)
```

5. Analysis Result:
Based on the results of the analysis, with 87,296 heterozygous SNPs and a total of 156,862 SNPs on the X chromosome, the ratio is 0.55. This relatively high ratio suggests that the individual is likely female.

## Additional Information

N/A

## Supervisors

- Svenja Mehringer

## Contributors

- Dzmitry Hramyka
- Sofya Shorzhina
- Sofiia Ilnytska
- Maksim Danilchyk
