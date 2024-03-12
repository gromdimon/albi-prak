# Algorithmische Bioinformatik: Praktikum

The practical part of the course "Algorithmische Bioinformatik" (Algorithmic Bioinformatics) at the Freie Universität Berlin.

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

Since the FastQC was already installed on the server, we could run it using the following command:

```sh
fastqc /group/albi-praktikum2023/analysis/gruppe_3/aufgabe01/test_gruppe3_1.fastq.gz
```

...

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
samtools flagstat aligned_reads_sorted.bam
```

This was the output of the command:

```sh
100 + 0 in total (QC-passed reads + QC-failed reads)
100 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
100 + 0 mapped (100.00% : N/A)
100 + 0 primary mapped (100.00% : N/A)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (N/A : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
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

The outputs were correspondingly:
    
```sh
pass
```

## Additional Information

N/A

## Supervisors

- Svenja Mehringer

## Contributors

- Dzmitry Hramyka
- Sofia Schorzina
- Sofia Ilyina
- Maxim Danilchik