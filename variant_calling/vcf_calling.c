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
    // char markDupCmd[256];
    // sprintf(markDupCmd, "/group/albi-praktikum2023/software/gatk MarkDuplicates -I %s.bam -O %s_marked.bam -M %s_metrics.txt", chr, chr, chr);

    char indexBamFileCmd[256];
    sprintf(indexBamFileCmd, "samtools index %s_marked.bam", chr);

    char callVarCmd[256];
    sprintf(callVarCmd, "/group/albi-praktikum2023/software/gatk HaplotypeCaller -R /group/albi-praktikum2023/data/NCBI-reference/GRCh38.fna -I %s_marked.bam -O %s_raw.vcf", chr, chr);

    char filterSNPsCmd[256];
    sprintf(filterSNPsCmd, "/group/albi-praktikum2023/software/gatk VariantFiltration -R /group/albi-praktikum2023/data/NCBI-reference/GRCh38.fna -V %s_raw.vcf -O %s_filtered.vcf --filter-expression \"QD < 2.0 || FS > 60.0 || MQ < 40.0\" --filter-name \"SNP_FILTER\"", chr, chr);

    // Execute commands
    // system(markDupCmd);
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

