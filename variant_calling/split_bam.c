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
