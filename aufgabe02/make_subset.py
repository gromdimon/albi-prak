import gzip
from Bio import SeqIO

def subset_fastq_gz(input_file, output_file, num_records):
    """
    Create a subset of a .fastq.gz file.

    Parameters:
    input_file (str): Path to the input .fastq.gz file.
    output_file (str): Path where the subset .fastq.gz file will be saved.
    num_records (int): Number of records to include in the subset.
    """
    print("Open files")
    with gzip.open(input_file, "rt") as in_handle:  # Open input file in read text mode
        with gzip.open(output_file, "wt") as out_handle:  # Open output file in write text mode
            records = (record for i, record in enumerate(SeqIO.parse(in_handle, "fastq")) if i < num_records)
            count = SeqIO.write(records, out_handle, "fastq")
    print(f"Done. Subset created with {count} records.")

# Example usage
input_fastq_gz = "/group/albi-praktikum2023/data/1000-Genome-Project/gruppe3/gruppe3_1.fastq.gz"
output_fastq_gz = "/group/albi-praktikum2023/analysis/gruppe_3/aufgabe02/subset_output_file.fastq.gz"
subset_fastq_gz(input_fastq_gz, output_fastq_gz, 1000)
