import allel
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Load the multi-sample VCF
vcf_path = 'merged_samples_no_header.vcf'
callset = allel.read_vcf(vcf_path)

# Extract genotype data
gt = allel.GenotypeArray(callset['calldata/GT'])

# Initialize a matrix to hold the Hamming distances
num_samples = gt.shape[1]
distance_matrix = np.zeros((num_samples, num_samples), dtype=int)

# Calculate the Hamming distances
for i in range(num_samples):
    for j in range(i + 1, num_samples):
        # Calculate the Hamming distance between samples i and j
        dist = np.sum(gt[:, i, :] != gt[:, j, :])
        distance_matrix[i, j] = dist
        distance_matrix[j, i] = dist  # The matrix is symmetric

# Sample names
samples = callset['samples']

# Save sample names to a text file
with open('sample_names.txt', 'w') as f:
    for sample in samples:
        f.write("%s\n" % sample)

# Plot the heat map
sns.heatmap(distance_matrix, annot=True, fmt="d", xticklabels=samples, yticklabels=samples)
plt.title('Hamming Distance Matrix')
plt.savefig('hamming_distance_heatmap.png')

