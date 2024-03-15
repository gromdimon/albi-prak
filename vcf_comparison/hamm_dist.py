import allel
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load VCF
vcf_path = 'merged_samples_no_header.vcf'
callset = allel.read_vcf(vcf_path, fields=['samples', 'calldata/GT'])

# Extract genotypes and sample names
gt_array = allel.GenotypeArray(callset['calldata/GT'])
samples = callset['samples']

# Initialize an empty matrix for Hamming distances
num_samples = len(samples)
hamming_matrix = np.zeros((num_samples, num_samples), dtype=float)

# Compute Hamming distance
for i in range(num_samples):
    for j in range(num_samples):
        if i != j:
            # Flatten genotype arrays for comparison
            gt_i = gt_array[:, i].to_n_alt(fill=-1)
            gt_j = gt_array[:, j].to_n_alt(fill=-1)
            
            # Compute Hamming distance
            diff = np.sum(gt_i != gt_j)
            hamming_matrix[i, j] = diff / gt_array.shape[0]
        else:
            hamming_matrix[i, j] = 0

# Round the numbers in the matrix to 5 digits
hamming_matrix = np.round(hamming_matrix, 5)

# Save the samples and Hamming distance matrix
np.savetxt('samples.txt', samples, fmt='%s')
np.savetxt('hamming_distance_matrix.txt', hamming_matrix, fmt='%0.5f')

# Create a DataFrame for better labeling in the heatmap
df_hamming = pd.DataFrame(hamming_matrix, index=samples, columns=samples)

# Plot the heatmap
plt.figure(figsize=(20, 20))
sns.heatmap(df_hamming, annot=True, fmt="0.5f", cmap='coolwarm', linewidths=.5)
plt.title('Hamming Distance Matrix Heatmap')
plt.tight_layout()

# Save the heatmap as an image
plt.savefig('hamming_distance_heatmap.png')
plt.close()
