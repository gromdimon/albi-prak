import numpy as np
import matplotlib.pyplot as plt
from sklearn.manifold import MDS

# Load the distance matrix
distance_matrix = np.loadtxt('hamming_distance_matrix.txt')

# Load sample names
with open('sample_names.txt', 'r') as f:
    samples = [line.strip() for line in f.readlines()]


# Initialize MDS model
mds = MDS(n_components=2, dissimilarity='precomputed', random_state=42)

# Fit the model & transform data
mds_coords = mds.fit_transform(distance_matrix)  # Shape: (n_samples, 2)

# Plot the results
plt.figure(figsize=(10, 8))
plt.scatter(mds_coords[:, 0], mds_coords[:, 1])

# Annotate points with sample names
for i, sample in enumerate(samples):
    plt.text(mds_coords[i, 0], mds_coords[i, 1], sample)

plt.title('2D MDS Projection of Hamming Distances')
plt.xlabel('MDS Dimension 1')
plt.ylabel('MDS Dimension 2')
plt.grid(True)
plt.savefig('distance_matrix_mds.png')
plt.show()
