import random
from sklearn.cluster import KMeans
import numpy as np

### This doesn't work because KMeans in numpy does not allow self-defined distances
# Just here for future referral, delete later

def kmeans_clustering(dataset, dm, id_list):
    motifs_map = dataset.mmm
    motif_ids = list(motifs_map.keys())
    data = []

    for motif_id in motif_ids:
        motif = motifs_map[motif_id]
        pssm = motif.pssm
        print("Matrix here\n")
        print(pssm)
        print("Checkpoint A\n")
        print(pssm[0])
        print(type(pssm[0]))
        # Here you should implement a method to convert your PSSM into a fixed-size vector
        # CUstomize input vector into K-means for each pssm
        
        scores = []
        for i in range(len(pssm)):  # pos_dict is a dictionary for each position in the motif
            scores = scores + pssm[i]
        
        data.append(scores)
        print("Checkpoint B\n")
        print(len(scores))

    print("Checkpoint C\n")
    # print(data)
    print(type(data))
    # Convert list of vectors to a 2D NumPy array
    data_array = np.array(data)

    # K-means clustering
    num_clusters = 3  # for example, set the desired number of clusters
    kmeans = KMeans(n_clusters=num_clusters, random_state=0).fit(data_array)

    # Output the clusters
    clusters = {i: [] for i in range(num_clusters)}
    for i, label in enumerate(kmeans.labels_):
        clusters[label].append(motif_ids[i])
    
    return clusters