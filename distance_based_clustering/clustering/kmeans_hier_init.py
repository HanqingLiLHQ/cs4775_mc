import random
import numpy as np

class TreeNode(object):
    """
    Data structure for the hierarchical tree. 
    This object is either an inner node or a leaf. 
    If it is an inner node:
    both left and right attributes should be another treenode object.
    height > 0. name should be "inner#". mn (member number)
    If it is a leaf:
    Both left and right attribute should be None.
    Height = 0. name should be a specific motif id. 
    mn = 1
    """
    def __init__(self, name, left_child, right_child, height, mn):
        self.left = left_child
        self.right = right_child
        self.name = name
        self.mn = mn
        self.height = height
        # if this is an innernode: 
        if self.left != None: 
            assert self.height > 0
            assert self.left.height <= self.height
            assert self.right.height <= self.height


# This algorithm will first build a tree represented by class "TreeNode", then traverse the treenode
# for heights from highest to the lowest until it finds k-1 nodes. After neglecting these top n-1 nodes, 
# the rest tree nodes just represents how to group the motifs into clusters. 

def wrap_dm(idlist, dm):
    """
    helper functions to first wrap the motif idlist into a leaf list, and to convert the distance matrix into a distance map. 
    Precondition: 
        idlist: a list of motif ids. 
        dm: a np array representing the distance matrix. it should have a shape of len(idlist) * len(idlist)
    Returns:
        nodemap: a map of treenodes: treenode.name -> treenode
        distance_map: a map of treenode.name -> treenode.name -> distance
    """
    nodemap = {id : TreeNode(id, None, None, 0, 1) for id in idlist}
    distance_map = {}
    for i in range(len(idlist)):
        motif_id1 = idlist[i]
        id1_map = {}
        for j in range(len(idlist)):
            id1_map[idlist[j]] = dm[i][j]
        distance_map[motif_id1] = id1_map
    return nodemap, distance_map

def shortest(dm):
    """
    helper function that finds the shortest distance between two different nodes in a distance map.
    Returns the shortest distance as well as the two TreeNode.names
    Precondition:
        dm: a map of treenode.name -> treenode.name -> distance
    """
    ret = (None, None, np.inf)
    for key1 in dm.keys():
        for key2 in dm[key1].keys():
            if key1 != key2 and dm[key1][key2] < ret[2]:
                ret = [key1, key2, dm[key1][key2]]
    return ret

def tree_construction(dm, idlist, dc = "complete"):
    """
    Returns a TreeNode object denoting the current tree
    Precondition: 
        idlist: a list of motif ids. 
        dm: a np array representing the distance matrix. it should have a shape of len(idlist) * len(idlist)
        dc: the way for distance update method. Should be one of "complete", "single", "UPGMA", "WPGMA:
    """
    nodemap, dm = wrap_dm(idlist,dm)
    # need to build treenode for len(list(nodemap.keys()))-1 times
    for iter in range(len(list(nodemap.keys()))-1):
        m1, m2, sd = shortest(dm)
        
        # create the new node and update the nodemap
        new_node_name = iter
        left_node = nodemap[m1]
        right_node = nodemap[m2]
        lmn = left_node.mn
        rmn = right_node.mn
        mn = lmn + rmn
        height = sd / 2
        new_node = TreeNode(name = new_node_name, left_child = left_node, right_child = right_node, height = height, mn = mn)
        nodemap[new_node_name] = new_node 
        del nodemap[m1]
        del nodemap[m2]

        # update the distance matrix
        del dm[m1]
        del dm[m2]
        new_node_distances = {new_node_name : 0}
        for m in dm.keys():
            distance_l = dm[m][m1]
            distance_r = dm[m][m2]
            del dm[m][m1]
            del dm[m][m2]
            if dc == "single":
                new_distance = min(distance_l, distance_r)
            elif dc == "complete":
                new_distance = max(distance_l, distance_r)
            elif dc == "UPGMA":
                new_distance = (distance_l * lmn + distance_r * rmn ) / mn
            elif dc == "WPGMA":
                new_distance = (distance_l + distance_r) / 2
            dm[m][new_node_name] = new_distance
            new_node_distances[m] = new_distance
        dm[new_node_name] = new_node_distances
    assert len(list(nodemap.keys())) == 1
    return list(nodemap.values())[0]
            
        

def cut_tree(tree, nclusters):
    """
    Cut the given tree object to give nclusters
    Precondition: 
        tree is a TreeNode object
        nclusters is an integer >= 1
    """
    height_map = {tree.height : [tree]}
    remaining_steps = nclusters - 1
    while remaining_steps > 0:
        # find the highest node and split it into two child tree nodes
        # there might still be a need to handle the really rare case that actually two nodes are at the same height
        highest = max(list(height_map.keys()))
        highest_node = height_map[highest][0]
        if len(height_map[highest]) > 1:
            height_map[highest].pop(0)
        else:
            del height_map[highest]
        left_node = highest_node.left
        right_node = highest_node.right
        if left_node.height not in height_map: 
            height_map[left_node.height] = [left_node]
        else:
            height_map[left_node.height].append(left_node)
        if right_node.height not in height_map: 
            height_map[right_node.height] = [right_node]
        else:
            height_map[right_node.height].append(right_node)
        remaining_steps -= 1
    node_clusters = []
    for height in height_map.keys():
        node_clusters += height_map[height]
    assert len(node_clusters) == nclusters
    # now it is just to get the motif names from the TreeNode objects
    def degradenode(treenode):
        """
        helper function to degrade the nodes from the treeNode objects
        """
        if treenode.left == None:
            return [treenode.name]
        return degradenode(treenode.left) + degradenode(treenode.right)
    clusters = []
    for node in node_clusters:
        clusters.append(degradenode(node))
    return clusters
    

def hierarchical_clustering(nclusters, dm, idlist, distance_update = "UPGMA"):
    """
    Grow a tree, then cut it!
    Preconditions:
        nclusters: an integer >= 1
        dm: np array for the distance matrix
        idlist: the motif ids
        distance_update: one of ["complete", "single", "UPGMA", "WPGMA"]
    """
    return cut_tree(tree_construction(dm, idlist, dc = distance_update), nclusters)


def select_initial_centroids(num_clusters, dm):
    # Randomly choose 'num_clusters' different indices to serve as initial centroids
    indices = list(range(len(dm)))
    random.shuffle(indices)
    centroids_indices = indices[:num_clusters]
    return centroids_indices

def assign_clusters(centroids_indices, dm):
    clusters = [[] for _ in range(len(centroids_indices))]
    for i, distances_to_centroids in enumerate(dm):
        # Find the closest centroid for each motif
        closest_centroid = np.argmin([dm[i][c] for c in centroids_indices])
        clusters[closest_centroid].append(i)
    return clusters

def find_new_centroids(clusters, dm):
    new_centroids_indices = []
    for cluster in clusters:
        # Find the motif in the cluster that has the minimum average 
        # distance to all other motifs in the cluster
        if cluster:
            min_avg_distance = float('inf')
            centroid_index = None
            for motif_index in cluster:
                avg_distance = np.mean([dm[motif_index][other_index] for other_index in cluster])
                if avg_distance < min_avg_distance:
                    min_avg_distance = avg_distance
                    centroid_index = motif_index
            new_centroids_indices.append(centroid_index)
    return new_centroids_indices

def kmeans_motifs(dm, id_list, num_clusters, max_iter=100):
    # Initialize centroids
    init_clusters_with_ids = hierarchical_clustering(num_clusters, dm, id_list)
    init_clusters = [[id_list.index(id) for id in ids] for ids in init_clusters_with_ids]
    clusters = init_clusters
    centroids_indices = []
    
    for iteration in range(max_iter):
        # Find new centroids within the current clusters
        new_centroids_indices = find_new_centroids(clusters, dm)
        
        # Assign clusters based on the current centroids
        clusters = assign_clusters(new_centroids_indices, dm)
        
        # Check for convergence (if centroids do not change)
        if set(new_centroids_indices) == set(centroids_indices):
            break
        
        centroids_indices = new_centroids_indices

    # Transform the indices to the actual motif ids
    clusters_with_ids = [[id_list[index] for index in cluster] for cluster in clusters]
    return clusters_with_ids