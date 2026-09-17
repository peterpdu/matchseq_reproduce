import torch
import numpy as np

from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist, squareform
from scipy.cluster.hierarchy import linkage, cophenet

def compute_consensus_matrix(W_array):
    """Computes the consensus matrix from W array"""
    n_iter, n_samples, k = W_array.shape
    consensus = np.zeros((n_samples, n_samples))
    for W in W_array:
        labels = np.argmax(W, axis=1)
        connectivity = (labels[:, None] == labels[None, :]).astype(float)
        consensus += connectivity
    return consensus / n_iter


def compute_cophenetic_correlation(consensus_matrix):
    """Computes the Cophenetic Correlation Coefficient."""
    dist_matrix = 1.0 - consensus_matrix
    np.fill_diagonal(dist_matrix, 0.0)
    condensed_dist = squareform((dist_matrix + dist_matrix.T) / 2.0, checks=False)
    Z = linkage(condensed_dist, method='average')
    coph_corr, _ = cophenet(Z, condensed_dist)
    return coph_corr


def calculate_stability(pool_of_components, k):
    """
    Clusters the aggregated components and calculates Silhouette Coefficient.
    
    Args:
        pool_of_components: (n_runs * k, n_features) array of all W columns
        k: The rank (number of clusters expected)
        
    Returns:
        silhouette: float, stability score (-1 to 1)
        centroids: (k, n_features) cluster centers
        labels: cluster assignments
    """
    if k < 2:
        return 0.0, np.mean(pool_of_components, axis=0, keepdims=True), np.zeros(len(pool_of_components))

    kmeans = KMeans(n_clusters=k, init='k-means++', max_iter=1000, tol=1e-4, random_state=23)
    labels = kmeans.fit_predict(pool_of_components)
    centroids = kmeans.cluster_centers_
    
    try:
        sil_score = silhouette_score(pool_of_components, labels)
    except ValueError:
        print('All points are in one cluster')
        sil_score = 0.0
        
    return sil_score, centroids, labels


def select_best_run_from_W(run_storage, centroids):
    """
    Finds the single run that is closest to the K-Means centroids.
    """
    best_score = float('inf')
    best_idx = -1
    
    for idx, run in enumerate(run_storage):
        W_run = run['W'].T # Shape (k, m)
        
        dists = cdist(W_run, centroids, metric='euclidean')
        row_ind, col_ind = linear_sum_assignment(dists)
        
        total_dist = dists[row_ind, col_ind].sum()
        
        if total_dist < best_score:
            best_score = total_dist
            best_idx = idx
            
    return best_idx


def select_best_run_from_H(run_storage, centroids):
    """
    Selects the single run whose H matrix is closest to the consensus centroids.
    """
    best_score = float('inf')
    best_idx = -1
    
    for idx, run in enumerate(run_storage):
        H_run = run['H'] # Shape (k, n_features)
        
        # Calculate distance between this run's components and consensus centroids
        dists = cdist(H_run, centroids, metric='euclidean')
        
        # Optimal assignment to match components to centroids
        row_ind, col_ind = linear_sum_assignment(dists)
        
        total_dist = dists[row_ind, col_ind].sum()
        
        if total_dist < best_score:
            best_score = total_dist
            best_idx = idx
            
    return best_idx

