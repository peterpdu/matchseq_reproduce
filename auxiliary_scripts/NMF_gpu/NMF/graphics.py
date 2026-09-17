import numpy as np
from scipy.spatial.distance import pdist, squareform

import matplotlib.pyplot as plt

def plot_all_clustered_components(pool_of_components, labels, k):
    """
    Plots all components from all runs, grouped by their K-Means cluster assignment.
    Creates k subplots.
    """
    fig, axes = plt.subplots(k, 1, figsize=(10, 3*k), sharex=True)
    if k == 1: axes = [axes]
    
    colors = plt.cm.viridis(np.linspace(0, 1, k))
    
    print(f"\nVisualizing component stability for Rank k={k}...")
    
    for i in range(k):
        ax = axes[i]
        cluster_indices = np.where(labels == i)[0]
        
        if len(cluster_indices) == 0:
            ax.set_title(f"Cluster {i} (Empty)")
            continue
            
        cluster_data = pool_of_components[cluster_indices]
        
        ax.plot(cluster_data.T, color=colors[i], alpha=0.2, linewidth=1)
        
        cluster_mean = np.mean(cluster_data, axis=0)
        ax.plot(cluster_mean, color='black', linestyle='--', linewidth=2, label='Cluster Mean')
        
        ax.set_title(f"Component Cluster {i+1} (Contains {len(cluster_indices)} vectors from all runs)")
        ax.set_ylabel("Value")
        ax.legend(loc='upper right')
        
    axes[-1].set_xlabel("Feature Index")
    plt.suptitle(f"Stability Analysis: All Runs Clustered (Rank {k})", fontsize=14)
    plt.tight_layout(rect=[0, 0.03, 1, 0.97])
    plt.show()


def plot_sweep_metrics(results, left_key='avg_error', right_key='silhouette'):
    ranks = sorted(results.keys())
    left_values = [results[k][left_key] for k in ranks]
    right_values = [results[k][right_key] for k in ranks]
    
    fig, ax1 = plt.subplots(figsize=(10, 6))
    color = 'tab:red'
    ax1.set_xlabel('Rank (k)')
    ax1.set_ylabel(left_key, color=color)
    ax1.plot(ranks, left_values, color=color, marker='o', label=left_key)
    
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.grid(True, alpha=0.3)

    ax2 = ax1.twinx() 
    color = 'tab:blue'
    ax2.set_ylabel(right_key, color=color)
    ax2.plot(ranks, right_values, color=color, marker='s', linestyle='--', label=right_key)
    
    ax2.tick_params(axis='y', labelcolor=color)
    plt.title('Rank Selection')
    fig.tight_layout() 
    plt.show()
