import torch
import numpy as np

from tqdm import tqdm

from .utils import calculate_stability, select_best_run_from_W, select_best_run_from_H
from .utils import compute_consensus_matrix, compute_cophenetic_correlation
from .nmf import StandardNMF

def run_rank_sweep(X_data, rank_range, n_runs=10, device='cpu', model_cls=StandardNMF):
    """
    Performs parameter sweep using the specified model class (SemiNMF or StandardNMF).
    """
    X_tensor = torch.tensor(X_data, dtype=torch.float32).to(device)
    
    final_results = {}
    
    for k in rank_range:
        print(f"--- Sweeping Rank k={k} using {model_cls.__name__} ---")
        storage = []
        errors = []
        
        # perform n runs
        for r in tqdm(range(n_runs)):
            model = model_cls(n_components=k, max_iter=200, device=device)
            W, H, error = model.fit(X_tensor)
            
            storage.append({
                'W': W.detach().cpu().numpy(),
                'H': H.detach().cpu().numpy(),
                'error': error
            })
            errors.append(error)

        # aggregate all W components for stability analysis
        all_vectors = [run['W'].T for run in storage]
        W_pool = np.vstack(all_vectors) # (n_runs * k, features)
        all_vectors = [run['H'].T for run in storage]
        H_array = np.stack(all_vectors, axis=0) # (n_runs, features, k)
        
        sil_score, centroids, labels = calculate_stability(W_pool, k)
        consensus = compute_consensus_matrix(H_array)
        ccc = compute_cophenetic_correlation(consensus)

        # pick run closest to centroid
        best_idx = select_best_run_from_W(storage, centroids)
        best_model = storage[best_idx]
        
        avg_error = np.mean(errors)
        std_error = np.std(errors)
        
        print(f"  > Avg Error: {avg_error:.4f} | Silhouette (W)={sil_score:.4f} | Cophenetic (H)={ccc:.4f}")
        
        final_results[k] = {
            'W': best_model['W'],
            'H': best_model['H'],
            'centroids': centroids,
            'all_components': W_pool,
            'labels': labels,
            'run_index': best_idx,
            'avg_error': avg_error,
            'std_error': std_error,
            'silhouette': sil_score,
            'H_consensus': consensus,
            'H_ccc': ccc,
            'H_array': H_array,
        }
        
    return final_results

