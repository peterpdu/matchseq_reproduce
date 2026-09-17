import torch
import numpy as np
import pandas as pd

class StandardNMF:
    """
    Non-negative matrix factorization using multiplicative updates
    """
    def __init__(self, n_components, max_iter=500, tol=1e-4, device='cuda' if torch.cuda.is_available() else 'cpu'):
        self.k = n_components
        self.max_iter = max_iter
        self.tol = tol
        self.device = device
        self.W = None
        self.H = None

    def fit(self, X):
        """
        returns W, H, and final reconstruction error.
        """
        if not isinstance(X, torch.Tensor):
            X = torch.tensor(X, dtype=torch.float32).to(self.device)
            
        if torch.min(X) < 0:
            raise ValueError("data contains negatives")

        m, n = X.shape
        
        # init with uniform random
        self.W = torch.rand(m, self.k, device=self.device, dtype=X.dtype)
        self.H = torch.rand(self.k, n, device=self.device, dtype=X.dtype)
        
        epsilon = 1e-10 
        prev_error = float('inf')
        final_error = float('inf')
        
        for i in range(self.max_iter):
            # W <- W * (X H^T) / (W H H^T)
            Ht = self.H.T
            numerator_W = X @ Ht
            denominator_W = self.W @ (self.H @ Ht) + epsilon
            self.W = self.W * (numerator_W / denominator_W)
            
            # H <- H * (W^T X) / (W^T W H)
            Wt = self.W.T
            numerator_H = Wt @ X
            denominator_H = (Wt @ self.W) @ self.H + epsilon
            self.H = self.H * (numerator_H / denominator_H)

            # convergence check
            if i % 10 == 0 or i == self.max_iter - 1:
                recon = self.W @ self.H
                final_error = torch.norm(X - recon, p='fro').item()
                
                if abs(prev_error - final_error) < self.tol:
                    break
                prev_error = final_error

        return self.W, self.H, final_error
    
    def transform(self, X):
        """
        Computes the scores (W) for new data X using the learned H.
        
        X: (n_new_samples, n_features)
        returns: W_new (n_new_samples, n_components)
        """
        if self.H is None:
            raise RuntimeError("Model must be fit before transform.")
            
        if not isinstance(X, torch.Tensor):
            X = torch.tensor(X, dtype=torch.float32).to(self.device)
            
        m, n = X.shape
        if n != self.H.shape[1]:
            raise ValueError(f"H has {self.H.shape[1]} features, input has {n}")

        W_new = torch.rand(m, self.k, device=self.device, dtype=X.dtype)
        epsilon = 1e-10
        
        Ht = self.H.T
        HHt = self.H @ Ht
        
        for i in tqdm(range(self.max_iter)): 
            numerator = X @ Ht
            denominator = W_new @ HHt + epsilon
            W_new = W_new * (numerator / denominator)
            
        return W_new

