import numpy as np
import scipy.sparse as sp
import pandas as pd

def rwr_diffusion(adjacency_matrix: np.ndarray, 
                                  spreading_coeff: float = 0.2, 
                                  n_iter: int = 100, 
                                  fc_vector: pd.Series = None) -> pd.Series:
    """
    Vectorized RWR diffusion over a graph.
    """
    if fc_vector is None:
        raise ValueError("fc_vector must be provided.")

    # Degree normalization: compute D^{-1}A
    D = np.sum(adjacency_matrix, axis=1)
    D_inv = np.diag(1.0 / np.where(D == 0, 1, D))  # avoid division by zero
    W = D_inv @ adjacency_matrix

    # Optionally use sparse matrices for speed
    W = sp.csr_matrix(W)

    # Convert fold-change vector to numpy
    p0 = fc_vector.values.astype(float)
    p = p0.copy()

    for _ in range(n_iter):
        p_new = ((1 - spreading_coeff) * p0) + (spreading_coeff * (W @ p))
        # Convergence check (optional)
        if np.allclose(p_new, p, atol=1e-9):
            print("Converged.")
            break
        p = p_new

    print("RWR diffusion iterations completed.")
    return pd.Series(p, index=fc_vector.index)


import numpy as np
import pandas as pd
from numpy.linalg import norm

def rwr_validate(A, fc_vector, spreading_coeff, n_iter, RWR_fn,
                 k=50, B=30, noise_std=0.05):
    """
    A: adjacency matrix
    fc_vector: original seed fold-change vector
    spreading_coeff: alpha
    n_iter: iterations for RWR
    RWR_fn: function handle to your RWR implementation
    k: top-k for stability
    B: bootstraps
    noise_std: noise added to fc_vector to test robustness
    """

    # 1. --- RUN MAIN RWR ---
    pred = RWR_fn(adjacency_matrix=A,
                  fc_vector=fc_vector,
                  spreading_coeff=spreading_coeff,
                  n_iter=n_iter)

    pred = np.array(pred).flatten()
    p = pred / pred.sum()   # normalize

    # 2. --- CONVERGENCE ---  
    # Re-run with one less iteration to estimate Δ
    pred_prev = RWR_fn(adjacency_matrix=A,
                       fc_vector=fc_vector,
                       spreading_coeff=spreading_coeff,
                       n_iter=max(n_iter - 1, 1))
    pred_prev = np.array(pred_prev).flatten()
    convergence = norm(pred - pred_prev, 2)

    # 3. --- CONTRAST ---  
    idx_sorted = np.argsort(pred)[::-1]
    top = pred[idx_sorted[:k]].mean()
    rest = pred[idx_sorted[k:]].mean()
    contrast = top / (rest + 1e-12)

    # 4. --- STABILITY UNDER BOOTSTRAP ---
    topk_sets = []
    for b in range(B):
        # Add slight noise to input FCs
        noisy_fc = fc_vector + np.random.normal(0, noise_std, size=fc_vector.shape)

        pred_b = RWR_fn(adjacency_matrix=A,
                        fc_vector=noisy_fc,
                        spreading_coeff=spreading_coeff,
                        n_iter=n_iter)

        pred_b = np.array(pred_b).flatten()
        idx = np.argsort(pred_b)[::-1][:k]
        topk_sets.append(set(idx))

    # pairwise Jaccard stability
    overlaps = []
    for i in range(B):
        for j in range(i + 1, B):
            s = len(topk_sets[i] & topk_sets[j]) / len(topk_sets[i] | topk_sets[j])
            overlaps.append(s)

    stability = np.mean(overlaps)

    # return all metrics
    return {
        "convergence": convergence,
        "contrast": contrast,
        "stability": stability,
        "predicted_fc": pred
    }

