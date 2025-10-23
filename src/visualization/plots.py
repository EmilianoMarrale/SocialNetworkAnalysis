import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import powerlaw

def plot_dist(G):
    # Use the non-deprecated method
    M = nx.to_scipy_sparse_array(G, dtype=float)

    M[M > 0] = 1  # Binarize the adjacency matrix
    indegrees = np.array(M.sum(axis=0)).flatten()
    degree_counts = np.bincount(indegrees.astype(int))
    
    # Remove zero counts to avoid log-scale errors
    nonzero_idx = degree_counts > 0
    degree_values = np.arange(len(degree_counts))[nonzero_idx]
    degree_counts_nonzero = degree_counts[nonzero_idx]
    
    # Fit powerlaw only if there are enough unique positive values
    if np.sum(degree_counts_nonzero > 0) >= 2:
        fit = powerlaw.Fit(degree_counts_nonzero, discrete=True, verbose=False)
    else:
        fit = None
        print("Not enough unique positive values to fit a power law.")
    
    fig = plt.figure(figsize=(16, 6))
    
    # Plot Degree Distribution
    plt.subplot(1, 3, 1)
    plt.plot(degree_values, degree_counts_nonzero, 'b.')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Degree')
    plt.ylabel('P(k)')
    
    # Plot CDF
    plt.subplot(1, 3, 2)
    if fit:
        fit.plot_cdf()
    else:
        plt.text(0.5, 0.5, "CDF not available", ha='center', va='center')
    plt.xlabel("Degree")
    plt.ylabel('CDF')
    
    # Plot CCDF
    plt.subplot(1, 3, 3)
    if fit:
        fit.plot_ccdf()
    else:
        plt.text(0.5, 0.5, "CCDF not available", ha='center', va='center')
    plt.xlabel('Degree')
    plt.ylabel('CCDF')
    
    plt.tight_layout()
    plt.show()
