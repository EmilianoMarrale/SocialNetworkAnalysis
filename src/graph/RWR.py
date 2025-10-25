import networkx as nx
import numpy as np
import pandas as pd

""" 
    Custom Module for performing Random Walk with Restart (RWR) on graphs in order to check wether fold change diffusion process yields good results.
    This module provides functions to perform RWR on a given graph and analyze the results.
"""


def random_walk_with_restart(adjacency_matrix: np.ndarray, restart_probability: float = 0.2, n_iter: int = 100, fc_vector: pd.Series = pd.Series(), diffusion_weight: float = 0.5) -> pd.Series:
    """Performs diffusion using Random Walk with Restart (RWR) on the given graph.

    Args:
        adjacency_matrix (np.ndarray): The adjacency matrix representation of the input graph to perform diffusion on.
        restart_probability (float, optional): The probability of restarting the walk, thus concluding the current interation. Defaults to 0.2.
        n_iter (int, optional): The number of diffusions to perform. The diffusion process ends once restart probability is met or all nodes have been visited. Defaults to 100.
        fc_vector (pd.Series, optional): Fold-Change vector, contains for each node its fold-change value. Defaults to pd.Series().

    Returns:
        pd.Series: The vector containing the diffused fold-change values for each node.
    """
    for iter in range(n_iter):
        diffusion_starting_node_index = randomize_starting_node(fc_vector)
        fc_vector = perform_diffusion(set(), diffusion_starting_node_index, adjacency_matrix, restart_probability, fc_vector, diffusion_weight)
    return fc_vector



def perform_diffusion(visited_nodes: set, node_to_visit: int, adjacency_matrix: np.ndarray, restart_probability: float, fc_vector: pd.Series, diffusion_weight: float) -> pd.Series:
    """Performs diffusion recursively from the given node to its neighbors until the path is exhausted or the restart probability is met.

    Args:
        visited_nodes (set): A set of nodes that have already been visited.
        node_to_visit (int): The index of the node to visit.
        adjacency_matrix (np.ndarray): The adjacency matrix representation of the graph.
        restart_probability (float): The probability of restarting the walk.
        fc_vector (pd.Series): The fold-change vector.
        diffusion_weight (float): The weight of the diffusion process.

    Returns:
        pd.Series: The updated fold-change vector after diffusion.
    """
    
    print("\n-------------------------------------------------- \n")
    # Transfer fold-change value to neighbors
    print("Selected node to visit:", node_to_visit)
    node_neighbors_vector = adjacency_matrix[:, node_to_visit].squeeze()
    print("Fold change vector: ", list(fc_vector))
    print("Selected node neighbors vector:", node_neighbors_vector)
    diffusion_result = np.dot(fc_vector.T, np.dot(node_neighbors_vector, diffusion_weight))
    print("Diffusion result: ", diffusion_result)
    fc_vector[node_to_visit] += diffusion_result
    print("Updated fold change vector: ", list(fc_vector))
    visited_nodes.add(node_to_visit)
    print("Updated set of Visited nodes:", visited_nodes)

    # Select the next node to visit and recursively perform diffusion
    node_neighbors_indices = get_neighbors_indices(node_neighbors_vector)
    print("Selected node neighbors indices:", node_neighbors_indices)
    node_neighbors_indices = (set(node_neighbors_indices) - visited_nodes)
    print("Unvisited neighbors indices:", node_neighbors_indices)
    print("\n-------------------------------------------------- \n")

    if len(node_neighbors_indices) == 0: # I have finished the nodes i can visit on this path.
        print("\n-------------------------------------------------- \n")
        print("\nNO MORE NEIGHBORS TO VISIT, ENDING DIFFUSION\n")
        print("\n-------------------------------------------------- \n")
        return fc_vector
    
    if np.random.rand() < restart_probability:
        print("\n-------------------------------------------------- \n")
        print("\nRESTARTING DIFFUSION\n")
        print("\n-------------------------------------------------- \n")
        return fc_vector
    
    next_node = np.random.choice(list(node_neighbors_indices))
    return perform_diffusion(visited_nodes, next_node, adjacency_matrix, restart_probability, fc_vector, diffusion_weight)        
    

def randomize_starting_node(fc_vector: pd.Series) -> int:
    """Randomly selects a starting node with a non-zero fold-change value for the diffusion process.

    Args:
        fc_vector (pd.Series): The fold-change vector.

    Returns:
        int: The index of the randomly selected starting node.
    """
    fc_vector = pd.Series({idx: value for idx, value in fc_vector.items() if value != 0})
    node_index = np.random.randint(0, len(fc_vector))
    return node_index


def get_neighbors_indices(node_column: np.ndarray) -> list[int]:
    """Returns the indices of the neighboring nodes for a given node.

    Args:
        node_column (np.ndarray): The column vector representing the node's connections.
        node_index (int): The index of the node for which to find neighbors.

    Returns:
        list[int]: A list of indices representing the neighboring nodes.
    """
    neighbors_indices = list(np.where(node_column > 0)[0])
    return neighbors_indices