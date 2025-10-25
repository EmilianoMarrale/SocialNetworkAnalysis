import json
from src.data_collection import stringdb_api
import networkx as nx
from config import config
from src.visualization.plots import plot_dist
from src.data_collection.classes import GeneNode
import pandas as pd

def build_graph_from_api(seed_gene_list: list[str], graph_name: str):
    """Build a gene interaction graph from an external API and save it as a GraphML file.
    Args:
        seed_gene_list (list[dict]): List of gene dictionaries to build the graph.
        graph_name (str): Name of the graph file (without extension).
    Returns:
        nx.Graph: The constructed gene interaction graph.
    """
    graph = nx.Graph()

    for seed_gene in seed_gene_list:
        print(f"Processing seed gene: {seed_gene['name']}")
        seed_node = GeneNode.from_dict(seed_gene)
        
        # Add seed node with attributes unpacked
        graph.add_node(seed_node.name, **seed_node.to_dict()) # I always want to add the seed node since it has differential expression info.
        print("added seed node", seed_node.name, seed_node.fold_change, seed_node.p_value, seed_node.padj)

        first_degree_neighbors = stringdb_api.get_string_interaction_partners(seed_node.name)
        for first_neighbor in first_degree_neighbors:
            first_neighbor_node = GeneNode.from_dict(first_neighbor)
            if not graph.has_node(first_neighbor_node.name): # Avoid duplicating nodes, we might overwrite some seed node otherwise.
                graph.add_node(first_neighbor_node.name, **first_neighbor_node.to_dict())
            graph.add_edge(seed_node.name, first_neighbor_node.name, weight=float(first_neighbor_node.combined_score))

            second_degree_neighbors = stringdb_api.get_string_interaction_partners(first_neighbor_node.name)
            for second_neighbor in second_degree_neighbors:
                second_neighbor_node = GeneNode.from_dict(second_neighbor)
                if not graph.has_node(second_neighbor_node.name): # Avoid duplicating nodes, we might overwrite some seed node otherwise.
                    graph.add_node(second_neighbor_node.name, **second_neighbor_node.to_dict())
                graph.add_edge(first_neighbor_node.name, second_neighbor_node.name, weight=float(second_neighbor_node.combined_score))

    print(f"Gene graph: {graph.number_of_nodes()} nodes, {graph.number_of_edges()} edges")
    nx.write_graphml(graph, config.GRAPHS_DIR / f"{graph_name}.graphml")

    return graph



def get_gene_graph(gene_list: list[dict], graph_name: str) -> nx.Graph:

    """Load a graph from a GraphML file if exists, otherwise build the graph integrating first and second level neighbors of gene_list using STRING api calls.
    
    Args:
        gene_list (list[dict]): List of gene dictionaries to build the graph if not found.
        graph_name (str): Name of the graph file (without extension).
    Returns:
        nx.Graph: The loaded or newly built graph.

    """
    graph_path = config.GRAPHS_DIR / f"{graph_name}.graphml"
    if not graph_path.exists():
        print(f"Graph file {graph_path} does not exist.")
        return build_graph_from_api(gene_list, graph_name)

    graph = nx.read_graphml(graph_path)

    convert_edge_weights_to_float(graph)

    print(f"Loaded graph '{graph_name}': {graph.number_of_nodes()} nodes, {graph.number_of_edges()} edges")
    return graph



def get_connected_seed_genes_graph(seed_genes: pd.DataFrame, graph_name: str) -> nx.Graph:

    graph_path = config.GRAPHS_DIR / f"{graph_name}.graphml"
    if graph_path.exists():
        graph = nx.read_graphml(graph_path)
        convert_edge_weights_to_float(graph)
        print(f"Loaded graph '{graph_name}': {graph.number_of_nodes()} nodes, {graph.number_of_edges()} edges")
        return graph

    graph = nx.Graph()

    for seed_gene in seed_genes.to_dict(orient="records"):
        print(f"Processing seed gene: {seed_gene['name']}")
        seed_node = GeneNode.from_dict(seed_gene)
        
        # Add seed node with attributes unpacked
        graph.add_node(seed_node.name, **seed_node.to_dict()) # I always want to add the seed node since it has differential expression info.
        print("added seed node", seed_node.name, seed_node.fold_change, seed_node.p_value, seed_node.padj)

        first_degree_neighbors = stringdb_api.get_string_interaction_partners(seed_node.name)

        for first_neighbor in first_degree_neighbors:
            first_neighbor_node = GeneNode.from_dict(first_neighbor)

            if first_neighbor in seed_genes["name"].values:

                if not graph.has_node(first_neighbor_node.name): # Avoid duplicating nodes, we might overwrite some seed node otherwise.
                    graph.add_node(first_neighbor_node.name, **first_neighbor_node.to_dict())

                graph.add_edge(seed_node.name, first_neighbor_node.name, weight=float(first_neighbor_node.combined_score))

    print(f"Gene graph: {graph.number_of_nodes()} nodes, {graph.number_of_edges()} edges")
    nx.write_graphml(graph, config.GRAPHS_DIR / f"{graph_name}.graphml")
    return graph



def convert_edge_weights_to_float(graph: nx.Graph):
    """ When we load a graph from GraphML file, the edge weights are read as strings.
    This function converts all edge weights to float using Python’s dictionary comprehension, which is much faster than looping over the edges.
    This is necessary for numerical computations.

    Args:
        graph (nx.Graph): The input graph.
    """
    # Get all weights as a dictionary
    weights = nx.get_edge_attributes(graph, "weight")

    # Convert all values to float at once
    weights = {k: float(v) for k, v in weights.items()}

    # Set them back to the graph
    nx.set_edge_attributes(graph, weights, "weight")

def print_analysis(graph: nx.Graph) -> dict:
    """Print and returns basic analysis of the graph.

    Args:
        graph (nx.Graph): The input graph.

    Returns:
        dict: A dictionary containing various analysis metrics.
        
    """

    analysis = {}
    print(f"Graph has {graph.number_of_nodes()} nodes and {graph.number_of_edges()} edges.\n")
    degrees = [d for n, d in graph.degree()]
    if degrees:
        analysis["avg_degree"] = sum(degrees) / len(degrees)
        analysis["max_degree"] = max(degrees)
        analysis["min_degree"] = min(degrees)

        print("Degree analysis:")
        print(f"\tAverage degree: {analysis['avg_degree']:.2f}")
        print(f"\tMax degree: {analysis['max_degree']}")
        print(f"\tMin degree: {analysis['min_degree']}")

        analysis["num_cc"] = nx.number_connected_components(graph)
        analysis["largest_cc"] = max(nx.connected_components(graph), key=len)
        print("\nConnectedness analysis:")
        print(f"\tNumber of connected components: {analysis['num_cc']}")
        print(f"\tSize of largest connected component: {len(analysis['largest_cc'])}")
        print(f"\tIs the graph connected? {nx.is_connected(graph)}")

        analysis["lcc_diameter"] = nx.diameter(graph.subgraph(analysis["largest_cc"]))
        analysis["lcc_avg_shortest_path"] = nx.average_shortest_path_length(graph.subgraph(analysis["largest_cc"]))
        print("\nLargest Component Path analysis:")
        print(f"\tDiameter of largest connected component: {analysis['lcc_diameter']}")
        print(f"\tAverage shortest path length of largest connected component: {analysis['lcc_avg_shortest_path']}")

        analysis["avg_clustering"] = nx.average_clustering(graph)
        analysis["density"] = nx.density(graph)
        print("\nClustering and Density analysis:")
        print(f"\tAverage clustering coefficient: {analysis['avg_clustering']}")
        print(f"\tGraph density: {analysis['density']}\n")

        analysis["betweenness_centrality"] = nx.betweenness_centrality(graph)
        print("\nCentrality analysis:")
        print(f"\tAverage betweenness centrality: {sum(analysis['betweenness_centrality'].values()) / len(graph)}")
        print(f"\tMax betweenness centrality: {max(analysis['betweenness_centrality'].values())}")
        print(f"\tMin betweenness centrality: {min(analysis['betweenness_centrality'].values())}")

        analysis["closeness_centrality"] = nx.closeness_centrality(graph)
        print(f"\tAverage closeness centrality: {sum(analysis['closeness_centrality'].values()) / len(graph)}")
        print(f"\tMax closeness centrality: {max(analysis['closeness_centrality'].values())}")
        print(f"\tMin closeness centrality: {min(analysis['closeness_centrality'].values())}")

    else:
        print("The graph has no nodes.")

    return analysis
    
def get_sorted_degrees(graph: nx.Graph, direction: str = "desc"):
    """Get the degrees of all nodes in the graph, sorted in descending order.

    Args:
        graph (nx.Graph): The input graph.
        direction (str, optional): The direction to sort the degrees. Can be "asc" for ascending or "desc" for descending. Defaults to "desc".

    Returns:
        List[Tuple[str, int]]: A list of tuples where each tuple contains a node and its degree, sorted by degree in descending order.
    """
    if direction not in ["asc", "desc"]:
        raise ValueError("Direction must be 'asc' or 'desc'")

    reverse = direction == "desc"
    degrees = dict(graph.degree())
    sorted_degrees = sorted(degrees.items(), key=lambda x: x[1], reverse=reverse)
    return sorted_degrees

def get_nodes_with_certain_degree(graph: nx.Graph, degree: int | tuple[int]):
    """Get nodes with a certain degree from the graph.

    Args:
        graph (nx.Graph): The input graph.
        degree (int | tuple(int,int)): The degree or range of degrees to filter nodes.

    Returns:
        List[str]: A list of node names that match the specified degree.
    """
    if isinstance(degree, int):
        return [n for n, d in graph.degree() if d == degree]
    elif isinstance(degree, tuple) and len(degree) == 2:
        return [n for n, d in graph.degree() if degree[0] <= d <= degree[1]]
    else:
        raise ValueError("Degree must be an int or a tuple of two ints.")
    

def get_ER_graph(n_nodes :int, n_edges :int):
    """Generate an Erdős-Rényi (ER) random graph with a specified number of nodes and edges.

    Args:
        n_nodes (int): number of nodes
        n_edges (int): number of edges

    Returns:
        _type_: nx.Graph: The generated ER graph.
    """
    # The expected number of edges is m = p * (n * n-1) / 2, thus we solve for p and have the following
    p = (2 * n_edges) / (n_nodes * (n_nodes - 1))
    return nx.erdos_renyi_graph(n_nodes, p)


def get_BA_graph(n_nodes :int, n_edges_per_node :int):
    """Generate a Barabási-Albert (BA) random graph with a specified number of nodes and edges per new node.

    Args:
        n_nodes (int): number of nodes
        n_edges_per_node (int): number of edges to attach from a new node to existing nodes

    Returns:
        _type_: nx.Graph: The generated BA graph.
    """
    return nx.barabasi_albert_graph(n_nodes, n_edges_per_node)