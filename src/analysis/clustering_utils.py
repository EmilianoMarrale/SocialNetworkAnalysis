from cdlib import NodeClustering
import numpy as np
import matplotlib.pyplot as plt
from src.data_collection import stringdb_api

def get_partitioning_metadata(partitioning: NodeClustering, background_genes: list, verbose: bool = False, debug: bool = False) -> dict:
    """Extracts metadata from a NodeClustering object.

    Args:
        communities (NodeClustering): The clustering result.
        verbose (bool): If True, prints detailed information.

    Returns:
        dict: A dictionary containing metadata about the clustering.
    """
    metadata = {}

    metadata["num_communities"] = len(partitioning.communities)

    avg_community_size = np.mean([len(c) for c in partitioning.communities])
    metadata["avg_community_size"] = avg_community_size
    metadata["internal_edge_density"] = partitioning.internal_edge_density()
    metadata["average_internal_degree"] = partitioning.average_internal_degree()
    metadata["modularity_density"] = partitioning.modularity_density()
    metadata["conductance"] = partitioning.conductance()

    community_enrichments = []
    distinct_kegg_descriptions_count = {}
    distinct_function_descriptions_count = {}
    for i, community in enumerate(partitioning.communities):
        if len(community) <= 2:
            continue
        result = stringdb_api.get_gene_set_enrichment(community, background_genes)
        res_filtered = [r for r in result if r["category"] in ["KEGG", "Function"] and r["number_of_genes_in_background"] > 2]
        res_sorted = sorted(res_filtered, key=lambda x: (x['fdr'], -(x["number_of_genes"]/x["number_of_genes_in_background"])), reverse=False)
        best_kegg = next((r for r in res_sorted if r["category"] == "KEGG"), None)
        best_function = next((r for r in res_sorted if r["category"] == "Function"), None)
        community_enrichments.append({
            "community_index": i,
            "best_kegg": trim_enrichment(best_kegg) if best_kegg else None,
            "best_function": trim_enrichment(best_function) if best_function else None,
        })
        if best_kegg:
            descr = best_kegg["description"]
            if descr in distinct_kegg_descriptions_count:
                distinct_kegg_descriptions_count[descr] += 1
            else:
                distinct_kegg_descriptions_count[descr] = 1
        if best_function:
            descr = best_function["description"]
            if descr in distinct_function_descriptions_count:
                distinct_function_descriptions_count[descr] += 1
            else:
                distinct_function_descriptions_count[descr] = 1
        if debug and i == 3:
            break

    metadata["community_enrichments"] = community_enrichments
    metadata["distinct_kegg_descriptions_count"] = distinct_kegg_descriptions_count
    metadata["distinct_function_descriptions_count"] = distinct_function_descriptions_count

    if verbose:
        print_partitioning_metadata(metadata)

    return metadata

def trim_enrichment(enrichment: dict) -> dict:
    """Trims unnecessary fields from enrichment results and adds the enrichment ratio.

    Args:
        enrichment (dict): The enrichment result.

    """
    trimmed_enrichment = {"category": enrichment["category"],
                          "term": enrichment["term"],
                          "number_of_genes": enrichment["number_of_genes"],
                          "number_of_genes_in_background": enrichment["number_of_genes_in_background"],
                          "enrichment_ratio": enrichment["number_of_genes"] / enrichment["number_of_genes_in_background"],
                          "p_value": enrichment["p_value"],
                          "fdr": enrichment["fdr"],
                          "description": enrichment["description"]}

    return trimmed_enrichment


def plot_partitioning_enrichment(metadata: dict) -> None:
    """Plots the enrichment ratios of the best KEGG and Function terms for each community.

    Args:
        metadata (dict): The metadata dictionary.
    """
    plt.subplot(1,2,0)
    plt.plot(metadata["distinct_kegg_descriptions_count"].keys(), metadata["distinct_kegg_descriptions_count"].values())

    plt.subplot(1,2,1)
    plt.plot(metadata["distinct_function_descriptions_count"].keys(), metadata["distinct_function_descriptions_count"].values())
    plt.show()
        

def print_partitioning_metadata(metadata: dict) -> None:
    """Prints the metadata of a partitioning.

    Args:
        metadata (dict): The metadata dictionary.
    """
    print("Partitioning Metadata:")
    print(f"\tNumber of communities: {metadata['num_communities']}")
    print(f"\tAverage community size: {metadata['avg_community_size']:.2f}")
    print(f"\tInternal edge density: {metadata['internal_edge_density']}")
    print(f"\tAverage internal degree: {metadata['average_internal_degree']}")
    print(f"\tModularity density: {metadata['modularity_density']}")
    print(f"\tConductance: {metadata['conductance']}")
