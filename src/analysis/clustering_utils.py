from cdlib import NodeClustering
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from src.data_collection import stringdb_api

def get_partitioning_metadata(partitioning: NodeClustering, background_genes: list, verbose: bool = False, debug: bool = False) -> dict:
    """Extracts metadata from a NodeClustering object.

    Args:
        communities (NodeClustering): The clustering result.
        background_genes (list): List of background genes for enrichment analysis. (The whole network in my case)
        verbose (bool): If True, prints detailed information.

    Returns:
        dict: A dictionary containing metadata about the clustering.
                - num_communities: Number of communities detected.
                - avg_community_size: Average size of the communities.
                - internal_edge_density: Internal edge density of the partitioning.
                - average_internal_degree: Average internal degree of the partitioning.
                - modularity_density: Modularity density of the partitioning.
                - conductance: Conductance of the partitioning.
                - community_enrichments: List of best enrichment results for each community, in the following form:
                    {
                        "community_index": int,
                        "best_kegg": dict or None,
                        "best_function": dict or None
                    }
                    where best_kegg and best_function are dictionaries containing enrichment results, in the following form:
                    {
                        "category": str, # either "KEGG" or "Function"
                        "term": str, # the enriched term
                        "number_of_genes": int, # number of genes in the community associated with the term
                        "number_of_genes_in_background": int, # number of genes in the background associated with the term
                        "enrichment_ratio": float, # ratio of number_of_genes to number_of_genes_in_background
                        "p_value": float, # p-value of the enrichment
                        "fdr": float, # false discovery rate of the enrichment
                        "description": str # description of the enriched term
                    }                   
                - distinct_kegg_descriptions_count: dict containing the count of distinct KEGG descriptions across communities.
                - distinct_function_descriptions_count: dict containing the count of distinct Function descriptions across communities.
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
        if debug and i == 50:
            break
        print("Community: ", i , "/", len(partitioning.communities), end="\r")

    metadata["community_enrichments"] = community_enrichments
    metadata["distinct_kegg_descriptions_count"] = distinct_kegg_descriptions_count
    metadata["distinct_function_descriptions_count"] = distinct_function_descriptions_count

    if verbose:
        print_partitioning_metadata(metadata)

    return metadata

def average_enrichment_coverage(metadata: dict) -> tuple[float, float]:
    """Calculates the average enrichment coverage across all communities.

    Args:
        metadata (dict): The metadata dictionary.

    Returns:
        (float, float): The average enrichment coverage for KEGG and Function terms respectively.
    """
    kegg_coverage_sum = 0
    kegg_num = 0

    function_coverage_sum = 0
    function_num = 0
    for com_enrich in metadata["community_enrichments"]:
        if com_enrich["best_kegg"] is not None:
            kegg_coverage_sum += com_enrich["best_kegg"]["enrichment_ratio"]
            kegg_num += 1
        if com_enrich["best_function"] is not None:
            function_coverage_sum += com_enrich["best_function"]["enrichment_ratio"]
            function_num += 1

    avg_kegg_coverage = kegg_coverage_sum / kegg_num if kegg_num > 0 else 0.0
    avg_function_coverage = function_coverage_sum / function_num if function_num > 0 else 0.0

    return avg_kegg_coverage, avg_function_coverage

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
    plt.subplots(figsize=(12, 6))

    plt.bar(metadata["distinct_kegg_descriptions_count"].keys(),
            metadata["distinct_kegg_descriptions_count"].values(),
            color='skyblue', label='KEGG Descriptions')
    plt.bar(metadata["distinct_function_descriptions_count"].keys(),
            metadata["distinct_function_descriptions_count"].values(),
            color='salmon', label='Function Descriptions')
    plt.xlabel('Descriptions')
    plt.ylabel('Count')
    plt.title('Distinct KEGG and Function Descriptions Count')
    plt.xticks(rotation=45, ha='right')
    plt.legend()
    plt.tight_layout()
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

def score_partitionings(df: pd.DataFrame, weights: dict) -> pd.DataFrame:
    """Scores partitionings based on provided weights for each metric.

    Args:
        df (pd.DataFrame): DataFrame containing partitioning metrics.
        weights (dict): Weights for each metric.

    Returns:
        pd.DataFrame: DataFrame with an additional 'score' column.
    """
    # Normalize metrics
    normalized_df = df.copy()
    for column in df.columns:
        if column != "Method":
            col_min = df[column].min()
            col_max = df[column].max()
            normalized_df[column] = (df[column] - col_min) / (col_max - col_min)

    # For metrics where lower is better → invert
    for col in ['n_clusters', 'max_conductance', 'n_distinct_kegg_terms', 'n_distinct_go_terms']:
        normalized_df[col] = 1 - normalized_df[col]

    # Calculate scores
    scores = []
    for index, row in normalized_df.iterrows():
        score = 0
        for metric, weight in weights.items():
            score += row[metric] * weight
        scores.append(score)
    
    normalized_df["score"] = scores
    return normalized_df