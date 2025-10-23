from importlib.metadata import metadata
import os
import pandas as pd
import numpy as np
from config.config import DEG_FILE, RAW_BULK_RNA_SEQ_FILE, PROCESSED_BULK_RNA_SEQ_FILE

from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

""" Differentially expressed genes (DEGs) are all the genes that show a statistically significant difference in expression levels between two conditions or groups.
    DEGs can be either upregulated or downregulated:
    - Upregulated Genes: These are genes that show an increase in expression levels in one condition compared to another. For example, if a gene is more active in diseased tissue compared to healthy tissue, it is considered upregulated in the disease state.
    - Downregulated Genes: Conversely, these are genes that show a decrease in expression levels in one condition compared to another. For instance, if a gene is less active in treated cells compared to untreated cells, it is considered downregulated in the treated state.
    Identifying DEGs is crucial for understanding the molecular mechanisms underlying various biological processes and diseases. They can provide insights into pathways that are activated or suppressed in response to specific conditions, treatments, or environmental factors. 
    Researchers often use DEGs to identify potential biomarkers for diseases, therapeutic targets, and to gain a deeper understanding of cellular functions and regulatory networks.
"""


def preprocess_bulk_rna_seq_data(padj_threshold :float, log2fc_threshold: float) -> pd.DataFrame:
    """ Preprocess the RNA-seq data and return a dataframe containing the significant differentially expressed genes

    Args:
        padj_threshold (float): The adjusted p-value threshold for significance.
        log2fc_threshold (float): The log2 fold change threshold for significance. |l2fc| = 1 means a 2-fold change (2x increase or decrease).

    Returns:
        pd.DataFrame: Dataframe containing the significant differentially expressed genes
    """
    if os.path.exists(PROCESSED_BULK_RNA_SEQ_FILE):
        results = pd.read_csv(PROCESSED_BULK_RNA_SEQ_FILE)
        return results[(results["padj"] <= padj_threshold) & (abs(results["log2FoldChange"]) > log2fc_threshold)]

    df = pd.read_excel(RAW_BULK_RNA_SEQ_FILE)
    columns_to_drop = df.columns.difference(['name', 'Young.1.liver_COUNT', 'Young.2.liver_COUNT', 'Young.3.liver_COUNT', 'Aging.1.liver_COUNT', 'Aging.2.liver_COUNT', 'Aging.3.liver_COUNT', '302b.1.liver_COUNT', '302b.2.liver_COUNT', '302b.3.liver_COUNT'])
    df.drop(columns=columns_to_drop, axis=1, inplace=True)
    dft = df.T
    new_header = dft.iloc[0]   # first row as header
    dft = dft[1:]               # remove first row from data
    dft.columns = new_header   # assign new header

    # Build the metadata dataframe needed for DESeq2
    metadata = pd.DataFrame(index=np.arange(dft.shape[0]), columns=["sample", "condition"])
    metadata["sample"] = dft.index
    metadata["group"] = (["Young"] * 3 + ["Aging"] * 3 + ["302b"] * 3)
    metadata["condition"] = (["Young.1.liver_COUNT"] + ["Young.2.liver_COUNT"] + ["Young.3.liver_COUNT"] +
                            ["Aging.1.liver_COUNT"] + ["Aging.2.liver_COUNT"] + ["Aging.3.liver_COUNT"] +
                            ["302b.1.liver_COUNT"] + ["302b.2.liver_COUNT"] + ["302b.3.liver_COUNT"])
    metadata.set_index("sample", inplace=True)

    # filter columns with counts less than 10 in all samples
    filtered_df = dft.loc[:, (dft >= 10).any(axis=0)]
    # Remove nan columns
    filtered_df = filtered_df.loc[:, ~filtered_df.columns.duplicated(keep='first')]

    inference = DefaultInference()
    dds = DeseqDataSet(counts=filtered_df, 
                    metadata=metadata,
                    design ="~group",
                    ref_level={"group": "Young"},
                    refit_cooks=True,
                    inference=inference,
                    )
    dds.deseq2()
    ds_B_vs_A = DeseqStats(dds, contrast=["group", "302b", "Aging"], inference=inference)
    ds_B_vs_A.summary()
    results = ds_B_vs_A.results_df
    results.to_csv(PROCESSED_BULK_RNA_SEQ_FILE, index=True)
    # Filter as specified in the paper and return the results
    return results[(results["padj"] <= padj_threshold) & (abs(results["log2FoldChange"]) > log2fc_threshold)]


def get_downregulated_genes(df :pd.DataFrame) -> list:
    """
    Extract downregulated genes from the DEG data.

    Returns:
    - list: A DataFrame of downregulated genes.
    """
    return df[df["log2FoldChange"] < 0]

def get_upregulated_genes(df :pd.DataFrame) -> list:
    """
    Extract upregulated genes from the DEG data.

    Returns:
    - list: A DataFrame of upregulated genes.
    """
    return df[df["log2FoldChange"] > 0]

def get_downregulated_genes_tissue(tissue: str = None) -> list:
    """
    Extract downregulated genes from the DEG data.

    Returns:
    - list: A list of downregulated genes.
    """
    xl_file = pd.read_excel(DEG_FILE, sheet_name=None)
    downregulated_genes = []
    if tissue == "liver":
        downregulated_genes.append(xl_file["liver"]["Aging-Exos vs Aging-pbs down"])
    elif tissue == "skin":
        downregulated_genes.append(xl_file["skin"]["Aging-Exos vs Aging-pbs down"])
    else:
        downregulated_genes.append(xl_file["liver"]["Aging-Exos vs Aging-pbs down"])
        downregulated_genes.append(xl_file["skin"]["Aging-Exos vs Aging-pbs down"])

    downregulated_genes = pd.concat(downregulated_genes).dropna().unique().tolist()
    return downregulated_genes

# def get_upregulated_genes(tissue: str = None) -> list:
#     """
#     Extract upregulated genes from the DEG data.

#     Returns:
#     - list: A list of upregulated genes.
#     """
#     xl_file = pd.read_excel(DEG_FILE, sheet_name=None)
#     upregulated_genes = []
#     if tissue == "liver":
#         upregulated_genes.append(xl_file["liver"]["Aging-Exos vs Aging-pbs up"])
#     elif tissue == "skin":
#         upregulated_genes.append(xl_file["skin"]["Aging-Exos vs Aging-pbs up"])
#     else:
#         upregulated_genes.append(xl_file["liver"]["Aging-Exos vs Aging-pbs up"])
#         upregulated_genes.append(xl_file["skin"]["Aging-Exos vs Aging-pbs up"])

#     upregulated_genes = pd.concat(upregulated_genes).dropna().unique().tolist()
#     return upregulated_genes
