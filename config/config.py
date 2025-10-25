# config/config.py

import os
from pathlib import Path

# === PROJECT ROOT ===
# Ottiene la root del progetto (parent di config/)
PROJECT_ROOT = Path(__file__).parent.parent

# === DATA PATHS ===
DATA_DIR = PROJECT_ROOT / 'data'
RAW_DATA_DIR = DATA_DIR / 'raw'
PROCESSED_DATA_DIR = DATA_DIR / 'processed'
GRAPHS_DIR = DATA_DIR / 'graphs'

# File specifici
DEG_FILE = RAW_DATA_DIR / 'DEG-TS6.xlsx' # Effects of hESC-Exo treatment on cell proliferation in vivo. scRNA-seq on the liver and skin tissues isolated from aging-Exos and aging-PBS mice
NETWORK_FILE = RAW_DATA_DIR / 'network_data.csv'

RAW_BULK_RNA_SEQ_FILE = RAW_DATA_DIR / 'Bulk-RNA-seq-TS1.xlsx'
aging_vs_302_PROCESSED_BULK_RNA_SEQ_FILE = PROCESSED_DATA_DIR / 'aging_vs_302_significant_DEGs.csv'
young_vs_aging_PROCESSED_BULK_RNA_SEQ_FILE = PROCESSED_DATA_DIR / 'young_vs_Aging_significant_DEGs.csv'
WHOLE_DOWNREGULATED_GRAPH_FILE_NAME = GRAPHS_DIR / 'downregulated_graph'
WHOLE_UPREGULATED_GRAPH_FILE_NAME = GRAPHS_DIR / 'upregulated_graph'
LIVER_DOWNREGULATED_GRAPH_FILE_NAME = GRAPHS_DIR / 'liver_downregulated_graph'
LIVER_UPREGULATED_GRAPH_FILE_NAME = GRAPHS_DIR / 'liver_upregulated_graph'

# === RESULTS PATHS ===
RESULTS_DIR = PROJECT_ROOT / 'results'
FIGURES_DIR = RESULTS_DIR / 'figures'
METRICS_DIR = RESULTS_DIR / 'metrics'

# === API CONFIGURATION ===
STRING_API_URL = "https://version-12-0.string-db.org/api"

TSV_OUTPUT_FORMAT = "tsv-no-header"
JSON_OUTPUT_FORMAT = "json"
IMAGE_OUTPUT = "highres_image"

INTERACTION_METHOD = "interaction_partners"
ENRICHMENT_METHOD = "enrichment"
FIGURE_ENRICHMENT_METHOD = "enrichmentfigure"

HUMAN_SPECIES = "9606"  # Homo sapiens (human)
MOUSE_SPECIES = "10090"  # Mus musculus (mouse)

STRING_INTERACTION_PARTNERS_ENDPOINT = "/".join([STRING_API_URL, TSV_OUTPUT_FORMAT, INTERACTION_METHOD])
STRING_ENRICHMENT_ENDPOINT = "/".join([STRING_API_URL, JSON_OUTPUT_FORMAT, ENRICHMENT_METHOD])
STRING_FIGURE_ENRICHMENT_ENDPOINT = "/".join([STRING_API_URL, IMAGE_OUTPUT, FIGURE_ENRICHMENT_METHOD])


# === GRAPH PARAMETERS ===
GRAPH_TYPE = 'undirected'  # o 'directed'
MIN_NODE_DEGREE = 1
WEIGHT_THRESHOLD = 0.5

# === ANALYSIS PARAMETERS ===
# Centrality
CENTRALITY_METRICS = ['degree', 'betweenness', 'closeness', 'eigenvector']

# Community Detection
COMMUNITY_ALGORITHM = 'louvain'  # louvain, girvan_newman, label_propagation
COMMUNITY_RESOLUTION = 1.0

# Link Prediction
LINK_PREDICTION_METHODS = ['common_neighbors', 'adamic_adar', 'preferential_attachment']
LINK_PREDICTION_TEST_SIZE = 0.2

# === VISUALIZATION PARAMETERS ===
FIGURE_DPI = 300
FIGURE_FORMAT = 'png'
NODE_SIZE_RANGE = (50, 500)
EDGE_WIDTH_RANGE = (0.5, 3)

# === FUNCTION TO CREATE DIRECTORIES ===
def create_directories():
    """Crea tutte le directory necessarie se non esistono"""
    directories = [
        RAW_DATA_DIR,
        PROCESSED_DATA_DIR,
        GRAPHS_DIR,
        FIGURES_DIR,
        METRICS_DIR
    ]
    for directory in directories:
        directory.mkdir(parents=True, exist_ok=True)
    print("✓ Directories created successfully")

# === FUNCTION TO VALIDATE PATHS ===
def validate_data_files():
    """Verifica che i file di dati esistano"""
    required_files = [DEG_FILE]
    
    missing_files = []
    for file_path in required_files:
        if not file_path.exists():
            missing_files.append(str(file_path))
    
    if missing_files:
        raise FileNotFoundError(
            f"Missing required data files:\n" + "\n".join(missing_files)
        )
    print("✓ All required data files found")

# Crea le directory all'import (opzionale)
if __name__ == "__main__":
    create_directories()
    print(f"Project root: {PROJECT_ROOT}")
    print(f"DEG file: {DEG_FILE}")