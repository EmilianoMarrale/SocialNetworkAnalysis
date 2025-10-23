import requests
import config.config as config


def get_string_interaction_partners(ids: list | str) -> str:
    """Fetch interaction partners from STRING API.

    Args:
        ids (list | str): Gene IDs to query.

    Returns:
        str: A list of interaction partners with their combined confidence score.
    """

    interaction_partners = []
    params = {
        "identifiers": "%0d".join(ids) if isinstance(ids, list) else ids,
        "species": config.MOUSE_SPECIES,
        "required_score": 700,  # Soglia di confidenza
    }

    response = requests.post(config.STRING_INTERACTION_PARTNERS_ENDPOINT, data=params)

    # if response.status_code != 200:
    #     print(f"Error fetching data from STRING API: {response.status_code} - {response.text} - {ids}. Retrying with MOUSE taxonomy.")
    #     params["species"] = config.MOUSE_SPECIES
    #     response = requests.post(config.STRING_REQUEST_ENDPOINT_URL, data=params)

    if response.status_code != 200:
        print(f"Error fetching data from STRING API: {response.status_code} - {response.text} - {ids}. Giving up.")
        return interaction_partners
    
    for line in response.text.strip().split("\n"):
        l = line.strip().split("\t")
        if len(l) < 6:
            continue
        partner_name = l[3]
        combined_score = l[5]
        interaction_partners.append({"name": partner_name, "combined_score": combined_score})

    return interaction_partners

def get_gene_set_enrichment(ids: list | str, background_ids: list | str) -> dict:
    """Fetch gene set enrichment from STRING API.

    Args:
        ids (list | str): Gene IDs to query.
        background_ids (list | str): Background gene IDs.

    Returns:
        dict: Enrichment results.
    """

    params = {
        "identifiers": "%0d".join(ids) if isinstance(ids, list) else ids,
        "background_identifiers": "%0d".join(background_ids) if isinstance(background_ids, list) else background_ids,
        "species": config.MOUSE_SPECIES,
        "caller_identity": "network_analysis_project",
    }

    response = requests.post(config.STRING_ENRICHMENT_ENDPOINT, data=params)

    if response.status_code != 200:
        print(f"Error fetching enrichment data from STRING API: {response.status_code} - {response.text} - {ids}.")
        return None

    return response.json()


def get_enrichment_figure(ids: list | str, category: list | str) -> bytes:
    """Fetch enrichment figure from STRING API.

    Args:
        ids (list | str): Gene IDs to query.
        category (list | str): Enrichment categories.

    Returns:
        bytes: Enrichment figure in bytes.
    """

    params = {
        "identifiers": "%0d".join(ids) if isinstance(ids, list) else ids,
        "category": "%0d".join(category) if isinstance(category, list) else category,
        "species": config.MOUSE_SPECIES,
        "caller_identity": "network_analysis_project",
    }

    response = requests.post(config.STRING_FIGURE_ENRICHMENT_ENDPOINT, data=params)

    if response.status_code != 200:
        print(f"Error fetching enrichment figure from STRING API: {response.status_code} - {response.text} - {ids}.")
        return None

    return response.content