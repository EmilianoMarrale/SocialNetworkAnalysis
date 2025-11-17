import requests
import config.config as config


def get_signor_interaction_partners(ids: list | str) -> list[dict]:
    """Fetch interaction partners from SIGNOR API.

    Args:
        ids (list | str): Gene IDs to query.

    Returns:
        str: A list of interaction partners with all the related informations.
    """
    result_headers = ["ENTITYA", "TYPEA", "IDA", "DATABASEA", "ENTITYB", "TYPEB", "IDB", "DATABASEB", "EFFECT", 
                      "MECHANISM", "RESIDUE", "SEQUENCE", "TAX_ID", "CELL_DATA", "TISSUE_DATA", "MODULATOR_COMPLEX", "TARGET_COMPLEX",
                        "MODIFICATIONA", "MODASEQ", "MODIFICATIONB", "MODBSEQ", "PMID", "DIRECT", "NOTES", "ANNOTATOR", "SENTENCE", "SIGNOR_ID", "SCORE"]
    
    interaction_partners = []
    params = {
        "proteins": "%0d".join(ids) if isinstance(ids, list) else ids,
        "organism": config.MOUSE_SPECIES,
        "type": "connect"
    }

    response = requests.get(config.SIGNOR_INTERACTION_PARTNERS_ENDPOINT, params=params)

    if response.status_code != 200:
        print(f"Error fetching data from SIGNOR API: {response.status_code} - {response.text} - {ids}. Giving up.")
        return interaction_partners
    
    if response.status_code == 200 and response.text == 'No result found.':
        params = {
            "proteins": "%0d".join(ids) if isinstance(ids, list) else ids,
            "organism": config.HUMAN_SPECIES,
            "type": "connect"
        }
        response = requests.get(config.SIGNOR_INTERACTION_PARTNERS_ENDPOINT, params=params)
        
    if response.status_code == 200 and response.text == 'No result found.':
        print("No interaction partners found for this specific protein ", ids, "\n")
        return interaction_partners
    
    for line in response.text.strip().split("\n"):
        request_result = line.strip().split("\t")
        if(len(request_result) < len(result_headers)):
            print("Invalid response format, number mismatch between headers and data: ", line)
            continue
        result_dict = {}
        for i in range(len(request_result)):
            result_dict[result_headers[i]] = request_result[i]

        interaction_partners.append(result_dict)
        
    filtered_partners = [interaction_partner for interaction_partner in interaction_partners if float(interaction_partner["SCORE"]) > 0.5]

    return filtered_partners