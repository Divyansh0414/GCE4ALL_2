import requests

def search_pubchem_for_cids(smiles_file_path, output_cid_file_path):
    """
    Search PubChem for CIDs of SMILES strings from the output file.
    
    Args:
        smiles_file_path (str): Path to the text file containing SMILES strings.
        output_cid_file_path (str): Path to the output text file where CIDs will be saved.
    """
    with open(smiles_file_path, 'r') as smiles_file, open(output_cid_file_path, 'w') as output_file:
        for line in smiles_file:
            if ':' in line:
                smiles = line.split(':')[1].strip()
                try:
                    response = requests.get(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/cids/TXT")
                    if response.status_code == 200:
                        cids = response.text.strip()
                        output_file.write(f"SMILES: {smiles}, CIDs: {cids}\n")
                        print(f"SMILES: {smiles}, CIDs: {cids}")
                    else:
                        output_file.write(f"SMILES: {smiles}, CIDs: Not found\n")
                        print(f"SMILES: {smiles}, CIDs: Not found")
                except Exception as e:
                    output_file.write(f"SMILES: {smiles}, CIDs: Failed to retrieve - {e}\n")
                    print(f"Failed to retrieve CID for SMILES: {smiles} - {e}")

# Example usage
smiles_file_path = '/mnt/data/chemical_structures_smiles.txt'  # Replace with your input SMILES file
output_cid_file_path = '/mnt/data/chemical_structures_cids.txt'  # Replace with your desired output CID file
search_pubchem_for_cids(smiles_file_path, output_cid_file_path)
