import os
import re
from rdkit import Chem
from rdkit.Chem import Draw

def draw_smiles_with_rdkit(smiles_file_path, output_folder):
    """
    Draw chemical structures using RDKit by reading SMILES strings from a file and export them with the compound name.
    
    Args:
        smiles_file_path (str): Path to the text file containing SMILES strings.
        output_folder (str): Path to the folder where the images will be saved.
    """
    os.makedirs(output_folder, exist_ok=True)
    drawn_compounds = set()  # Keep track of compounds already drawn

    with open(smiles_file_path, 'r') as smiles_file:
        for line in smiles_file:
            if 'SMILES:' in line and 'CIDs:' in line and 'Not found' not in line:
                smiles = line.split('SMILES:')[1].split(', CIDs:')[0].strip()
                cid = line.split('CIDs:')[1].strip()

                # Skip entries where CID is '0'
                if cid == '0':
                    print(f"Skipping SMILES with CID 0: {smiles}")
                    continue

                # Skip duplicates
                if cid in drawn_compounds:
                    print(f"Skipping duplicate compound with CID {cid}")
                    continue

                # Clean up SMILES to remove invalid characters or trailing symbols
                smiles = re.sub(r'[.]+$', '', smiles)

                try:
                    compound_name = f"Compound_{cid}"
                    mol = Chem.MolFromSmiles(smiles)
                    if mol:
                        image_path = os.path.join(output_folder, f"{compound_name}.png")
                        Draw.MolToFile(mol, image_path)
                        drawn_compounds.add(cid)  # Mark this compound as drawn
                        print(f"Saved structure for {compound_name} at {image_path}")
                    else:
                        print(f"Invalid SMILES: {smiles}")
                except Exception as e:
                    print(f"Failed to draw structure for SMILES: {smiles} - {e}")

# Example usage
smiles_file_path = 'chemical_structures_cids.txt'  # Replace with your SMILES output file
output_folder = '/home/divy/Work/GCE4ALL/RDkit_structures'  # Replace with your desired output folder
draw_smiles_with_rdkit(smiles_file_path, output_folder)
