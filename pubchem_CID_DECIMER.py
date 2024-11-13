import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
import os
import re

# Function to sanitize filename
def sanitize_filename(filename):
    """Sanitizes a filename by removing or replacing invalid characters."""
    return re.sub(r'[<>:"/\\|?*]', '', filename)[:100]

# Function to fetch chemical data from PubChem
def fetch_chemical_data(cid):
    """Fetches the chemical data from PubChem given a CID, including the SMILES."""
    try:
        compound = pcp.Compound.from_cid(cid)
        smiles = compound.isomeric_smiles  # SMILES representation of the compound
        return smiles
    except Exception as e:
        print(f"Error fetching chemical data for PubChem CID {cid}: {e}")
        return None

# Function to draw and save a structure centered on the amide carbon
def draw_and_save_structure(smiles, output_name, output_dir):
    """Draws a chemical structure from a SMILES string, centering on the amide carbon."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"Invalid SMILES string for compound {output_name}.")
        return

    # Create the specified output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Sanitize the name to make it a valid filename
    sanitized_name = sanitize_filename(output_name)

    # Define the file path
    file_path = os.path.join(output_dir, f"{sanitized_name}.png")

    # Draw the molecule
    drawer = rdMolDraw2D.MolDraw2DCairo(500, 500)  # Image size 500x500 pixels
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    drawer.WriteDrawingText(file_path)

    print(f"Image saved as: {file_path}")

# Process the input file
def process_smiles_file(input_file_path, output_dir):
    """Reads the input file, fetches chemical data, and draws the structures."""
    with open(input_file_path, 'r') as input_file:
        for line in input_file:
            if ':' in line:
                try:
                    # Extract name and CID
                    name, cid_str = line.split(':')[0].strip(), line.split('CIDs:')[1].strip()
                    cid = int(cid_str)

                    # Fetch SMILES from PubChem
                    smiles = fetch_chemical_data(cid)
                    if smiles:
                        # Draw and save the structure
                        draw_and_save_structure(smiles, name, output_dir)
                except ValueError:
                    print(f"Failed to parse line: {line}")

# Example usage
input_file_path = 'chemical_structures_cids_2.txt'  # Input file containing SMILES and CIDs
output_directory = 'output_images'  # Directory to save images
process_smiles_file(input_file_path, output_directory)
