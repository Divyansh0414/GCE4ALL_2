#!/bin/bash

# Installation script for Conda, DECIMER, and DECIMER Segmentation

# Step 1: Install Miniconda if not already installed
# This step checks if conda is installed on the system. If not, it downloads and installs Miniconda.
if ! command -v conda &> /dev/null
then
    echo "Conda not found. Installing Miniconda..."
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda
    export PATH="$HOME/miniconda/bin:$PATH"
    echo "Miniconda installed."
elif [ -d "$HOME/miniconda" ]; then
    echo "Miniconda installation detected. Updating Miniconda..."
    bash Miniconda3-latest-Linux-x86_64.sh -b -u -p $HOME/miniconda
    export PATH="$HOME/miniconda/bin:$PATH"
fi

# Step 2: Create and activate a conda environment for DECIMER
# This step creates a new conda environment named DECIMER_ENV and activates it.
if conda env list | grep -q "DECIMER_ENV"; then
    echo "Conda environment DECIMER_ENV already exists. Activating it..."
else
    conda create --name DECIMER_ENV python=3.10 -y
fi
source $HOME/miniconda/bin/activate DECIMER_ENV

# Step 3: Install necessary Python packages
# This step installs the necessary Python packages such as decimer, pdf2image, Pillow, numpy, and argparse.
pip install --upgrade pip
pip install decimer pdf2image Pillow numpy argparse

# Step 4: Clone and install DECIMER-Image-Segmentation
# Clone the DECIMER-Image-Segmentation repository, install it, and install poppler for PDF processing.
if [ -d "DECIMER-Image-Segmentation" ]; then
    echo "DECIMER-Image-Segmentation directory already exists. Skipping cloning..."
else
    git clone https://github.com/Kohulan/DECIMER-Image-Segmentation
fi
cd DECIMER-Image-Segmentation
pip install .
conda install -c conda-forge poppler -y
cd ..

# Step 5: Install TensorFlow and Keras compatible versions
# Uninstall any existing versions of TensorFlow and Keras, then install specific compatible versions.
if pip show tensorflow &> /dev/null; then
    echo "Uninstalling existing TensorFlow..."
    pip uninstall tensorflow -y
fi
if pip show keras &> /dev/null; then
    echo "Uninstalling existing Keras..."
    pip uninstall keras -y
fi

# Attempt to install a compatible version of TensorFlow and Keras
if ! pip install "tensorflow>=2.10,<2.13" "keras>=2.10,<2.13"; then
    echo "Compatible TensorFlow or Keras versions not found. Attempting to install the latest available version compatible with DECIMER."
    pip install tensorflow keras
fi

# Python script for extracting chemical structures and predicting SMILES
# Create a Python script that extracts chemical structures from a PDF and predicts SMILES strings.
cat <<EOF > decimer_smiles_extraction.py
#!/usr/bin/env python

"""
This script extracts chemical structures from a given PDF file and predicts their SMILES strings.
It uses DECIMER and DECIMER Segmentation to perform these tasks.

Steps:
1. Convert each page of the PDF into images.
2. Segment chemical structures from the images.
3. Predict SMILES strings for each segmented structure.
"""

from pdf2image import convert_from_path
from PIL import Image
import os
import io
import numpy as np
from DECIMER import predict_SMILES
from decimer_segmentation import segment_chemical_structures_from_file
import argparse

def extract_images_from_pdf(pdf_path, output_folder):
    """
    Convert each page of a PDF file into individual image files.
    
    Args:
        pdf_path (str): Path to the input PDF file.
        output_folder (str): Folder where the extracted images will be saved.
    
    Returns:
        list: A list of paths to the extracted image files.
    """
    # Ensure the output folder exists
    os.makedirs(output_folder, exist_ok=True)
    
    # Convert each page of the PDF to images
    images = convert_from_path(pdf_path)
    
    # Save each page as an individual image file
    image_paths = []
    for i, image in enumerate(images):
        img_path = os.path.join(output_folder, f"page_{i + 1}.png")
        image.save(img_path, 'PNG')
        image_paths.append(img_path)
    
    return image_paths

def convert_images_to_smiles(image_paths, output_txt_path):
    """
    Convert segmented chemical structures from images into SMILES strings and save them to a text file.
    
    Args:
        image_paths (list): A list of paths to the extracted image files.
        output_txt_path (str): Path to the output text file where SMILES strings will be saved.
    """
    with open(output_txt_path, 'w') as output_file:
        for i, img_path in enumerate(image_paths):
            try:
                # Segment the structures on the page using DECIMER segmentation
                segmented_structures = segment_chemical_structures_from_file(img_path, expand=True)
                
                for j, structure_img in enumerate(segmented_structures):
                    # Convert numpy array to PIL image
                    if isinstance(structure_img, np.ndarray):
                        structure_img = Image.fromarray(structure_img)
                    
                    # Convert each structure to a format compatible with DECIMER
                    img_byte_arr = io.BytesIO()
                    structure_img.save(img_byte_arr, format='PNG')
                    img_byte_arr.seek(0)

                    # Predict SMILES using DECIMER
                    smiles = predict_SMILES(img_byte_arr)
                    
                    # Write the identifier and SMILES to the file
                    output_file.write(f"Page_{i + 1}_Structure_{j + 1}: {smiles}\n")
                    
                    print(f"Converted Page {i + 1} Structure {j + 1} to SMILES: {smiles}")
            except Exception as e:
                print(f"Failed to process structures on page {i + 1}: {e}")
                output_file.write(f"Page_{i + 1}: Failed to convert structures\n")

def main():
    """
    Main function to parse arguments and run the extraction and conversion processes.
    """
    parser = argparse.ArgumentParser(description="Extract chemical structures from a PDF and predict SMILES strings.")
    parser.add_argument('--pdf_path', type=str, required=True, help='Path to the input PDF file.')
    parser.add_argument('--output_folder', type=str, required=True, help='Folder to store extracted images.')
    parser.add_argument('--output_txt_path', type=str, required=True, help='Path to the output text file with predicted SMILES strings.')
    
    args = parser.parse_args()
    
    # Step 1: Extract images from the original PDF
    image_paths = extract_images_from_pdf(args.pdf_path, args.output_folder)
    
    # Step 2: Convert each segmented structure to SMILES and save to a text file
    convert_images_to_smiles(image_paths, args.output_txt_path)

if __name__ == "__main__":
    main()
EOF

chmod +x decimer_smiles_extraction.py

# Provide usage instructions
echo "Installation complete. You can now use the script as follows:"
echo "./decimer_smiles_extraction.py --pdf_path <path_to_pdf> --output_folder <output_folder> --output_txt_path <output_txt_file>"