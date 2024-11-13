This Bash Script should install and compile the DECIMER system as well as the DECIMER image segmentation system to work on your local system. 
It will also install an script that can be called upon to extract out the SMILES of structures present on an PDF file. 

## INSTALLATION: 
    Execute the script to install all necessary components and run the Python script:
     `   ./Install.sh    `


## After running the installation, you can use the generated Python script (decimer_smiles_extraction.py) with specific parameters:

    ./decimer_smiles_extraction.py --pdf_path <path_to_pdf> --output_folder <output_folder> --output_txt_path <output_txt_file>

    ### Replace <path_to_pdf>, <output_folder>, and <output_txt_file> with the appropriate paths for your PDF, the folder to store images, and the output file to store SMILES predictions, respectively.

#### Example: 
    ./decimer_smiles_extraction.py --pdf_path /path/to/myfile.pdf --output_folder /path/to/output_images --output_txt_path /path/to/output_smiles.txt
