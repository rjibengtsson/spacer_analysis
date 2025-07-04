import subprocess
import os
from pathlib import Path
import typing as t
from typing import Optional
from Bio import SeqIO
import gffutils
from dataclasses import dataclass
import yaml



def get_storage_path() -> Path:
    """
    This function is used to define the storage path for genome files.
    Returns:
        Path: Path object pointing to the storage location.
    """
    
    storage_path = Path("/spacers_db/genomes")
    
    return storage_path


def generate_species_yaml(output_dir: Path, biosampleaccn: str, species: str) -> Path:
    """
    Generates a YAML file for the species with the given name.
    
    Args:
        species (str): The species name to be included in the YAML file.
    
    Returns:
        species YAML file path (Path) if successful, None otherwise.
    """
    try:
        # Define the path for the species YAML file


        species_yaml_file =  output_dir / f"{biosampleaccn}_submol.yaml"

        yaml_content = {
            'organism': {
                'genus_species': species,
                }
            }
        # Write the YAML content to the file
        with open(species_yaml_file, 'w') as yaml_file:
            yaml.dump(yaml_content, yaml_file, default_flow_style=False)

        return str(species_yaml_file)

    except Exception as e:
        print(f"Error generating species YAML: {e}")
        return None



def generate_yaml_config(fasta_file: Path, species: str) -> Path:
    """
    Generates a YAML configuration file for PGAP (Prokaryotic Genome Annotation Pipeline).
    
    Args:
        output_dir (Path): The directory where the YAML file will be saved.
        biosampleaccn (str): The biosample accession number.
        species (str): The species name for the genome.
    
    Returns:
        yaml file path (Path) if successful, None otherwise.
    """
    try:
        # Check if the fasta file exists
        file = Path(fasta_file)
        output_dir = file.parent

        if not Path(fasta_file).exists():
            print(f"FASTA file {fasta_file} does not exist.")
            return None

        bioaccession = os.path.basename(fasta_file).split('.')[0]       

        # Generate species YAML file
        species_yaml_file = generate_species_yaml(output_dir, bioaccession, species)
        if species_yaml_file is None:
            print("Failed to generate species YAML file.")
            return None

        yaml_content = {
            'fasta': {
                'class': 'File',
                'location': fasta_file
            },
            'submol': {
                'class': 'File',
                'location': species_yaml_file
                }
        }
    
        # Write the YAML content to the file
        yaml_file_path = output_dir / f"{bioaccession}.yaml"
        with open(yaml_file_path, 'w') as yaml_file:
            yaml.dump(yaml_content, yaml_file, default_flow_style=False)

        return yaml_file_path

    except Exception as e:
        print(f"Error generating YAML config: {e}")
        return None



def run_pgap(fasta_file: Path, species: str) -> Path:
    """
    Runs the PGAP (Prokaryotic Genome Annotation Pipeline) on the downloaded FASTA files.
    
    Args:
        fasta_file (Path): The path to the FASTA file to be annotated.
    
    Returns:
        File path of the PGAP output directory.
    """
    try:
        bioaccession = os.path.basename(fasta_file).split('.')[0]

        file = Path(fasta_file)
        output_dir_pgap = file.parent / f"{bioaccession}_pgap_output"

        # Generate the YAML configuration file
        yaml_file = generate_yaml_config(fasta_file, species)
        if yaml_file is None:
            print("Failed to generate YAML configuration file.")
            return None

        # Construct the PGAP command
        command = [
            "/home/unimelb.edu.au/rbengtsson/work/spacer_analysis/tools/pgap-2025-05-06.build7983/pgap.py", "-r",
            "-o", str(output_dir_pgap),
            str(yaml_file), "--taxcheck"
        ]

        print(f"Running PGAP with command: {' '.join(command)}")
        
        # Execute the command
        subprocess.run(command, check=True  )

        # return output_dir
    except subprocess.CalledProcessError as e:
        print(f"Error running PGAP: {e}")
        return None