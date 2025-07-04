import subprocess
import os, sys
import typing as t
from pathlib import Path
from dataclasses import dataclass



def _shell_cmd(cmd: str) -> None:
    try:
        print(f"Running command: {cmd}")
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Error running command:\n{result.stderr}", file=sys.stderr)
            sys.exit(result.returncode)
        return result
    except Exception as e:
        print(f"An error occurred while running the command: {e}", file=sys.stderr)
        sys.exit(1)



@dataclass
class CrisprDetection:
    """
    Class to handle prediction of CRIPSR array.

    Attributes:
        id (str): Unique identifier for the CRISPR detection instance.
    
    Methods:
        run_crispridentify(fasta_file: Path, output_dir: Path) -> t.Optional[Path]:
            Runs CRISPRidentify on the provided FASTA file and saves the results in the specified output directory.
    
    Args:
        fasta_file (Path): The path to the FASTA file containing the genome.
        output_dir (Path): The directory where the CRISPRidentify results will be saved.
    """


    @staticmethod
    def strand_to_str(strand: str) -> str:
        if strand == "Forward":
            return "Forward"
        elif strand == "Reversed":
            return "Reverse"
        else:
            return "?"


    @classmethod
    def run_crispridentify(cls, fasta_file: Path, output_dir: Path):
        """
        Runs CRISPRidentify on the provided FASTA file and saves the results in the specified output directory.
        
        Args:
            fasta_file (Path): The path to the FASTA file containing the genome.
            output_dir (Path): The directory where the CRISPRidentify results will be saved.
        Returns:
            Path: Path to the CRISPRidentify results file if successful, None otherwise.
        """

        bioaccession = os.path.basename(fasta_file).split('.')[0]  

        # Ensure the output directory exists
        output_folder = output_dir / f"{bioaccession}_crispridentify"


        # Run CRISPRidentify command
        _shell_cmd(
        f"python3 /home/unimelb.edu.au/rbengtsson/work/spacer_analysis/tools/CRISPRidentify-1.2.1/CRISPRidentify.py "
        f"--cas True "
        f"--fasta_report True "
        f"--file {fasta_file} "
        f"--result_folder {output_folder}"
        )


        





    # @staticmethod
    # def run_casfinder(fasta_file: Path, output_dir: Path) -> t.Optional[Path]:
    #     """
    #     Runs CasFinder on the provided FASTA file and saves the results in the specified output directory.

    #     Args:
    #         fasta_file (Path): The path to the FASTA file containing the genome.
    #         output_dir (Path): The directory where the CasFinder results will be saved.

    #     Returns:
    #         CasFinder results file path (Path) if successful, None otherwise.
    #     """    
    #     try:
    #         # Ensure the output directory exists
    #         output_dir.mkdir(parents=True, exist_ok=True)