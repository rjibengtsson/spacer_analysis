import subprocess
import os
from pathlib import Path
import pandas as pd
import typing as t
from typing import Optional
from Bio import SeqIO
import gffutils
from dataclasses import dataclass


def get_fasta_from_ncbi(ftppath_genbank: str, biosampleaccn, output_dir: Path) -> t.Optional[Path]:
    """
    Fetches the FASTA file for a given ftppath_genbank from NCBI assembly database.

    Args:
        ftppath_genbank (str): The genbank path to fetch fasta file from.
        biosampleaccn (str): The biosample accession number.
        output_dir (Path): The directory where the FASTA file will be saved.

    Returns:
        fasta file path (Path) if successful, None otherwise.
    """
    try:
        fasta_file_path = f"{ftppath_genbank}/{ftppath_genbank.split('/')[-1]}_genomic.fna.gz"
        fasta_file_path  = fasta_file_path.replace("GCA", "GCF")  # Convert Genbank to RefSeq

        # Construct the wget command to fetch the fasta file
        command = [
            "wget",
            "-A", "*.fna.gz",  # accept only .fna.gz files
            "-P", str(output_dir),  # specify output directory
            fasta_file_path
        ]

        # Execute the command
        subprocess.run(command, check=True)

        # change the file name to biosample accession number
        # This is to ensure that the file name is consistent with the biosample accession number
        # and can be easily identified later.
        new_file_path = output_dir / f"{biosampleaccn}.fna.gz"
        os.rename(output_dir / os.path.basename(fasta_file_path), new_file_path)

        # Return the path to the downloaded fasta file
        return new_file_path

    except subprocess.CalledProcessError as e:
        print(f"Error fetching fasta file: {e}")
        return None
    

def get_gff_from_ncbi(ftppath_genbank: str, biosampleaccn: str, output_dir: Path) -> t.Optional[Path]:
    """
    Fetches the GFF file for a given ftppath_genbank from NCBI assembly database.
    
    Args:
        ftppath_genbank (str): The genbank path to fetch gff file from.
        biosampleaccn (str): The biosample accession number.
        output_dir (Path): The directory where the GFF file will be saved.
    
    Returns:
        gff file path (Path) if successful, None otherwise.
    """
    try:
        gff_file_path = f"{ftppath_genbank}/{ftppath_genbank.split('/')[-1]}_genomic.gff.gz"
        gff_file_path = gff_file_path.replace("GCA", "GCF")  # Convert Genbank to RefSeq

        # Construct the wget command to fetch the gff file
        command = [
            "wget",
            "-A", "*.gff.gz",  # accept only .gff.gz files
            "-P", str(output_dir),  # specify output directory
            gff_file_path
        ]

        # Execute the command
        subprocess.run(command, check=True)

        # change the file name to biosample accession number
        # This is to ensure that the file name is consistent with the biosample accession number
        # and can be easily identified later.
        new_file_path = output_dir / f"{biosampleaccn}.gff.gz"
        os.rename(output_dir / os.path.basename(gff_file_path), new_file_path)

        # Return the path to the downloaded gff file
        return new_file_path

    except subprocess.CalledProcessError as e:
        print(f"Error fetching gff file: {e}")
        return None



def get_gbff_from_ncbi(ftppath_genbank: str, biosampleaccn: str, output_dir: Path) -> t.Optional[Path]:
    """
    Fetches the GenBank flat file (GBFF) for a given ftppath_genbank from NCBI assembly database.
    
    Args:
        ftppath_genbank (str): The genbank path to fetch gbff file from.
        biosampleaccn (str): The biosample accession number.
        output_dir (Path): The directory where the GBFF file will be saved.
    
    Returns:
        gbff file path (Path) if successful, None otherwise.
    """
    try:
        gbff_file_path = f"{ftppath_genbank}/{ftppath_genbank.split('/')[-1]}_genomic.gbff.gz"
        gbff_file_path = gbff_file_path.replace("GCA", "GCF")  # Convert Genbank to RefSeq

        # Construct the wget command to fetch the gbff file
        command = [
            "wget",
            "-A", "*.gbff.gz",  # accept only .gbff.gz files
            "-P", str(output_dir),  # specify output directory
            gbff_file_path
        ]

        # Execute the command
        subprocess.run(command, check=True)

        # change the file name to biosample accession number
        # This is to ensure that the file name is consistent with the biosample accession number
        # and can be easily identified later.
        new_file_path = output_dir / f"{biosampleaccn}.gbff.gz"
        os.rename(output_dir / os.path.basename(gbff_file_path), new_file_path)

        # Return the path to the downloaded gbff file
        return new_file_path

    except subprocess.CalledProcessError as e:
        print(f"Error fetching gbff file: {e}")
        return None


def get_fasta_from_allthebacteria(biosampleaccn: str, output_dir: Path) -> t.Optional[Path]:
    """
    Fetches the FASTA file for a given sample_accession from AWS.
    
    Args:
        biosampleaccn (str): Biosample accession to fetch fasta file from.
        output_dir (Path): The directory where the FASTA file will be saved.

    Returns:
        fasta file path (Path) if successful, None otherwise.
    """
    try:
        # Ensure the output directory exists
        output_dir.mkdir(parents=True, exist_ok=True)

        S3_url = f"https://allthebacteria-assemblies.s3.eu-west-2.amazonaws.com/{biosampleaccn}.fa.gz"
        
        # Construct the wget command to fetch the fasta file
        command = [
            "wget",
            "-P", str(output_dir),  # specify output directory
            S3_url
        ]

        # Execute the command
        subprocess.run(command, check=True)

        # change prefix of file name
        new_file_path = output_dir / f"{biosampleaccn}.fna.gz"
        os.rename(output_dir / f"{biosampleaccn}.fa.gz", new_file_path)
        
        # Return the path to the downloaded fasta file
        return new_file_path
    
    except subprocess.CalledProcessError as e:
        print(f"Error fetching fasta file from allthebacteria: {e}")
        return None
    

    

def get_contig(fasta_file: Path, contig_id: str, output_dir: Path) -> Path:
    """
    Extracts a specific contig from the FASTA file.
    
    Args:
        fasta_file (Path): The path to the FASTA file.
        contig_id (str): The ID of the contig to extract.
        output_dir (Path): The directory where the contig file will be saved.
    
    Returns:
        Path: The path to the extracted contig file.
    """
    try:
        # Ensure the output directory exists
        output_dir = fasta_file.parent / "contigs"
        output_dir.mkdir(parents=True, exist_ok=True)

        contig_file_path = output_dir / f"{contig_id}.fna"

        # Extract the contig using SeqIO
        with open(fasta_file, "r") as fasta_handle, open(contig_file_path, "w") as contig_handle:
            for record in SeqIO.parse(fasta_handle, "fasta"):
                if record.id == contig_id:
                    SeqIO.write(record, contig_handle, "fasta")
                    break

        return contig_file_path
    except Exception as e:
        print(f"Error extracting contig: {e}")
        return None




@dataclass
class CasContig:
    """
    Class to handle the Cas contig data.
    
    Attributes:
        contig_id (str): The ID of the contig.
        contig_seq (str): The sequence of the contig.
        cas_type (str): The type of the Cas system.
        orientation (str): The orientation of the Cas gene.
        start (int): The start position of the Cas gene.
        end (int): The end position of the Cas gene.
    """

    contig_id: Optional[str] = None
    cas_type: Optional[str] = None
    locus_tag: Optional[str] = None
    orientation: Optional[str] = None
    start: Optional[int] = None
    end: Optional[int] = None


    @staticmethod
    def strand_to_str(strand: str) -> str:
        if strand == "+":
            return "Forward"
        elif strand == "-":
            return "Reverse"
        else:
            return "?"


    @classmethod
    def get_cas_info(cls, gbff_file: Path, gff_file: Path, cas_type="cas13b") -> t.Optional['CasContig']:
        """
        Retrieve Cas information from GenBank and GFF file.
        
        Args:
            cas_type (str): The type of the Cas system to search for.
            gbff_file (Path): Path to the GenBank file.
            gff_file (Path): Path to the GFF file.
        
        Returns:
            CasContig var: The class CasContig variable.
        """

        # Parse GenBank file, set the locus_tag and cas_type
        for record in SeqIO.parse(gbff_file, "genbank"):
            for f in record.features:
                if f.type == "CDS" and (f.qualifiers.get("gene") == [cas_type]):
                    # Create CasContig object
                    cls = CasContig(
                        cas_type=f.qualifiers.get("gene", [""])[0],
                        locus_tag= f.qualifiers.get("locus_tag", [""])[0])
          
        # If locus_tag is not found, return None
        if not cls.locus_tag:
            return None

        db = gffutils.create_db(
                str(gff_file),
                dbfn='mydb.db',
                force=True,
                keep_order=True,
                merge_strategy='merge')

        db = gffutils.FeatureDB('mydb.db')

        # Search all features with matching locus_tag
        for feature in db.all_features():
            if feature.featuretype == "CDS" and feature.attributes.get("locus_tag", [None])[0] == cls.locus_tag:
                # Create CasContig object
                return CasContig(
                    cas_type=cls.cas_type,
                    locus_tag=cls.locus_tag,
                    contig_id=feature.chrom,
                    orientation=cls.strand_to_str(feature.strand),
                    start=feature.start,
                    end=feature.end)
                            
        return None



        

    @classmethod
    def get_cas_contig(cls, contig_id: str, fasta_file: Path) -> Path:
        """
        Extracts the contig sequence from the FASTA file based on the contig ID.
        
        Args:
            fasta_file (Path): The path to the FASTA file.
        
        Returns:
            contig fasta file path (Path) if successful, None otherwise.
        """
        bioaccession = os.path.basename(fasta_file).split('.')[0]
        output_dir = fasta_file.parent      
        contig_file = output_dir / f"{bioaccession}_{contig_id}.fna"

        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.id == contig_id:
                # Write the contig sequence to a new FASTA file
                SeqIO.write(record, contig_file, "fasta")
                return contig_file
    

