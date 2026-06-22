import os, sys
import subprocess
from pathlib import Path
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import typing as t
from typing import Optional
from dataclasses import dataclass
from concurrent.futures import ProcessPoolExecutor
import modules.analysis_utils as analysis_utils
import modules.database_utils as db_utils


def generate_sequence_file(sequence_id: str, sequence: str, output_dir: Path) -> Path:
    """
    Format the spacer sequence and write it to a file.
    """
    output_file = output_dir / f"{sequence_id}.fasta"

    with open(output_file, 'w') as f:
        SeqIO.write(
            SeqIO.SeqRecord(
                Seq(sequence),
                id=sequence_id,
                description=""),f,"fasta")

    return output_file


def blast_search(blast_db: Path, spacerSeq_file: Path, sequence_id: str, output_dir: Path, pident: float, cov_perc: float) -> Optional[Path]:
    """
    Perform a BLAST search for the given spacer sequence.
    """
    output_file = output_dir / f"{sequence_id}_blast_results.txt"

    blast_cmd = [
        "blastn",
        "-task", "blastn-short",
        "-query", f"{spacerSeq_file}",
        "-db", f"{blast_db}",
        "-max_target_seqs", "10",
        "-perc_identity", f"{pident}",
        "-qcov_hsp_perc", f"{cov_perc}",
        "-outfmt", "6",
        "-out", f"{output_file}"
    ]

    try:
        print(f"Running BLAST command: {' '.join(blast_cmd)}")
        subprocess.run(blast_cmd, check=True)
        print("BLAST search completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error during BLAST search: {e}")
    
    # clean up the fasta file after BLAST search
    rm_cmd = ["rm", spacerSeq_file]
    subprocess.run(rm_cmd)

    if os.path.getsize(output_file) > 0:
        return output_file
    else:
        rm_cmd = ["rm", output_file]
        subprocess.run(rm_cmd)
        print(f"No results found for {sequence_id}. BLAST output file removed.")
        return None



def process_sequence(blast_db, row_data, pident: float = 90, cov_perc: float = 95):
    row, output_dir = row_data
    # if len(row['sequence']) == 30:
    sequence_id = row['sequence_id']
    sequence = row['sequence']
    fasta_file = generate_sequence_file(sequence_id, sequence, output_dir)
    blast_search(blast_db, fasta_file, sequence_id, output_dir, pident=pident, cov_perc=cov_perc)


@dataclass
class PhageBlast:
    """
    Class to handle BLAST searches for phage sequences.
    """
    sequence_id: Optional[str] = None
    spacerhost: Optional[str] = None
    phage_id: Optional[str] = None
    pident: Optional[float] = None
    # spacerlength: Optional[int] = None
    algn_length: Optional[int] = None
    mismatch: Optional[int] = None
    gapopen: Optional[int] = None
    phagestart: Optional[int] = None
    phageend: Optional[int] = None
    phagelength: Optional[int] = None
    evalue: Optional[float] = None
    bitscore: Optional[float] = None
    phagecg: Optional[float] = None
    phagetaxonomy: Optional[str] = None
    completeness: Optional[str] = None
    phagehost: Optional[str] = None
    lifestyle: Optional[str] = None
    phagesource: Optional[str] = None

    @classmethod
    def retrieve_phage_info(cls, phage_id: str) -> t.Dict[str, t.Any]:
        """
        Retrieve phage information from the metadata file.
        """
        phage_metadata_df = pd.read_csv('/phagescope_db/phagescope/phagescope_phage_meta_data.tsv', sep='\t', header=0, dtype=str)
        
        result = phage_metadata_df[phage_metadata_df['Phage_ID'] == phage_id]
        return result.iloc[0].to_dict() if not result.empty else {}
        


    @classmethod
    def from_blast_line(cls, line: str) -> 'PhageBlast':
        """
        Create a PhageBlast instance from a BLAST output line.
        """
        fields = line.strip().split('\t')
        phage_instance = PhageBlast(
            sequence_id=fields[0],
            phage_id=fields[1],
            pident=float(fields[2]),
            algn_length=int(fields[3]),
            mismatch=int(fields[4]),
            gapopen=int(fields[5]),
            phagestart=int(fields[8]),
            phageend=int(fields[9]),
            phagelength=int(cls.retrieve_phage_info(fields[1]).get('Length')),
            evalue=float(fields[10]),
            bitscore=float(fields[11]),
            phagecg=float(cls.retrieve_phage_info(fields[1]).get('GC_content')),
            phagetaxonomy=cls.retrieve_phage_info(fields[1]).get('Taxonomy'),
            completeness=cls.retrieve_phage_info(fields[1]).get('Completeness'),
            phagehost=cls.retrieve_phage_info(fields[1]).get('Host'),
            lifestyle=cls.retrieve_phage_info(fields[1]).get('Lifestyle'),
            phagesource=cls.retrieve_phage_info(fields[1]).get('Phage_source')
        )

        return phage_instance
