"""
This module contains functions to extract data from the CrisprCas Atlas SQL database.
"""

import json
import uuid
import os, sys
import typing as t
from pathlib import Path
import psycopg2
import pandas as pd
import typing as t
from typing import Optional
import subprocess
from sqlalchemy import create_engine, text
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from dataclasses import dataclass



def get_db_credentials() -> t.Dict[str, str]:
    """
    Retrieves database credentials from a JSON file.
    """
    credentials_path = "/home/unimelb.edu.au/rbengtsson/work/spacer_analysis/modules/db_credentials.json"
    with open(credentials_path, "r") as f:
        credentials = json.load(f)
    
    return credentials


def extract_data_from_database(query: str, DATABASE: str = 'crisprcas_atlas') -> pd.DataFrame:
    """
    Connects to the PostgreSQL database and executes the provided SQL query.

    Parameters:
    - query: The SQL query to execute.
    - DATABASE: The name of the database to connect to.


    Returns:
    A pandas DataFrame containing the results of the query.
    """
    credentials = get_db_credentials()
    creds = credentials['db_credentials']

    try:
        # Establish a connection to the database
        connection = psycopg2.connect(
            host=creds['host'],
            database=DATABASE,
            user=creds['user'],
            password=creds['password']
        )
    
        # Execute the query and fetch the results into a DataFrame
        df = pd.read_sql_query(query, connection)
        
        return df

    except Exception as e:
        print(f"An error occurred while connecting to the database or executing the query: {e}")
        return pd.DataFrame()  # Return an empty DataFrame in case of error

    finally:
        if 'connection' in locals():
            connection.close()



def get_terminase_features(phage_ids: tuple, output_dir: Path) -> Path:
    """
    Fetch terminase features for a list of phage IDs.

    Args:
        phage_ids: list of phage IDs to query
        output_dir: directory to save the output file

    Returns:
        Path to the output file containing terminase features
    """
    query = text("""
        SELECT fp.phageid, f.featureid, f.start, f.end, f.strand, f.product 
        FROM feature f
        JOIN featureid_phageid fp ON f.featureid = fp.featureid
        WHERE fp.phageid IN :phage_ids
        AND (
            f.product ILIKE '%terminase%'
            OR f.product ~* '\mterl\M'
        )
    """)

    database = "phagescope"
    credentials = get_db_credentials()
    creds = credentials['db_credentials']
    
    engine = create_engine(
        f"postgresql+psycopg2://{creds['user']}:{creds['password']}"
        f"@{creds['host']}:{creds['port']}/{database}"
    )

    with engine.connect() as conn:
        result = pd.read_sql(query, conn, params={"phage_ids": phage_ids})

    output_file = output_dir / "terminase_features.csv"
    result.to_csv(output_file, index=False)
    
    return output_file




def blast_search_batch(chunk_df: pd.DataFrame, blast_db: Path, output_dir: Path,
                       pident: float, cov_perc: float, num_threads: int = 4) -> Optional[Path]:
    """
    Perform a BLAST search for a batch (chunk) of spacer sequences.
    Writes all sequences in the chunk to a single multi-FASTA file,
    runs blastn with -num_threads, and returns the output path.
    """
    uuid_str = str(uuid.uuid4())
    fasta_output = output_dir / f"batch_{uuid_str}.fasta"
    blastn_output = output_dir / f"batch_{uuid_str}_blastn.txt"

    spacer_idx = chunk_df.groupby("operon_id").cumcount()
    records = [
        SeqRecord(Seq(row.spacer_sequence), id=f"{row.operon_id}_sp{spacer_idx.loc[row.Index]}", description="")
        for row in chunk_df.itertuples()
    ]

    with open(fasta_output, 'w') as f:
        SeqIO.write(records, f, "fasta")

    blast_cmd = [
        "blastn", "-task", "blastn-short",
        "-query", str(fasta_output),
        "-db", str(blast_db),
        "-num_threads", str(num_threads),
        "-max_target_seqs", "10",
        "-perc_identity", str(pident),
        "-qcov_hsp_perc", str(cov_perc),
        "-outfmt", "6",
        "-out", str(blastn_output)
    ]

    try:
        print(f"Running BLAST for batch {uuid_str} ({len(chunk_df)} spacers)")
        subprocess.run(blast_cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"Error during BLAST search for batch {uuid_str}: {e}")
        fasta_output.unlink(missing_ok=True)
        return None
    finally:
        fasta_output.unlink(missing_ok=True)

    if blastn_output.exists() and os.path.getsize(blastn_output) > 0:
        return blastn_output
    else:
        blastn_output.unlink(missing_ok=True)
        return None




@dataclass
class PhageBlast:
    """
    Class to handle BLAST searches for phage sequences.
    """
    sequence_id: Optional[str] = None
    phage_id: Optional[str] = None
    pident: Optional[float] = None
    algn_length: Optional[int] = None
    mismatch: Optional[int] = None
    gapopen: Optional[int] = None
    phagestart: Optional[int] = None
    phageend: Optional[int] = None
    phagelength: Optional[int] = None
    evalue: Optional[float] = None
    bitscore: Optional[float] = None
    phage_cg: Optional[float] = None
    phage_taxonomy: Optional[str] = None
    phage_completeness: Optional[str] = None
    phage_host: Optional[str] = None
    phage_lifestyle: Optional[str] = None
    phage_source: Optional[str] = None

    @classmethod
    def retrieve_phage_info(cls, phage_id: str) -> t.Dict[str, t.Any]:
        """
        Retrieve phage information from the metadata file.
        """
        phage_metadata_df = pd.read_csv('/phagescope_db/phagescope_db/data/phagescope/phagescope_phage_meta_data.tsv', sep='\t', header=0, dtype=str)
        
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
            phage_cg=float(cls.retrieve_phage_info(fields[1]).get('GC_content')),
            phage_taxonomy=cls.retrieve_phage_info(fields[1]).get('Taxonomy'),
            phage_completeness=cls.retrieve_phage_info(fields[1]).get('Completeness'),
            phage_host=cls.retrieve_phage_info(fields[1]).get('Host'),
            phage_lifestyle=cls.retrieve_phage_info(fields[1]).get('Lifestyle'),
            phage_source=cls.retrieve_phage_info(fields[1]).get('Phage_source')
        )

        return phage_instance



def filter_results(input_file: Path, pident_threshold: float, cov_perc_threshold: float) -> pd.DataFrame:
    """
    Filter results based on various thresholds and selection criteria.
    """

    phage_completeness = 'High-quality'



    df = pd.read_csv(input_file)

    filtered_df = df[(df['pident'] >= pident_threshold) & (df['algn_length'] / df['phagelength'] * 100 >= cov_perc_threshold)]
    return filtered_df




def reorient_coordinates(old_pos, new_start, genome_length):
    """Convert old coordinate to reoriented coordinate."""
    return (old_pos - new_start) % genome_length