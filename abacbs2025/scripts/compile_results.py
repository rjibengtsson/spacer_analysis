import sys, os
from pathlib import Path
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq

# Add the directory containing your other modules
target_dir = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(target_dir))

from abacbs2025.scripts.generate_blastn_meta import process_blast_results
import modules.database_utils as db_utils
from modules.blast_utils import PhageBlast
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
import json
import typing as t
import subprocess
from typing import Optional
import glob


def remove_redundant_hits(df: pd.DataFrame) -> pd.DataFrame:
    """
    Remove redundant BLAST hits, keeping only the best hit per phage_id.
    """
    
    filtered_df = pd.DataFrame(columns=df.columns)

    # Step 1: keep only duplicated groups
    dup_groups = df[df.duplicated('sequence_id', keep=False)]

    # Step 2: Keep rows where 'bitscore' equals the group max
    max_rows = dup_groups[dup_groups['bitscore'] == dup_groups.groupby('sequence_id')['bitscore'].transform('max')]

    # Step 3: If multiple have same 'col_check', keep one (drop duplicates)
    result = max_rows.drop_duplicates(subset=['sequence_id', 'phagelength'], keep='first')
    
    # Step 4: Append non-duplicated rows from original df
    non_dup_rows = df[~df['sequence_id'].isin(dup_groups['sequence_id'])]
    result = pd.concat([result, non_dup_rows], ignore_index=True) 

    return result



def main():

    blastn_out_dir = Path("/spacers_db/blastn_results/")
    file_paths = glob.glob(f"{blastn_out_dir}/*_results.tsv", recursive=True)

    df_header = pd.read_csv(file_paths[0], header=0, sep='\t').columns.tolist()
    df_final = pd.DataFrame(columns=df_header)

    for file_path in file_paths:
        df = pd.read_csv(file_path, header=0, sep='\t')
        filtered_df = remove_redundant_hits(df)
        df_final = pd.concat([df_final, filtered_df], ignore_index=True)
        
    df_final.to_csv(f"/home/unimelb.edu.au/rbengtsson/work/spacer_analysis/abacbs2025/data/spacer_phageblastn_results_90id_95cov_23-10-2025.csv", index=False)

    




if __name__ == "__main__":
    main()