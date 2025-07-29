import os, sys
import uuid
import subprocess
import typing as t
from typing import Optional
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import sys, os
import pandas as pd
import json
from pathlib import Path
from dataclasses import dataclass
from concurrent.futures import ProcessPoolExecutor


# Add the directory containing your other modules
target_dir = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(target_dir))

import modules.analysis_utils as analysis_utils
import modules.database_utils as db_utils
from modules.blast_utils import PhageBlast
from modules.analysis_utils import PhageElement

home_dir = Path('/home/unimelb.edu.au/rbengtsson/work/spacer_analysis/')


def main():
    """
    Main function to execute the BLAST search.
    """
    home_dir = Path('/home/unimelb.edu.au/rbengtsson/work/spacer_analysis/')
    blastn_file = f"{home_dir}/phage/spacer_phageblastn_results.csv"

    blastn_df = db_utils.read_table_from_db(blastn_file)
    sequence_id_list = blastn_df['sequence_id'].unique().tolist()

    out_dir = Path(f"/spacers_db/blastn_results/")

    # PhageElement.get_intersection(sequence_id_list, blastn_file, out_dir)
    
    # run intersection in parallel
    with ProcessPoolExecutor(max_workers=4) as executor:
        executor.map(PhageElement.get_intersection, sequence_id_list, [blastn_file]*len(sequence_id_list), [out_dir]*len(sequence_id_list))


if __name__ == "__main__":
    main()