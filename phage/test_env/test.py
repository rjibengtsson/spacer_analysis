import os, sys
import uuid
import re
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
    home_dir = Path('/home/unimelb.edu.au/rbengtsson/work/spacer_analysis/')
    blastn_file = f"/spacers_db/blastn_results/spacer_phageblastn_results_85identity_95cov_01-08-2025-16-46-20.csv"
    blastn_df = db_utils.read_table_from_db(blastn_file)

    
    # fasta_file = f"/phagescope_db/phagescope/phagescope_combined.fasta"

    # species = 'Prevotella pectinovora'
    # outfile = home_dir / "phage/test_env/phage_genomes.fasta"

    # phage_ids = blastn_df[blastn_df['spacerhost'] == species]['phage_id'].unique().tolist()

    # with ProcessPoolExecutor(max_workers=4) as executor:
    #     record_list = list(executor.map(
    #         analysis_utils.get_genome_seq_from_fasta, 
    #         [fasta_file]*len(phage_ids), phage_ids
    #         ))

    # with open(f"{outfile}", "w") as output_handle:
    #     SeqIO.write(record_list, output_handle, "fasta")

#     for record in SeqIO.parse(phage_genome, "fasta"):
#         if record.id == phage_id:
#             start = int(127657) - 1
#             end = int(127848)
#             # extension_length = 200
#             # protospacer_seq = record.seq[start-extension_length:end+extension_length]
#             # record = SeqRecord(seq=protospacer_seq, id=phage_id, description=f"")
#             # record_list.append(record)
#             seq = record.seq[start:end]
#             print(seq.reverse_complement())


#     # with open(f"{species.replace(' ','_')}_{sequence_id}.fasta", "w") as output_handle:
#     #     SeqIO.write(record_list, output_handle, "fasta")


if __name__ == "__main__":
    main()