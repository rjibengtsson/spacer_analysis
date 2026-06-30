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


def main():
    home_dir = Path('/home/unimelb.edu.au/rbengtsson/work/spacer_analysis/')
    array_df = db_utils.read_table_from_db(f"{home_dir}/prelim/analysis/preliminary_crisprs.csv")
    spacer_df = db_utils.read_table_from_db(f"{home_dir}/prelim/analysis/all_spacers.csv")
    output_dir = Path(f"{home_dir}/phage/test_env/")

    array_filt = array_df.drop_duplicates(subset=['array_id'])
    array_filt = array_filt[
        (array_filt['category'] != 'Alternative') &
        (array_filt['spacer_len_var'] < 1) &
        (array_filt['filter'] == 'PASSED')
    ]

    filt_seqs = spacer_df[spacer_df['array_id'].isin(array_filt['array_id'])]

    seq_with_species = filt_seqs.merge(
        array_df[['array_id', 'speciesname', 'biosampleaccn']],
        on='array_id',
        how='left'        
    )

    species = 'Prevotella pectinovora'
    # sequence_id = '572d304e-117a-467b-90ac-2c3195ca24b1'

    seq_with_species_spacers = seq_with_species[(seq_with_species['type'] == 'spacer') &
                                                (seq_with_species['speciesname'] == species)]
    
    sequence_id_list = seq_with_species_spacers['sequence_id'].unique().tolist()

    blastn_file = f"/spacers_db/blastn_results/spacer_phageblastn_results_85identity_95cov_01-08-2025-16-46-20.csv"
    
    out_file = f"{output_dir}/protospacer_phage_intersection_{species}.csv"
    with ProcessPoolExecutor(max_workers=4) as executor:
        df_list = list(executor.map(
            PhageElement.get_intersection, 
            sequence_id_list, 
            [blastn_file]*len(sequence_id_list), 
            [output_dir]*len(sequence_id_list)
            ))
    
    combined_df = pd.concat(df_list, ignore_index=True)
    combined_df.to_csv(out_file, index=False)
        



if __name__ == "__main__":
    main()