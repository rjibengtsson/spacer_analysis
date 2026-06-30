import os, sys
from datetime import datetime
import subprocess
import typing as t
from typing import Optional
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import sys, os
import uuid
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
import modules.blast_utils as blast_utils
from modules.blast_utils import PhageBlast
from modules.analysis_utils import PhageElement




def main():
    """
    Main function to execute the BLAST search.
    """
    home_dir = Path('/home/unimelb.edu.au/rbengtsson/work/spacer_analysis/')
    array_df = db_utils.read_table_from_db(f"{home_dir}/prelim/analysis/preliminary_crisprs.csv")
    spacer_df = db_utils.read_table_from_db(f"{home_dir}/prelim/analysis/all_spacers.csv")
    output_dir = Path(f"/spacers_db/blastn_results")
    output_dir.mkdir(parents=True, exist_ok=True)

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

    seq_with_species_spacers = seq_with_species[seq_with_species['type'] == 'spacer']

    row_data_list = [(row, output_dir) for _, row in seq_with_species_spacers.iterrows()]

    blast_db = Path(f'/phagescope_db/phagescope/phagescope_db')

    pident = 85
    cov_perc = 95

    # Run blastn in parallel
    with ProcessPoolExecutor(max_workers=4) as executor:
        executor.map(blast_utils.process_sequence, [blast_db]*len(row_data_list), row_data_list, [pident]*len(row_data_list), [cov_perc]*len(row_data_list))


    df = pd.DataFrame(columns=[
        'sequence_id', 'spacerhost', 'phage_id', 'pident', 'spacerlength', 'algn_length',
        'mismatch', 'gapopen', 'phagestart', 'phageend', 'phagelength','evalue', 'bitscore',
        'phagecg', 'phagetaxonomy', 'phagecompleteness', 'phagehost',
        'phagelifestyle', 'phagesource'
    ])

    for index, row in seq_with_species_spacers.iterrows():
        sequence_id = row['sequence_id']
        speciesname = row['speciesname']
        blast_file = output_dir / f"{sequence_id}_blast_results.txt"
        if not blast_file.exists():
            pass
        else:
            for line in open(blast_file, 'r'):
                phage_instance = PhageBlast.from_blast_line(line)
                phage_instance.spacerhost = speciesname
                phage_instance.spacerlength = len(row['sequence'])
                new_df = pd.DataFrame([phage_instance.__dict__])
                df = pd.concat([df, new_df], ignore_index=True)

    timestamp = datetime.now().strftime('%d-%m-%Y-%H-%M-%S')
    df.to_csv(f"{output_dir}/spacer_phageblastn_results_{pident}identity_{cov_perc}cov_{timestamp}.csv")

if __name__ == "__main__":
    main()