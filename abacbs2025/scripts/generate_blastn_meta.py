import sys, os
from pathlib import Path
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq

# Add the directory containing your other modules
target_dir = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(target_dir))

import modules.database_utils as db_utils
from modules.blast_utils import PhageBlast
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
import json
import typing as t
import subprocess
from typing import Optional
import glob


def get_species_from_operon_id(operon_id: str) -> Optional[str]:
    cas13b_summary_file = "/home/unimelb.edu.au/rbengtsson/work/spacer_analysis/abacbs2025/data/crispr_cas_atlas_Cas13b_summary_edit.csv"
    metadata_df = pd.read_csv(cas13b_summary_file, header=0)
    match = metadata_df[metadata_df['operon_id'] == operon_id]
    if not match.empty:
        return match.iloc[0]['species']
    return None

def retrieve_operon_ids(blast_file_path: Path) -> str:
    basename = os.path.basename(blast_file_path)
    operon_id = basename.replace("_blast_results.txt", "")
    return operon_id



def process_blast_results(blast_result: str, blastn_out_dir: Path):
    operon_id = retrieve_operon_ids(Path(blast_result))

    df = pd.DataFrame(columns = ['sequence_id', 'spacerhost', 'phage_id', 'pident', 'algn_length',
       'mismatch', 'gapopen', 'phagestart', 'phageend', 'phagelength',
       'evalue', 'bitscore', 'phagecg', 'phagetaxonomy', 'completeness',
       'phagehost', 'lifestyle', 'phagesource'])

    for line in open(blast_result, 'r'):
        phage_instance = PhageBlast.from_blast_line(line)
        phage_instance.spacerhost = get_species_from_operon_id(operon_id)
        new_df = pd.DataFrame([phage_instance.__dict__])
        df = pd.concat([df, new_df], ignore_index=True)

    df.to_csv(f"{blastn_out_dir}/{operon_id}_results.tsv", sep='\t', index=False)




def main():

    blastn_out_dir = Path("/spacers_db/blastn_results/")
    file_paths = glob.glob(f"{blastn_out_dir}/*_blast_results.txt", recursive=True)

    # run in parallel
    with ProcessPoolExecutor(max_workers=4) as executor:
        executor.map(process_blast_results, file_paths, [blastn_out_dir]*len(file_paths))


if __name__ == "__main__":
    main()