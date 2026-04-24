import sys, os
from pathlib import Path
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq

# Add the directory containing your other modules
target_dir = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(target_dir))

import modules.database_utils as db_utils
import modules.blast_utils as blast_utils
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
import json
import typing as t
import subprocess
from typing import Optional



def retrieve_spacers(json_file: str, operon_id: str) -> t.List[str]:
    with open(json_file, "r") as f:
        records = [json.loads(line) for line in f]
        for rec in records:
            if rec.get("operon_id") == operon_id:
                spacers = rec.get("crispr", {})[0].get("crispr_spacers", [])
                cleaned = list(filter(None, spacers))
    return cleaned

def filter_spacers_by_length(spacers: t.List[str]) -> t.List[str]:
    return [s for s in spacers if len(s) == 30]



def generate_sequence_file(sequence_id: str, sequence_list: list, output_dir: Path) -> Path:
    """
    Format the spacer sequence and write it to a file.
    """
    output_file = output_dir / f"{sequence_id}.fasta"

    counter = 1
    with open(output_file, 'w') as f:
        for sequence in sequence_list:
            SeqIO.write(
                SeqIO.SeqRecord(
                    Seq(sequence),
                    id=f"{sequence_id}_{counter}",
                    description=""), f, "fasta")
            counter += 1

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


def blastn_sequence(blast_db, 
                     spacer_sequences: list, 
                     operon_id: str,
                     output_dir: Path,
                     pident: float = 90, 
                     cov_perc: float = 95):
    fasta_file = generate_sequence_file(operon_id, spacer_sequences, output_dir)
    blast_search(blast_db, fasta_file, operon_id, output_dir, pident=pident, cov_perc=cov_perc)


def process_sequence(blast_db, atlas_json_file: str, operon_id: str):
    spacers = retrieve_spacers(atlas_json_file, operon_id)
    filtered_spacers = filter_spacers_by_length(spacers)
    if filtered_spacers:
        blastn_sequence(blast_db, 
                         spacer_sequences=filtered_spacers, 
                         operon_id=operon_id,
                         output_dir=Path("/spacers_db/blastn_results/"),
                         pident=90, 
                         cov_perc=95)
    else:
        pass



def main():
    atlas_json_file = "/spacers_db/crispr-cas-atlas/crispr-cas-atlas-v1.0-VI.jsonl"
    cas13b_summary_file = "/home/unimelb.edu.au/rbengtsson/work/spacer_analysis/abacbs2025/data/crispr_cas_atlas_Cas13b_summary_edit.csv"
    blast_db = Path('/phagescope_db/phagescope/phagescope_db')

    df = pd.read_csv(cas13b_summary_file, header=0)
    operon_id_list = df['operon_id'].tolist()
    # Run blastn in parallel
    with ProcessPoolExecutor(max_workers=4) as executor:
        executor.map(process_sequence, [blast_db]*len(operon_id_list), [atlas_json_file]*len(operon_id_list), operon_id_list)
    




if __name__ == "__main__":
    main()