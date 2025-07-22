import sys
from pathlib import Path
import os
import uuid
from sqlalchemy import create_engine, text
from datetime import datetime

# Add the directory containing your other modules
target_dir = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(target_dir))


import modules.database_utils as db_utils
from modules.database_utils import Allthebacteria
import pandas as pd
import modules.genome_utils as genome_utils
from modules.genome_utils import CasContig
import modules.annotation as annotation
from modules.crispr_detecion import CrisprDetection
from modules.analysis_utils import ArrayCandidate
import modules.analysis_utils as analysis_utils
from Bio import SeqIO
from Bio.Seq import Seq
import logomaker


# Display all rows
pd.set_option('display.max_rows', None)

home_dir = Path("/home/unimelb.edu.au/rbengtsson/work/spacer_analysis")
array_df = db_utils.read_table_from_db(f"{home_dir}/prelim/analysis/preliminary_crisprs.csv")
spacer_df = db_utils.read_table_from_db(f"{home_dir}/prelim/analysis/all_spacers.csv")

# filter away duplicates
array_filt = array_df.drop_duplicates(subset=['array_id'])

# Distingush between categories
array_filt = array_filt[
    (array_filt['category'] != 'Alternative') &
    (array_filt['spacer_len_var'] < 1) &
    (array_filt['filter'] == 'PASSED')
]

filt_seqs = spacer_df[spacer_df['array_id'].isin(array_filt['array_id'])]

# Step 2 ─ bring the speciesname column across with a merge
seq_with_species = filt_seqs.merge(
    array_df[['array_id', 'speciesname', 'biosampleaccn']],   # only the columns we need
    on='array_id',
    how='left'                               # keep every spacer row, fill NaN if no match
)


species = "Porphyromonas gingivalis"
bioaccn = "SAMN03284263"
df_species = seq_with_species[
    (seq_with_species['speciesname'] == species) & 
    (seq_with_species['type'] == 'DR')]


samples = []

for index, row in df_species.iterrows():
    rna = str(row['sequence']).replace('T', 'U')
    if not rna.startswith('GUUGGA'):
        if row['biosampleaccn'] not in samples:
            samples.append(row['biosampleaccn'])
    else:
        pass


print(len(samples))

# for i in samples:
#     print(i)






# # Calculate the nucleotide content of each spacee
# df_nt_content = spacer_with_species

# for index, row in df_nt_content.iterrows():
#     spacer_seq = row['sequence']
#     nt_content = analysis_utils.nucleotide_content_calc(spacer_seq, nt='dna')
#     df_nt_content.loc[index, 'A'] = nt_content.get('A')
#     df_nt_content.loc[index, 'C'] = nt_content.get('C')
#     df_nt_content.loc[index, 'G'] = nt_content.get('G')
#     df_nt_content.loc[index, 'T'] = nt_content.get('T')


# repeat_lengths = repeats_with_species['sequence'].apply(len)
# repeats_with_species['repeat_len'] = repeat_lengths

# print(repeats_with_species)

# # Select only repeats with length 36
# repeats_filt = repeats_with_species[repeats_with_species['repeat_len'] == 36]









# home = Path("/home/unimelb.edu.au/rbengtsson/work/spacer_analysis")
# array_df = db_utils.read_table_from_db(f"{home}/prelim/analysis/preliminary_crisprs.csv")
# spacer_df = db_utils.read_table_from_db(f"{home}/prelim/analysis/all_spacers.csv")

# def calculate_gc_content(sequence: str) -> float:
#     """
#     Calculate the GC content of a DNA sequence.
    
#     Args:
#         sequence (str): The DNA sequence.
        
#     Returns:
#         float: The GC content as a percentage.
#     """
#     seq = Seq(sequence)
#     gc_count = seq.count("G") + seq.count("C")
#     return (gc_count / len(seq)) * 100 if len(seq) > 0 else 0.0


# fasta = sys.argv[1]

# for record in SeqIO.parse(fasta, "fasta"):
#     # Get the sequence from the record
#     sequence = str(record.seq)
    
#     # Calculate GC content
#     gc_content = calculate_gc_content(sequence)
    
#     # Print the result
#     print(f"GC content of {record.id}: {gc_content:.2f}%")


# if __name__ == "__main__":
#     main()