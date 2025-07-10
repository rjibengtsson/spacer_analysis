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


home = Path("/home/unimelb.edu.au/rbengtsson/work/spacer_analysis")
array_df = db_utils.read_table_from_db(f"{home}/prelim/analysis/preliminary_crisprs.csv")
spacer_df = db_utils.read_table_from_db(f"{home}/prelim/analysis/all_spacers.csv")

# filter away duplicates
array_filt = array_df.drop_duplicates(subset=['array_id'])

# Distingush between categories
array_filt = array_filt[
    (array_filt['category'] != 'Alternative') &
    (array_filt['spacer_len_var'] < 1) &
    (array_filt['filter'] == 'PASSED')
]

# Step 1 ─ keep only spacer rows whose array_id survived your earlier filtering
good_spacers = spacer_df[
    spacer_df['array_id'].isin(array_filt['array_id']) &
    (spacer_df['type'] == 'spacer')
]

# Step 2 ─ bring the speciesname column across with a merge
spacer_with_species = good_spacers.merge(
    array_df[['array_id', 'speciesname']],   # only the columns we need
    on='array_id',
    how='left'                               # keep every spacer row, fill NaN if no match
)

# Calculate the nucleotide content of each spacee
df_nt_content = spacer_with_species.head(10)

for index, row in df_nt_content.iterrows():
    spacer_seq = row['sequence']
    targetmRNA = analysis_utils.transcribe_dna(spacer_seq)
    # print(f"target mRNA seq:\t{targetmRNA}")
    print(f"crRNA seq:\t\t{analysis_utils.reverse_complement(targetmRNA)}")
    # exit(0)
    # nt_content = analysis_utils.nucleotide_content_calc(spacer_seq)
    # df_nt_content.loc[index, 'A'] = nt_content.get('A')
    # df_nt_content.loc[index, 'C'] = nt_content.get('C')
    # df_nt_content.loc[index, 'G'] = nt_content.get('G')
    # df_nt_content.loc[index, 'T'] = nt_content.get('T')

# df_nt_content.to_csv("prelim/analysis/spacer_nt_content.csv", index=False)











# def retrieve_data_from_db(outdir: Path, query: str, database_name: str = 'cas13_bacterial_db') -> pd.DataFrame:
#     # connect to the database
#     creds = db_utils.get_credentials()['db_credentials']
#     database = database_name

#     pg_url = "postgresql+psycopg2://{0}:{1}@{2}:{3}/{4}".format(
#                     creds['user'], creds['password'], creds['host'], creds['port'], database)

#     engine = create_engine(pg_url, pool_pre_ping=True, echo=False)

#     # Retrieve data from the database
#     # with engine.connect() as conn:
#     with engine.connect() as conn:
#         df = pd.read_sql(text(query), conn.execution_options(stream_results=True))

#     timestamp = datetime.now().strftime('%d-%m-%Y-%H-%M-%S')

#     out_file = f"{outdir}/db_query_result_{timestamp}.csv"

#     # write dataframe to csv
#     df.to_csv(out_file, index=False)

#     return df


# def main():
#     outdir = Path(__file__).resolve().parent
#     # # query = """
#     # #     SELECT *
#     # #     FROM cas13b_crisprs
#     # #     WHERE biosampleaccn IN ('SAMN02463902', 'SAMEA4362427', 'SAMEA104307710',
#     # #            'SAMN03284263', 'SAMN03284264', 'SAMN03284261', 'SAMN03284262',
#     # #            'SAMN02463748', 'SAMN00216833', 'SAMN03197167', 'SAMEA104224811',
#     # #            'SAMN19926200', 'SAMN31809581', 'SAMN31809580'); 
#     # # """
#     # # query_results_df = retrieve_data_from_db(outdir, query)

#     samples_list = Path("array_id.list")
#     with open(samples_list, 'r') as f:
#         samples = [line.strip() for line in f if line.strip()]
        
#         df = db_utils.query_db_with_list(samples, 'spacer_table', 'array_id')
#         subset = df[df['type'] == 'spacer']
#         subset.to_csv(outdir / 'preliminary_spacer_seqs.csv', index=False)

#     #     # df.to_csv(outdir / 'biosampleaccn_query_results.csv', index=False)
    
#     # print(datetime.now().strftime("%d-%m-%Y-%H-%M-%S"))


# if __name__ == "__main__":
#     main()