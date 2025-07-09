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


def retrieve_data_from_db(outdir: Path, query: str, database_name: str = 'cas13_bacterial_db') -> pd.DataFrame:
    # connect to the database
    creds = db_utils.get_credentials()['db_credentials']
    database = database_name

    pg_url = "postgresql+psycopg2://{0}:{1}@{2}:{3}/{4}".format(
                    creds['user'], creds['password'], creds['host'], creds['port'], database)

    engine = create_engine(pg_url, pool_pre_ping=True, echo=False)

    # Retrieve data from the database
    # with engine.connect() as conn:
    with engine.connect() as conn:
        df = pd.read_sql(text(query), conn.execution_options(stream_results=True))

    timestamp = datetime.now().strftime('%d-%m-%Y-%H-%M-%S')

    out_file = f"{outdir}/db_query_result_{timestamp}.csv"

    # write dataframe to csv
    df.to_csv(out_file, index=False)

    return df


def main():
    outdir = Path(__file__).resolve().parent
    # # query = """
    # #     SELECT *
    # #     FROM cas13b_crisprs
    # #     WHERE biosampleaccn IN ('SAMN02463902', 'SAMEA4362427', 'SAMEA104307710',
    # #            'SAMN03284263', 'SAMN03284264', 'SAMN03284261', 'SAMN03284262',
    # #            'SAMN02463748', 'SAMN00216833', 'SAMN03197167', 'SAMEA104224811',
    # #            'SAMN19926200', 'SAMN31809581', 'SAMN31809580'); 
    # # """
    # # query_results_df = retrieve_data_from_db(outdir, query)

    samples_list = Path("array_id.list")
    with open(samples_list, 'r') as f:
        samples = [line.strip() for line in f if line.strip()]
        
        df = db_utils.query_db_with_list(samples, 'spacer_table', 'array_id')
        subset = df[df['type'] == 'spacer']
        subset.to_csv(outdir / 'preliminary_spacer_seqs.csv', index=False)

    #     # df.to_csv(outdir / 'biosampleaccn_query_results.csv', index=False)
    
    # print(datetime.now().strftime("%d-%m-%Y-%H-%M-%S"))


if __name__ == "__main__":
    main()