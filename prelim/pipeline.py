import sys
from pathlib import Path
import os
import uuid
import typing as t
from typing import Optional
from sqlalchemy import create_engine, text
from datetime import datetime

# Add the directory containing your other modules
target_dir = Path(__file__).resolve().parent.parent
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



def fetch_data(query_data: pd.DataFrame, outdir: Path) -> None:

    for index, row in query_data.iterrows():
        biosample_accession = row['biosampleaccn']
        ftp_path = row['ftppath_genbank']
        gff_file_path = genome_utils.get_gff_from_ncbi(ftp_path, biosample_accession, outdir)
        if gff_file_path is None:
            pass
        elif os.path.exists(gff_file_path):
            genome_utils.get_fasta_from_ncbi(ftp_path, biosample_accession, outdir)
            genome_utils.get_gbff_from_ncbi(ftp_path, biosample_accession, outdir)
            print(f"Fetched data for {biosample_accession}")



def cas13b_search(dataframe: pd, outdir: Path, cas_type: str = "cas13b") -> None:

    df_output = pd.DataFrame(columns=[
        'biosample_accession',
        'organism',
        'contig_id',
        'cas_gene',
        'locus_tag',
        'orientation',
        'start',
        'end'
    ])

    rows = []

    for index, row in dataframe.iterrows():
        biosample_accession = row['biosampleaccn']
        print(biosample_accession)
        gbff_file = outdir / f"{biosample_accession}.gbff"
        gff_file = outdir / f"{biosample_accession}.gff"
        fasta_file = outdir / f"{biosample_accession}.fna"
        cls = CasContig.get_cas_info(gbff_file, gff_file)

        if cls is not None:
            # write to pandas dataframe
            df_row = {
                'biosample_accession': biosample_accession,
                'organism': row['organism'],
                'contig_id': cls.contig_id,
                'cas_gene': cls.cas_type,
                'locus_tag': cls.locus_tag,
                'orientation': cls.orientation,
                'start': cls.start,
                'end': cls.end
                }
            rows.append(df_row)
            # Extract contig sequence
            contig_outdir = outdir / f"{biosample_accession}_contig"
            if not contig_outdir.exists():
                contig_outdir.mkdir(parents=True, exist_ok=True)
            CasContig.get_cas_contig(cls.contig_id, fasta_file, contig_outdir)
        else:
            print(f"No {cas_type} found for {biosample_accession} in {gbff_file} and {gff_file}")
            df_row = {
                'biosample_accession': biosample_accession,
                'organism': row['organism'],
                'contig_id': None,
                'cas_gene': None,
                'locus_tag': None,
                'orientation': None,
                'start': None,
                'end': None
            }
            rows.append(df_row)

    df_output = pd.DataFrame(rows)

    # Save the DataFrame to a CSV file
    output_file = outdir / "cas13b_contigs.csv"
    df_output.to_csv(output_file, index=False)

    return df_output



def run_crispr_prediction(dataframe: pd, outdir: Path) -> None:
    """
    Run CRISPR prediction on the given dataframe containing biosample accessions.
    The results will be saved in the specified output directory.
    """
    for index, row in dataframe.iterrows():
        if pd.isna(row['cas_gene']) or row['cas_gene'] == "":
            pass
        else:
            bioaccession = row['biosample_accession']
            fasta_file = outdir / f"{bioaccession}/{bioaccession}.fna"
            # print(fasta_file)
            CrisprDetection.run_crispridentify(fasta_file, outdir)



def upload_result(contig_file: Path, outdir: Path) -> None:
    # read the cas13b contigs file
    df = db_utils.read_table_from_db(contig_file)


    for index, row in df.iterrows():
        # check if the cas_gene is not NaN or empty
        if pd.isna(row['cas_gene']) or row['cas_gene'] == "":
            pass
        else:
            # Run data extraction
            bioaccession = row['biosample_accession']
            gbff_file = outdir / f"{bioaccession}.gbff"
            gff_file = outdir / f"{bioaccession}.gff"
            result_folder = Path(f"/{outdir}/{bioaccession}_crispridentify/{bioaccession}")

            candidates = ArrayCandidate.complie_crispr_arrays(result_folder, gbff_file, gff_file, outdir)
            if candidates is None:
                print(f"No candidates found for {bioaccession}. Skipping.")
                continue
            
            new_candidates_list = []    

            for c in candidates:
                candidate_new = ArrayCandidate.get_spacer_dr_seq(c, result_folder)
                new_candidates_list.append(ArrayCandidate.get_spacer_len_var(candidate_new))
            
            # # Filter candidates based on average DR and spacer lengths
            # # This will return a list of ArrayCandidate objects
            filter_candidates = ArrayCandidate.filter_candidates(new_candidates_list, avg_dr_len=36, avg_spacer_len=30)

            # Get the CasContig class instance for the given gbff and gff files
            # This will extract the contig_id and other relevant information
            cls = CasContig.get_cas_info(gbff_file, gff_file)

            array_df = db_utils.generate_array_table(filter_candidates, cls)
            spacer_df = db_utils.generate_spacer_table(filter_candidates)

            # Upload the dataframes to the SQL database
            db_utils.upload_arraytable_to_sql(array_df, "cas13_bacterial_db", "cas13b_crisprs")
            db_utils.upload_spacertable_to_sql(spacer_df, "cas13_bacterial_db", "spacer_table")



def main():

    ### Step 1: Set up the output directory
    # timestamp = datetime.now().strftime('%d-%m-%Y-%H-%M-%S')
    timestamp = "29-07-2025-16-47-04"
    outdir = Path(f"/spacers_db/prediction_results_{timestamp}")
    if not outdir.exists():
        outdir.mkdir(parents=True, exist_ok=True)
    

    # ## Step 2: Define the query to retrieve data from the database
    # # Define the query to retrieve data from the database
    # # set query
    # query = """
    #     SELECT n.biosampleaccn,
    #         n.organism,
    #         s.sequencingtechnology,
    #         s.genomecoverage,
    #         n.ftppath_genbank
    #     FROM ncbi_assembly AS n
    #     JOIN seq_platform AS s ON n.assemblyname = s.assemblyname
    #     WHERE n.biosampleaccn = 'SAMN03382583';
    # """

    # # Retrieve data from the database
    # query_results_df = retrieve_data_from_db(outdir, query)


    # ## Step 3: Fetch data from NCBI based on the query results
    # # Fetch data from NCBI based on the query results
    # fetch_data(query_results_df, outdir)

    # # Check if the outdir is empty
    # if not os.listdir(outdir):
    #     print("No data fetched from NCBI. Exiting.")
    #     return

    # # # unzip the files in the outdir
    # genome_utils.unzip_files(outdir)


    # ### Step 4: Search for cas13b in the retrieved data
    # # Search for cas13b in the retrieved data
    # # query_results_df = db_utils.read_table_from_db(f"/spacers_db/db_query_result.csv")
    # cas13b_df = cas13b_search(query_results_df, outdir)



    # ### Step 5: Run CRISPR prediction on the cas13b contigs
    # # # Place holder for running CRISPR prediction
    # # run_crispr_prediction(cas13b_df, outdir)


    # ### Step 6: Upload the results to the database
    data_file = outdir / "cas13b_contigs.csv"
    # data_file = outdir / "test.csv"
    upload_result(f"{data_file}", outdir)


    # # ### Step 7: Clean database to remove duplicates
    db_utils.db_clean_duplicates()
    
if __name__ == "__main__":
    main()