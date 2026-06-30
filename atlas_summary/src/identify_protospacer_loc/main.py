"""
Main script to extract Cas13b spacers from the CrisprCas Atlas database and map them
to Phagescope database to identify protospacer region.
"""

from concurrent.futures import ProcessPoolExecutor
import pandas as pd
from pathlib import Path
from utils import extract_data_from_database, blast_search_batch, PhageBlast, get_terminase_features, reorient_coordinates



def extract_cas13b_spacers(outfile_path: Path, file_name: str):
    # Query to extract Cas13b spacers from atlas database
    query = """
            SELECT s.operon_id, s.spacer_length, s.spacer_sequence
            FROM spacer s
            JOIN cas c
                ON s.operon_id = c.operon_id
            WHERE c.gene_name ILIKE 'cas13b%'
                AND c.evalue = 0
                AND (s.spacer_length >= 25 OR s.spacer_length <= 35);
        """

    # write the results to a CSV file
    output_file = outfile_path / file_name
    df = extract_data_from_database(query)
    df.to_csv(output_file, index=False)



def blast_chunk_to_phagescope(chunk_df: pd.DataFrame):
    """
    Perform batch BLAST search for a chunk of spacers against the Phagescope database.
    """
    blast_db = Path('/phagescope_db/phagescope_db/data/phagescope/phagescope_db')
    output_path = Path("/phagescope_db/crisprcas_atlas_blastn/results/")

    blast_search_batch(chunk_df, blast_db, output_path, pident=85, cov_perc=95, num_threads=4)



def _process_blast_file(blast_file: Path) -> pd.DataFrame:
    """Parse a single BLAST output file and return a DataFrame of results."""
    rows = []
    with open(blast_file, 'r') as f:
        for line in f:
            phage_instance = PhageBlast.from_blast_line(line)
            rows.append(phage_instance.__dict__)
    return pd.DataFrame(rows)


def generate_results_file(blast_output_dir: Path, result_output_dir: Path, max_workers: int = 4):
    """
    Generate a results file summarising the BLAST search results, including protospacer location and phage information.
    """
    blast_files = list(blast_output_dir.glob("batch_*_blastn.txt"))
    print(f"Processing {len(blast_files)} BLAST output files with {max_workers} workers")

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        dfs = list(executor.map(_process_blast_file, blast_files))

    df = pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()
    df.to_csv(f"{result_output_dir}/cas13b_spacer_blastn_results_raw.csv", index=False)




def main():
    output_path = Path("/home/unimelb.edu.au/rbengtsson/work/spacer_analysis/atlas_summary/data/raw/")
    file_name = "cas13b_spacers.csv"
    # extract_cas13b_spacers(output_path, file_name)

    #### Filter spacers length between 25 and 35 nt
    cas13b_spacers_raw = output_path / file_name
    cas13b_spacers_df = pd.read_csv(cas13b_spacers_raw)
    spacers_filt = cas13b_spacers_df[cas13b_spacers_df['spacer_length'].between(25, 35)]


    # #### Chunk spacers and run blastn in parallel (one chunk per worker)
    # chunk_size = 500
    # chunks = [
    #     spacers_filt.iloc[i:i + chunk_size]
    #     for i in range(0, len(spacers_filt), chunk_size)
    # ]
    # print(f"Processing {len(spacers_filt)} spacers in {len(chunks)} chunks of {chunk_size}")

    # # # TEST: run a single chunk of 10 spacers synchronously (no executor)
    # # test_chunk = spacers_filt.iloc[:10]
    # # blast_chunk_to_phagescope(test_chunk)

    # with ProcessPoolExecutor(max_workers=4) as executor:
    #     executor.map(blast_chunk_to_phagescope, chunks)
    

    ### Generate results file
    blast_output_dir = Path("/phagescope_db/crisprcas_atlas_blastn/results/")
    result_output_dir = Path("/home/unimelb.edu.au/rbengtsson/work/spacer_analysis/atlas_summary/data/raw/")
    # generate_results_file(blast_output_dir, result_output_dir)


    ### Get terminase features for phages with protospacer hits
    phage_hits_file = Path(result_output_dir) / "cas13b_spacer_blastn_results_raw.csv"
    phage_hits_df = pd.read_csv(phage_hits_file)
    phage_ids = tuple(phage_hits_df['phage_id'].unique())
    terminase_output_dir = Path("/home/unimelb.edu.au/rbengtsson/work/spacer_analysis/atlas_summary/outputs/")
    # phage_ids_tmp = tuple(['Han_2018_ERR1398082_NODE_13_length_96793_cov_101.246098'])
    # get_terminase_features(phage_ids, terminase_output_dir)
    

    ### Merge terminase and blastn results
    terminase_file = terminase_output_dir / "terminase_features.csv"
    terminase_df = pd.read_csv(terminase_file)

    blastn_file = result_output_dir / "cas13b_spacer_blastn_results_raw.csv"
    blastn_df = pd.read_csv(blastn_file)

    merged = blastn_df.merge(terminase_df, on='phage_id', how='left')
    merged = merged.dropna(subset=['product'])

    ### reorient protospacer location relative to terminase start
    merged['reoriented_phagestart'] = merged.apply(
        lambda row: reorient_coordinates(row['phagestart'], row['start'], row['phagelength']),
        axis=1
    )
    merged['reoriented_phageend'] = merged.apply(
        lambda row: reorient_coordinates(row['phageend'], row['start'], row['phagelength']),
        axis=1
    )

    ### save merged results with reoriented coordinates
    merged_output_file = result_output_dir / "cas13b_spacer_blastn_results.csv"
    merged.to_csv(merged_output_file, index=False)




if __name__ == "__main__":
    main()