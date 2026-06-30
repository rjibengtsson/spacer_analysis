from concurrent.futures import ProcessPoolExecutor
import pandas as pd
from pathlib import Path


def _process_blast_file(blast_file: Path) -> pd.DataFrame:
    
    df = pd.read_csv(blast_file)


    best = (
        df
        .groupby('sequence_id', as_index=False)
        .agg(best_bitscore=('bitscore', 'max'), best_evalue=('evalue', 'min'))
    )


    filtered_df = filtered_df.merge(best, on='sequence_id')
    filtered_df = filtered_df[
        (filtered_df['bitscore'] == filtered_df['best_bitscore']) &
        (filtered_df['evalue'] == filtered_df['best_evalue'])
    ].drop(columns=['best_bitscore', 'best_evalue'])


    return filtered_df


def main():
    result_output_dir = Path("/home/unimelb.edu.au/rbengtsson/work/spacer_analysis/atlas_summary/data/raw/")
    blastn_file = result_output_dir / "cas13b_spacer_blastn_results_raw.csv"    
    terminase_file = Path("/home/unimelb.edu.au/rbengtsson/work/spacer_analysis/atlas_summary/outputs/terminase_features.csv")
    
    blastn_df = pd.read_csv(blastn_file, header=0)
    terminase_df = pd.read_csv(terminase_file, header=0)

    merged = blastn_df.merge(terminase_df, on='phage_id', how='left')
    merged = merged.dropna(subset=['product'])
    print(len(merged))
    
    # filtered_results = _process_blast_file(blastn_file)
    # print(filtered_results.head())



if __name__ == "__main__":
    main()