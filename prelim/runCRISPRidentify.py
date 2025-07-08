import pandas as pd
import sys, os
from pathlib import Path

# Add the directory containing your other modules
target_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(target_dir))

from modules.crispr_detecion import CrisprDetection


def read_table_from_db(file: Path) -> pd.DataFrame:
    """
    Reads a table from a database and returns it as a pandas DataFrame.

    args:
        file (Path): Path to the database file.
    """
    df = pd.read_csv(file, header=0)

    return df



def run_crispr_prediction(dataframe: pd, outdir: Path) -> None:
    """
    Run CRISPR prediction on the given dataframe containing biosample accessions.
    The results will be saved in the specified output directory.
    """
    for _, row in dataframe.iterrows():
        cas_gene = row.get('cas_gene')
        if not cas_gene or pd.isna(cas_gene):
            continue

        bioaccession = row['biosample_accession']
        fasta_paths = [
            outdir / f"{bioaccession}_contig/{bioaccession}.fna",
            outdir / f"{bioaccession}.fna"
        ]

        for fasta_file in fasta_paths:
            if fasta_file.exists():
                CrisprDetection.run_crispridentify(fasta_file, outdir)


def main():
    
    # timestamp = datetime.now().strftime('%d-%m-%Y-%H-%M-%S')
    timestamp = "07-07-2025-11-30-52"
    outdir = Path(f"/spacers_db/prediction_results_{timestamp}")
    if not outdir.exists():
        outdir.mkdir(parents=True, exist_ok=True)

    # Load the cas13b contigs dataframe
    cas13b_df = read_table_from_db(f"{outdir}/cas13b_contigs.csv")

    # Run CRISPR prediction
    run_crispr_prediction(cas13b_df, outdir)


if __name__ == "__main__":
    main()