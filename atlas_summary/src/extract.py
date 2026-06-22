"""
This module contains functions to extract data from the CrisprCas Atlas SQL database.
"""

import json
import typing as t
from pathlib import Path
import psycopg2
import pandas as pd


def get_db_credentials() -> t.Dict[str, str]:
    """
    Retrieves database credentials from a JSON file.
    """
    credentials_path = "/home/unimelb.edu.au/rbengtsson/work/spacer_analysis/modules/db_credentials.json"
    with open(credentials_path, "r") as f:
        credentials = json.load(f)
    
    return credentials


def extract_data_from_database(query: str, DATABASE: str = 'crisprcas_atlas') -> pd.DataFrame:
    """
    Connects to the PostgreSQL database and executes the provided SQL query.

    Parameters:
    - query: The SQL query to execute.
    - DATABASE: The name of the database to connect to.


    Returns:
    A pandas DataFrame containing the results of the query.
    """
    credentials = get_db_credentials()
    creds = credentials['db_credentials']

    try:
        # Establish a connection to the database
        connection = psycopg2.connect(
            host=creds['host'],
            database=DATABASE,
            user=creds['user'],
            password=creds['password']
        )
    
        # Execute the query and fetch the results into a DataFrame
        df = pd.read_sql_query(query, connection)
        
        return df

    except Exception as e:
        print(f"An error occurred while connecting to the database or executing the query: {e}")
        return pd.DataFrame()  # Return an empty DataFrame in case of error

    finally:
        if 'connection' in locals():
            connection.close()



def main():
    # Example usage
    query = """
        SELECT operon_id, gene_name, hmm_name, evalue, score, truncated, length
        FROM cas
        WHERE gene_name ILIKE ANY (ARRAY['Cas9%', 'Cas12%', 'Cas13%']);
        """


    output_path = Path("/home/unimelb.edu.au/rbengtsson/work/spacer_analysis/atlas_summary/data/raw/")
    file_name = "class2_cas_raw.csv"


    # write the results to a CSV file
    output_file = output_path / file_name
    df = extract_data_from_database(query)
    df.to_csv(output_file, index=False)


if __name__ == "__main__":
    main()