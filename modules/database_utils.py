import pandas as pd
import uuid
from datetime import datetime
import typing as t
from dataclasses import dataclass
from pathlib import Path
from decimal import Decimal
from sqlalchemy import (create_engine, MetaData, Table, select, or_, inspect, 
                        Column, Integer, String, Text, Numeric, Float)
from sqlalchemy.dialects.postgresql import NUMERIC
from sqlalchemy import text as sql_text
import json
import os


@dataclass
class Allthebacteria:
    """
    Class to handle database operations for Allthebacteria.
    """
    
    @staticmethod
    def get_connection() -> create_engine:
        """
        Establishes a connection to the Allthebacteria database using credentials from a JSON file.
        """
        database = 'allthebacteria'
        credentials = get_credentials()
        creds = credentials['db_credentials']
        return create_engine(
            url="postgresql://{0}:{1}@{2}:{3}/{4}".format(
                creds['user'], creds['password'], creds['host'], creds['port'], database
            )
        )
    
    @staticmethod
    def retrieve_metadata(sample_accessions: list) -> pd.DataFrame:
        """
        Retrieves metadata for a given sample accession from the Allthebacteria database.

        args:
            sample_accession (list): The sample accession to retrieve metadata for.
        
        returns:
            pd.DataFrame: A DataFrame containing the metadata.
        """

        engine = Allthebacteria.get_connection()
        metadata = MetaData()
        assembly = Table("assembly", metadata, autoload_with=engine)
        stats = Table("assembly_stats", metadata, autoload_with=engine)
        
        filters = or_(*[assembly.c.sample_accession == accession for accession in sample_accessions])

        # SELECT columns from assembly and assembly_stats tables
        query = (
            select(
            assembly.c.sample_accession,
            assembly.c.run_accession,
            assembly.c.assembly_accession,
            assembly.c.asm_fasta_on_osf,
            stats.c.total_length,
            stats.c.number,
            stats.c.n50,
            stats.c.n50n
            )
            .select_from(assembly.join(stats, assembly.c.sample_accession == stats.c.sample_accession))
            .where(
            filters,
            assembly.c.asm_fasta_on_osf == '1'
            )
        )

        with engine.connect() as conn:
            result = conn.execute(query)
            df = pd.DataFrame(result.fetchall(), columns=result.keys())
        
        return df


    def check_SRA_df(SRA_file: Path) -> t.Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Checks the SRA DataFrame for sample accessions that are present in the Allthebacteria database.

        args:
            SRA_file (Path): Path to the SRA CSV file.
        
        returns:
            ra_df: A DataFrame indicating sample accessions that are present in Allthebacteria.
            atb_merged_df: A DataFrame containing metadata from Allthebacteria for the matching sample accessions.
        """
        
        sra_df = read_table_from_db(SRA_file)
        sample_accessions = sra_df["biosample"].to_list()
        atb_df = Allthebacteria.retrieve_metadata(sample_accessions)

        # Check if sample_accession is in the Allthebacteria database
        sra_df['in_allthebacteria'] = sra_df['biosample'].isin(atb_df['sample_accession'])
        
        # Select relevant columns from sra_df
        sra_selected = sra_df[['biosample', 'organism', 'instrument_model', 'title']]

        # Merge sra_df with atb_df on biosample/sample_accession, keeping only matches
        atb_merged_df = pd.merge(
            sra_selected,
            atb_df,
            right_on='sample_accession',
            left_on='biosample',
            how='inner'
        )

        atb_merged_df = atb_merged_df.drop(columns=['biosample'])

        return sra_df, atb_merged_df


@dataclass
class CasBacteriaDB:
    """
    Class to handle  database operations for Cas13 bacteria db.
    """

    @staticmethod
    def get_connection() -> create_engine:
        """
        Establishes a connection to the CasDatabase using credentials from a JSON file.
        """
        database = 'cas13_bacterial_db'
        credentials = get_credentials()
        creds = credentials['db_credentials']
        return create_engine(
            url="postgresql://{0}:{1}@{2}:{3}/{4}".format(
                creds['user'], creds['password'], creds['host'], creds['port'], database
            )
        )
    
    
    def query_assembly(self, sample_accession: str) -> pd.DataFrame:
        """
        Queries the assembly table for a given sample accession.

        args:
            sample_accession (str): The sample accession to query.
        
        returns:
            pd.DataFrame: A DataFrame containing the assembly data for the sample accession.
        """
        engine = self.get_connection()
        metadata = MetaData()
        assembly = Table("assembly", metadata, autoload_with=engine)

        query = select(assembly).where(assembly.c.sample_accession == sample_accession)

        with engine.connect() as conn:
            result = conn.execute(query)
            df = pd.DataFrame(result.fetchall(), columns=result.keys())
        
        return df



def get_connection(database_name) -> create_engine:
    """
    Establishes a connection to the Allthebacteria database using credentials from a JSON file.
    """
    database = database_name
    credentials = get_credentials()
    creds = credentials['db_credentials']
    return create_engine(
        url="postgresql://{0}:{1}@{2}:{3}/{4}".format(
            creds['user'], creds['password'], creds['host'], creds['port'], database
        )
    )



def get_credentials() -> t.Dict[str, t.Any]:
    """
    Retrieves database credentials from a JSON file.
    Expects a file named 'db_credentials.json' in the same directory as this script.
    """

    credentials_path = Path(os.path.dirname(__file__)) / "db_credentials.json"
    with open(credentials_path, "r") as f:
        credentials = json.load(f)
    return credentials



def read_table_from_db(file: Path) -> pd.DataFrame:
    """
    Reads a table from a database and returns it as a pandas DataFrame.

    args:
        file (Path): Path to the database file.
    """
    df = pd.read_csv(file, header=0)

    return df



def generate_array_table(candidates: t.List[t.Any], cls: t.Any) -> pd.DataFrame:
    """
    Generates a DataFrame from a list of ArrayCandidate instances.

    args:
        candidates (list): List of ArrayCandidate instances.
        cls (type): The class type of the ArrayCandidate instances.
    
    returns:
        pd.DataFrame: A DataFrame containing the array data.
    """
    array_table = pd.DataFrame(columns=['array_id', 'biosampleaccn', 'contig_id',
                                    'arraystart', 'arrayend', 'arraylen', 'avgrepeatlen', 
                                    'avgspacerlen','number_of_spacers', 'dist_to_cas13b', 
                                    'orientation', 'category', 'score', 'filter', 'date'])


    for candidate in candidates:
        new_row = {
            'array_id': candidate.instance_id,
            'biosampleaccn': candidate.biosample_accn,
            'contig_id': cls.contig_id,
            'arraystart': candidate.arry_start,
            'arrayend': candidate.arry_end,
            'arraylen': candidate.length,
            'avgrepeatlen': candidate.avg_repeat_length,
            'avgspacerlen': candidate.avg_spacer_length,
            'number_of_spacers': candidate.num_spacers,
            'dist_to_cas13b': candidate.dist_to_cas,
            'orientation': candidate.orientation,
            'category': candidate.category,
            'score': candidate.score,
            'filter': candidate.filter,
            'date': datetime.now().strftime("%d-%m-%Y")
        }
        # array_table = pd.concat([array_table, pd.DataFrame([new_row])], ignore_index=True)
        new_df = pd.DataFrame([new_row])
        if not array_table.empty and not array_table.isna().all().all():
            array_table = pd.concat([array_table, new_df], ignore_index=True)
        else:
            array_table = new_df
    
    # Convert the score column to Decimal (as strings to preserve exact digits)
    array_table["score"] = array_table["score"].astype(float)

    return array_table



def generate_spacer_table(candidates: t.List[t.Any]) -> pd.DataFrame:
    """
    Generates a DataFrame from a list of ArrayCandidate instances, specifically for spacers.

    args:
        candidates (list): List of ArrayCandidate instances.
    
    returns:
        pd.DataFrame: A DataFrame containing the spacer data.
    """
    spacer_table = pd.DataFrame(columns=['sequence_id', 'array_id', 'type', 'position', 'sequence'])

    # Append spacer table data
    for candidate in candidates:
    
        # Fill in spacer sequence
        pos_counter = 0
        for seq in candidate.spacerseq:
            pos_counter += 1
            new_row = {
                'sequence_id': str(uuid.uuid4()),  # Generate a unique ID for each spacer
                'array_id': candidate.instance_id,
                'type': 'spacer',  # Assuming all spacers are of type 'spacer'
                'position': pos_counter,
                'sequence': seq
            }
            # spacer_table = pd.concat([spacer_table, pd.DataFrame([new_spacer_row])], ignore_index=True)
            new_df = pd.DataFrame([new_row])
            if not spacer_table.empty and not spacer_table.isna().all().all():
                spacer_table = pd.concat([spacer_table, new_df], ignore_index=True)
            else:
                spacer_table = new_df


        # Fill in DR sequence
        pos_counter = 0
        for seq in candidate.drseq:
            pos_counter += 1
            new_row = {
                'sequence_id': str(uuid.uuid4()),  # Generate a unique ID for each spacer
                'array_id': candidate.instance_id,
                'type': 'DR',  # Assuming all spacers are of type 'spacer'
                'position': pos_counter,
                'sequence': seq
            }
            new_df = pd.DataFrame([new_row])
            if not spacer_table.empty and not spacer_table.isna().all().all():
                spacer_table = pd.concat([spacer_table, new_df], ignore_index=True)
            else:
                spacer_table = new_df        

    return spacer_table



def upload_arraytable_to_sql(df: pd.DataFrame, database_name: str, table_name: str) -> None:
    """
    Uploads a the Array table to a database table.

    args:
        df (pd.DataFrame): The DataFrame to upload.
        table_name (str): The name of the table in the database.
        db_url (str): The name of the database.
    """
    engine = get_connection(database_name)

    TABLE_NAME = table_name
    SCHEMA  = "public"

    
    insp = inspect(engine)


    # Check if the table exists
    if TABLE_NAME in insp.get_table_names(schema=SCHEMA):
        print(f"{SCHEMA}.{TABLE_NAME} already exists - skipping.")

        # Alter score column to NUMERIC(20,4) using raw SQL
        with engine.connect() as conn:
            conn.execute(sql_text(f"""
                ALTER TABLE {SCHEMA}.{TABLE_NAME}
                ALTER COLUMN score TYPE NUMERIC(20, 4)
                USING ROUND(score::numeric, 4);
            """))
    
    # Otherwise, create the table
    else:
        metadata = MetaData(schema=SCHEMA)
        
        array_table = Table(
            TABLE_NAME, metadata,
            Column("array_id", String(50), primary_key=True, nullable=False),
            Column("biosampleaccn", String, nullable=False),
            Column("contig_id", String, nullable=False),
            Column("arraystart", Integer, nullable=False),
            Column("arrayend", Integer, nullable=False),
            Column("arraylen", Integer, nullable=False),
            Column("avgrepeatlen", Integer, nullable=False),
            Column("avgspacerlen", Integer, nullable=False),
            Column("number_of_spacers", Integer, nullable=False),
            Column("dist_to_cas13b", Integer, nullable=False),
            Column("orientation", Text, nullable=False),
            Column("category", Text, nullable=False),
            Column("score", NUMERIC(20, 4), nullable=False),
            Column("filter", Text, nullable=False),
            Column("date", String(10), nullable=False)  # Date in 'dd-mm-yyyy' format
        )

        metadata.create_all(engine)          # creates only the defined table
        print(f"Created {SCHEMA}.{TABLE_NAME}.")

    # Upload the DataFrame to the database
    df.to_sql(
        TABLE_NAME,
        con=engine,
        if_exists='append',
        index=False,
        method=None,
        schema=SCHEMA,
        dtype={"score": NUMERIC(20, 4)}
    )
    print(f"Uploaded {len(df)} rows to {SCHEMA}.{TABLE_NAME}.")



def upload_spacertable_to_sql(df: pd.DataFrame, database_name: str, table_name: str) -> None:
    """
    Uploads Spacer table to a database table.

    args:
        df (pd.DataFrame): The DataFrame to upload.
        table_name (str): The name of the table in the database.
        db_url (str): The name of the database.
    """
    engine = get_connection(database_name)

    TABLE_NAME = table_name
    SCHEMA  = "public"

    
    insp = inspect(engine)


    # Check if the table exists
    if TABLE_NAME in insp.get_table_names(schema=SCHEMA):
        print(f"{SCHEMA}.{TABLE_NAME} already exists - skipping.")
    # Otherwise, create the table
    else:
        metadata = MetaData(schema=SCHEMA)
        
        array_table = Table(
            TABLE_NAME, metadata,
            Column("sequence_id", String(50), primary_key=True, nullable=False),
            Column("array_id", String(50), nullable=False),
            Column("type", Text, nullable=False),
            Column("position", Integer, nullable=False),
            Column("sequence", String, nullable=False),
        )

        metadata.create_all(engine)          # creates only the defined table
        print(f"Created {SCHEMA}.{TABLE_NAME}.")

    # Upload the DataFrame to the database
    df.to_sql(
        TABLE_NAME,
        con=engine,
        if_exists='append',
        index=False,
        method='multi',
        schema=SCHEMA
    )
    print(f"Uploaded {len(df)} rows to {SCHEMA}.{TABLE_NAME}.")