"""
Create/update SQL tables for storing CRISPR-Cas operon data.
"""

import json
import pandas as pd
import uuid
from datetime import datetime
import typing as t
from dataclasses import dataclass
from pathlib import Path
from decimal import Decimal
from sqlalchemy import (BigInteger, Boolean, ForeignKey, ForeignKeyConstraint, Identity, create_engine, MetaData, Table, select, or_, inspect, 
                        Column, Integer, String, Float, UniqueConstraint)
from sqlalchemy.dialects.postgresql import NUMERIC
from sqlalchemy import text as sql_text



def get_db_credentials() -> t.Dict[str, str]:
    """
    Retrieves database credentials from a JSON file.
    """
    credentials_path = "/home/unimelb.edu.au/rbengtsson/work/spacer_analysis/modules/db_credentials.json"
    with open(credentials_path, "r") as f:
        credentials = json.load(f)
    
    return credentials



def get_connection():
    """
    Establish a connection to the PostgreSQL database using SQLAlchemy.
    """
    credentials = get_db_credentials()
    creds = credentials['db_credentials']
    user = creds["user"]
    password = creds["password"]
    host = creds["host"]
    port = creds["port"]
    database = "crisprcas_atlas"

    return create_engine(
        url="postgresql://{0}:{1}@{2}:{3}/{4}".format(
            creds['user'], creds['password'], creds['host'], creds['port'], database
        )
    )



def create_summary_table() -> None:
    """
    Create the 'summary' table in the SQL database.
    """
    engine = get_connection()
    
    TABLE_NAME = "summary"
    SCHEMA  = "public"
    metadata = MetaData(schema=SCHEMA)

    insp = inspect(engine)

    # Check if the table already exists
    if TABLE_NAME in insp.get_table_names(schema=SCHEMA):
        print(f"{SCHEMA}.{TABLE_NAME} already exists - skipping.")
    else:
        print(f"Creating table {SCHEMA}.{TABLE_NAME}...")
        summary_table = Table(
            TABLE_NAME,
            metadata,
            Column("operon_id", String, primary_key=True, unique=True),
            Column("subtype", String),
            Column("subtype_score", String),
            Column("operon_length", Integer),
            Column("n_crispr", Integer),
            Column("n_spacers", Integer),
            Column("n_tracr", Integer),
            Column("n_cas", Integer),
            Column("n_genes", Integer)
        )
        metadata.create_all(engine)
        print(f"Table {SCHEMA}.{TABLE_NAME} created successfully.")


def create_metadata_table() -> None:
    """
    Create the 'metadata' table in the SQL database.
    """
    engine = get_connection()
    
    TABLE_NAME = "metadata"
    SCHEMA  = "public"
    metadata_obj = MetaData(schema=SCHEMA)

    insp = inspect(engine)

    # Check if the table already exists
    if TABLE_NAME in insp.get_table_names(schema=SCHEMA):
        print(f"{SCHEMA}.{TABLE_NAME} already exists - skipping.")
    else:
        print(f"Creating table {SCHEMA}.{TABLE_NAME}...")

        # Reflect the referenced summary table into the same MetaData
        Table("summary", metadata_obj, autoload_with=engine, schema=SCHEMA)

        summary_table = Table(
            TABLE_NAME,
            metadata_obj,
            Column("operon_id", String, ForeignKey(f"{SCHEMA}.summary.operon_id", ondelete="CASCADE"), primary_key=True),
            Column("source_db", String),
            Column("assembly_type", String),
            Column("biosample_id", String),
            Column("sample_name", String),
            Column("taxonomy", String),
            Column("biome", String)
        )
        metadata_obj.create_all(engine)
        print(f"Table {SCHEMA}.{TABLE_NAME} created successfully.")



def create_tracr_table() -> None:
    """
    Create the 'tracr' table in the SQL database.
    """
    engine = get_connection()
    
    TABLE_NAME = "tracr"
    SCHEMA  = "public"
    tracr_obj = MetaData(schema=SCHEMA)

    insp = inspect(engine)

    # Check if the table already exists
    if TABLE_NAME in insp.get_table_names(schema=SCHEMA):
        print(f"{SCHEMA}.{TABLE_NAME} already exists - skipping.")
    else:
        print(f"Creating table {SCHEMA}.{TABLE_NAME}...")

        # Reflect the referenced summary table into the same MetaData
        Table("summary", tracr_obj, autoload_with=engine, schema=SCHEMA)

        tracr_table = Table(
            TABLE_NAME,
            tracr_obj,
            Column("operon_id", String, ForeignKey(f"{SCHEMA}.summary.operon_id", ondelete="CASCADE"), primary_key=True),
            Column("cm_id", String),
            Column("evalue", Float),
            Column("truncated", String),
            Column("gene_overlap", Integer),
            Column("terminator", Integer),
            Column("confidence", String),
            Column("tracr", String)
        )
        tracr_obj.create_all(engine)
        print(f"Table {SCHEMA}.{TABLE_NAME} created successfully.")



def create_repeat_table() -> None:
    """
    Create the 'repeat' table in the SQL database.
    """
    engine = get_connection()
    
    TABLE_NAME = "repeat"
    SCHEMA  = "public"
    repeat_obj = MetaData(schema=SCHEMA)

    insp = inspect(engine)

    # Check if the table already exists
    if TABLE_NAME in insp.get_table_names(schema=SCHEMA):
        print(f"{SCHEMA}.{TABLE_NAME} already exists - skipping.")
    else:
        print(f"Creating table {SCHEMA}.{TABLE_NAME}...")

        # Reflect the referenced summary table into the same MetaData
        Table("summary", repeat_obj, autoload_with=engine, schema=SCHEMA)

        repeat_table = Table(
            TABLE_NAME,
            repeat_obj,
            # Upload script-generated PK (NOT uploaded from record)
            Column("id", String, primary_key=True),
            # FK to parent
            Column("operon_id", String, ForeignKey(f"{SCHEMA}.summary.operon_id", ondelete="CASCADE"), nullable=False),
            Column("array_idx", Integer, nullable=False),
            Column("number_spacers", Integer),
            Column("repeat_seq_length", Integer),
            Column("repeat_sequence", String),
            UniqueConstraint("operon_id", "array_idx", name="uq_repeat_operon_array")
        )
        repeat_obj.create_all(engine)
        print(f"Table {SCHEMA}.{TABLE_NAME} created successfully.")



def create_spacer_table() -> None:
    """
    Create the 'spacer' table in the SQL database.
    """
    engine = get_connection()
    
    TABLE_NAME = "spacer"
    SCHEMA  = "public"
    spacer_obj = MetaData(schema=SCHEMA)

    insp = inspect(engine)

    # Check if the table already exists
    if TABLE_NAME in insp.get_table_names(schema=SCHEMA):
        print(f"{SCHEMA}.{TABLE_NAME} already exists - skipping.")
    else:
        print(f"Creating table {SCHEMA}.{TABLE_NAME}...")

        # Reflect the referenced repeat table into the same MetaData
        Table("repeat", spacer_obj, autoload_with=engine, schema=SCHEMA)

        spacer_table = Table(
            TABLE_NAME,
            spacer_obj,
            # Upload script generated PK (NOT uploaded from record)
            Column("id", String, primary_key=True),
            # FK to parent
            Column("operon_id", String, nullable=False),
            Column("array_idx", Integer, nullable=False),
            Column("spacer_position", Integer, nullable=False),
            Column("spacer_length", Integer),
            Column("spacer_sequence", String),
            ForeignKeyConstraint(["operon_id", "array_idx"],
            [f"{SCHEMA}.repeat.operon_id", f"{SCHEMA}.repeat.array_idx"],
            ondelete="CASCADE"),
            UniqueConstraint("operon_id", "array_idx", "spacer_position", name="uq_spacer_operon_array_position")
        )
        spacer_obj.create_all(engine)
        print(f"Table {SCHEMA}.{TABLE_NAME} created successfully.")



def create_cas_table() -> None:
    """
    Create the 'cas' table in the SQL database.
    """
    engine = get_connection()
    
    TABLE_NAME = "cas"
    SCHEMA  = "public"
    cas_obj = MetaData(schema=SCHEMA)

    insp = inspect(engine)

    # Check if the table already exists
    if TABLE_NAME in insp.get_table_names(schema=SCHEMA):
        print(f"{SCHEMA}.{TABLE_NAME} already exists - skipping.")
    else:
        print(f"Creating table {SCHEMA}.{TABLE_NAME}...")

        # Reflect the referenced summary table into the same MetaData
        Table("summary", cas_obj, autoload_with=engine, schema=SCHEMA)

        cas_table = Table(
            TABLE_NAME,
            cas_obj,
            # Upload script-generated PK (NOT uploaded from record)
            Column("id", String, primary_key=True),
            # FK to parent
            Column("operon_id", String, ForeignKey(f"{SCHEMA}.summary.operon_id", ondelete="CASCADE"), nullable=False),
            Column("gene_name", String),
            Column("hmm_name", String),
            Column("evalue", Float),
            Column("score", Float),
            Column("truncated", String),
            Column("length", Integer),
            Column("protein", String),
            UniqueConstraint("operon_id", "gene_name", "hmm_name", "evalue", "score", name="uq_cas_operon_gene_hmm")
        )
        cas_obj.create_all(engine)
        print(f"Table {SCHEMA}.{TABLE_NAME} created successfully.")


def create_operonid_spacerid_table() -> None:
    """
    Create the 'operonid_spacerid' table in the SQL database.
    """
    engine = get_connection()
    
    TABLE_NAME = "operonid_spacerid"
    SCHEMA  = "public"
    mapping_obj = MetaData(schema=SCHEMA)

    insp = inspect(engine)

    # Check if the table already exists
    if TABLE_NAME in insp.get_table_names(schema=SCHEMA):
        print(f"{SCHEMA}.{TABLE_NAME} already exists - skipping.")
    else:
        print(f"Creating table {SCHEMA}.{TABLE_NAME}...")

        # Reflect the referenced summary and spacer tables into the same MetaData
        Table("summary", mapping_obj, autoload_with=engine, schema=SCHEMA)
        Table("spacer", mapping_obj, autoload_with=engine, schema=SCHEMA)

        mapping_table = Table(
            TABLE_NAME,
            mapping_obj,
            Column("operon_id", String, ForeignKey(f"{SCHEMA}.summary.operon_id", ondelete="CASCADE"), primary_key=True),
            Column("spacer_id", String, ForeignKey(f"{SCHEMA}.spacer.id", ondelete="CASCADE"), primary_key=True)
        )
        mapping_obj.create_all(engine)
        print(f"Table {SCHEMA}.{TABLE_NAME} created successfully.")


def main():
    for create_table in (
        create_summary_table,
        create_metadata_table,
        create_tracr_table,
        create_repeat_table,
        create_spacer_table,
        create_operonid_spacerid_table,
        create_cas_table    
    ):
        create_table()



if __name__ == "__main__":
    main()
