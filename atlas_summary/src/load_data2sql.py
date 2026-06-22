"""
ETL script to load CRISPR-Cas operon data from JSONL files into SQL database.
"""

from concurrent.futures import ProcessPoolExecutor
from itertools import islice
import pandas as pd
import uuid
from datetime import datetime
import typing as t
from dataclasses import dataclass
from pathlib import Path
from decimal import Decimal
from sqlalchemy import (create_engine, MetaData, Table, exists, insert, select, or_, inspect, 
                        Column, Integer, String, Text, Numeric, Float)
from sqlalchemy.dialects.postgresql import insert as pg_insert
from sqlalchemy.dialects.postgresql import NUMERIC
from sqlalchemy import text as sql_text
import json
import os


Tables = t.Dict[str, pd.DataFrame]


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


def iter_jsonl(file_path: str) -> list:
    """
    Read a JSONL file where each line is an operon record, and return
    a list of these records.
    """
    with open(file_path, "r", encoding="utf-8") as fh:
        for line_no, line in enumerate(fh, start=1):
            line = line.strip()
            if not line:
                continue
            record = json.loads(line)
            yield record



def extract_record_info(record: t.Dict) -> Tables:
    """
    Convert a single operon record into multiple dataframes for different components.
    """
    operon_id = record.get("operon_id")

    # Summary / metadata (one-row tables)
    summary_df  = pd.DataFrame([{"operon_id": operon_id, **(record.get("summary") or {})}])
    metadata_df = pd.DataFrame([{"operon_id": operon_id, **(record.get("metadata") or {})}])


    # tracr can be missing / None / {}
    tracr = record.get("tracr")
    tracr_df = (
        pd.DataFrame([{"operon_id": operon_id, **tracr}])
        if isinstance(tracr, dict) and tracr
        else pd.DataFrame()
    )

    # CRISPR arrays -> repeats + spacers (many rows tables)
    repeat_rows: t.List[t.Dict] = []
    spacer_rows: t.List[t.Dict] = []

    for array_idx, array in enumerate(record.get("crispr", []), start=1):
        repeat_seq = array.get("crispr_repeat")
        repeat_idx = str(uuid.uuid4())
        if repeat_seq is not None:
            repeat_rows.append(
                {   
                    "id" : repeat_idx, # generate a unique repeat_id for this repeat
                    "operon_id": operon_id,
                    "array_idx": array_idx,
                    "number_spacers": len(array.get("crispr_spacers")[:-1] or []), # count spacers for this repeat
                    "repeat_seq_length": len(repeat_seq),
                    "repeat_sequence": repeat_seq,
                }
            )

        spacers = (array.get("crispr_spacers") or [])[:-1]  # exclude last spacer
        for spacer_pos, spacer_seq in enumerate(spacers, start=1):
            spacer_rows.append(
                {   
                    "id": str(uuid.uuid4()), # generate a unique spacer_id for this spacer
                    "operon_id": operon_id,
                    "array_idx": array_idx, # link to the parent repeat
                    "spacer_position": spacer_pos,
                    "spacer_length": len(spacer_seq),
                    "spacer_sequence": spacer_seq,
                }
            )

    repeat_df = pd.DataFrame(repeat_rows)
    spacer_df = pd.DataFrame(spacer_rows)

    # cas (many rows tables)
    cas_rows: t.List[t.Dict] = []
    for cas in record.get("cas", []) or []:
        cas_rows.append(
            {   
                "id": str(uuid.uuid4()),
                "operon_id": operon_id,
                "gene_name": cas.get("gene_name"),
                "hmm_name": cas.get("hmm_name"),
                "evalue": float(cas.get("evalue") or 0),
                "score": float(cas.get("score") or 0),
                "truncated": cas.get("truncated", False),
                "length": int(cas.get("length") or 0),
                "protein": cas.get("protein", ""),
            }
        )
    cas_df = pd.DataFrame(cas_rows)

    return {
        "summary": summary_df, # parent table with business key as operon_id
        "metadata": metadata_df,
        "tracr": tracr_df,
        "repeat": repeat_df,
        "spacer": spacer_df,
        "cas": cas_df,
    }



def get_reflected_tables(engine):
    """
    Reflect the existing tables from the database to get their metadata.
    """
    md = MetaData(schema="public")
    return {
        "summary": Table("summary", md, autoload_with=engine),
        "metadata": Table("metadata", md, autoload_with=engine),
        "tracr": Table("tracr", md, autoload_with=engine),
        "repeat": Table("repeat", md, autoload_with=engine),
        "spacer": Table("spacer", md, autoload_with=engine),
        "cas": Table("cas", md, autoload_with=engine),
    }


def upload_tables_to_sql(engine, sql_tables, tables: Tables) -> None:
    """
    Upload one operon's extracted tables.
    Assumes tables['summary'] contains exactly one row.
    """

    summary_df = tables['summary']
    if summary_df.empty:
        raise ValueError("summary dataframe is empty")    

    # operon_id is the PK for summary table, but foreign keys in child tables
    operon_id = summary_df["operon_id"].iloc[0]
    # print(summary_df.iloc[0].to_dict())

    summary_tbl = sql_tables["summary"]
    metadata_tbl = sql_tables["metadata"]
    tracr_tbl = sql_tables["tracr"]
    repeat_tbl = sql_tables["repeat"]
    spacer_tbl = sql_tables["spacer"]
    cas_tbl = sql_tables["cas"]

    with engine.begin() as conn:
        # 1) upsert summary table (parent table)
        # use operon_id as UNIQUE key, and get record.id for FK in child tables
        summary_payload = summary_df.drop(columns=["operon_id"]).iloc[0].to_dict()

        stmt = (
            pg_insert(summary_tbl)
            .values(operon_id=operon_id, **summary_payload)
            .on_conflict_do_update(
                index_elements=[summary_tbl.c.operon_id],
                set_=summary_payload,
            )
            .returning(summary_tbl.c.operon_id)
        )
        operon_id = conn.execute(stmt).scalar_one()


        # 2) upload 1:1 tables (metadata/tracr) — upsert on operon_id
        # this assumes metadata/tracr tables have a operon_id column as FK
        for name, tbl in [("metadata", metadata_tbl), ("tracr", tracr_tbl)]:
            df = tables[name]
            if df is None or df.empty:
                continue

            # Insert FK into payload for upsert
            payload = df.iloc[0].to_dict()
            payload["operon_id"] = operon_id

            stmt = (
                pg_insert(tbl)
                .values(**payload)
                .on_conflict_do_update(
                    index_elements=[tbl.c.operon_id],
                    set_=payload,
                )
            )
            conn.execute(stmt)        


        # 3) upload 1:many tables (repeat/spacer/cas) — delete existing rows for operon_id, then insert new rows
        # this assumes these tables have operon_id and appropriate UNIQUE constraints
        # repeat unique: (operon_id, array_idx)
        # spacer unique: (operon_id, repeat_id, spacer_idx)
        # cas unique: (operon_id, cas_idx)

        # repeat (1:many)
        repeat_df = tables.get("repeat")
        if repeat_df is not None and not repeat_df.empty:
            repeat_df = repeat_df.copy()
            repeat_df["operon_id"] = operon_id

            rows = repeat_df.to_dict(orient="records")
            stmt = pg_insert(repeat_tbl).values(rows)
            stmt = stmt.on_conflict_do_update(
                index_elements=[repeat_tbl.c.operon_id, repeat_tbl.c.array_idx],
                set_={
                    "number_spacers": stmt.excluded.number_spacers,
                    "repeat_seq_length": stmt.excluded.repeat_seq_length,
                    "repeat_sequence": stmt.excluded.repeat_sequence,
                },
            )
            conn.execute(stmt)

        # spacer (1:many)
        spacer_df = tables.get("spacer")
        if spacer_df is not None and not spacer_df.empty:
            spacer_df = spacer_df.copy()
            spacer_df["operon_id"] = operon_id

            rows = spacer_df.to_dict(orient="records")
            stmt = pg_insert(spacer_tbl).values(rows)
            stmt = stmt.on_conflict_do_update(
                index_elements=[
                    spacer_tbl.c.operon_id,
                    spacer_tbl.c.array_idx,
                    spacer_tbl.c.spacer_position,
                ],
                set_={
                    "spacer_length": stmt.excluded.spacer_length,
                    "spacer_sequence": stmt.excluded.spacer_sequence,
                },
            )
            conn.execute(stmt)

        # cas (1:many)
        cas_df = tables.get("cas")
        if cas_df is not None and not cas_df.empty:
            cas_df = cas_df.copy()
            cas_df["operon_id"] = operon_id

            rows = cas_df.to_dict(orient="records")
            stmt = pg_insert(cas_tbl).values(rows)
            stmt = stmt.on_conflict_do_update(
                index_elements=[
                    cas_tbl.c.operon_id,
                    cas_tbl.c.gene_name,
                    cas_tbl.c.hmm_name,
                ],
                set_={
                    "evalue": stmt.excluded.evalue,
                    "score": stmt.excluded.score,
                    "truncated": stmt.excluded.truncated,
                    "length": stmt.excluded.length,
                    "protein": stmt.excluded.protein,
                },
            )
            conn.execute(stmt)

        print(f"upserted; operon.id={operon_id}. Children loaded.")


def build_update_set(stmt, table, exclude=("operon_id",)):
    return {
        c.name: getattr(stmt.excluded, c.name)
        for c in table.columns
        if c.name not in exclude
    }


def upload_chunk_to_sql(engine, sql_tables, tables: Tables) -> None:
    """
    Bulk upsert a whole chunk of dataframes.
    Assumes operon_id is present in all tables.    
    """

    summary_tbl = sql_tables["summary"]
    metadata_tbl = sql_tables["metadata"]
    tracr_tbl = sql_tables["tracr"]
    repeat_tbl = sql_tables["repeat"]
    spacer_tbl = sql_tables["spacer"]
    cas_tbl = sql_tables["cas"]

    with engine.begin() as conn:
        # summary
        summary_df = tables["summary"]
        if summary_df is not None and not summary_df.empty:
            rows = summary_df.to_dict(orient="records")
            stmt = pg_insert(summary_tbl).values(rows)
            stmt = stmt.on_conflict_do_update(
                index_elements=[summary_tbl.c.operon_id],
                set_=build_update_set(stmt, summary_tbl),
            )
            conn.execute(stmt)

        # metadata
        metadata_df = tables["metadata"]
        if metadata_df is not None and not metadata_df.empty:
            rows = metadata_df.to_dict(orient="records")
            stmt = pg_insert(metadata_tbl).values(rows)
            stmt = stmt.on_conflict_do_update(
                index_elements=[metadata_tbl.c.operon_id],
                set_=build_update_set(stmt, metadata_tbl),
            )
            conn.execute(stmt)

        # tracr
        tracr_df = tables["tracr"]
        if tracr_df is not None and not tracr_df.empty:
            rows = tracr_df.to_dict(orient="records")
            stmt = pg_insert(tracr_tbl).values(rows)
            stmt = stmt.on_conflict_do_update(
                index_elements=[tracr_tbl.c.operon_id],
                set_=build_update_set(stmt, tracr_tbl),
            )
            conn.execute(stmt)

        # repeat
        repeat_df = tables["repeat"]
        if repeat_df is not None and not repeat_df.empty:
            rows = repeat_df.to_dict(orient="records")
            stmt = pg_insert(repeat_tbl).values(rows)
            stmt = stmt.on_conflict_do_update(
                index_elements=[repeat_tbl.c.operon_id, repeat_tbl.c.array_idx],
                set_={
                    "number_spacers": stmt.excluded.number_spacers,
                    "repeat_seq_length": stmt.excluded.repeat_seq_length,
                    "repeat_sequence": stmt.excluded.repeat_sequence,
                },
            )
            conn.execute(stmt)

        # spacer
        spacer_df = tables["spacer"]
        if spacer_df is not None and not spacer_df.empty:
            rows = spacer_df.to_dict(orient="records")
            stmt = pg_insert(spacer_tbl).values(rows)
            stmt = stmt.on_conflict_do_update(
                index_elements=[
                    spacer_tbl.c.operon_id,
                    spacer_tbl.c.array_idx,
                    spacer_tbl.c.spacer_position,
                ],
                set_={
                    "spacer_length": stmt.excluded.spacer_length,
                    "spacer_sequence": stmt.excluded.spacer_sequence,
                },
            )
            conn.execute(stmt)

        # cas
        cas_df = tables["cas"]
        if cas_df is not None and not cas_df.empty:
            rows = cas_df.to_dict(orient="records")
            stmt = pg_insert(cas_tbl).values(rows)
            stmt = stmt.on_conflict_do_update(
                index_elements=[
                    cas_tbl.c.operon_id,
                    cas_tbl.c.gene_name,
                    cas_tbl.c.hmm_name,
                    cas_tbl.c.evalue,
                    cas_tbl.c.score,
                ],
                set_={
                    "evalue": stmt.excluded.evalue,
                    "score": stmt.excluded.score,
                    "truncated": stmt.excluded.truncated,
                    "length": stmt.excluded.length,
                    "protein": stmt.excluded.protein,
                },
            )
            conn.execute(stmt)

    print(f"upserted; operon.id={summary_df['operon_id'].tolist()}. Children loaded.")     
