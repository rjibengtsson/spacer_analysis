"""
Script to upload data in parallel to PostgreSQL. 
This is used to speed up the upload of large files.
"""

import load_data2sql 
import json
import io
import os
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
from itertools import islice
import typing as t
import pandas as pd
import logging



def chunked_iterator(iterator, chunk_size: int):
    """
    Yield lists of records of size chunk_size.
    """
    iterator = iter(iterator)
    while True:
        chunk = list(islice(iterator, chunk_size))
        if not chunk:
            break
        yield chunk



Tables = t.Dict[str, pd.DataFrame]

def process_record_chunk(records: t.List) -> load_data2sql.Tables:
    """
    Process a chunk of operon records and return one dataframe per table.
    """
    chunk_tables = {
        "summary": [],
        "metadata": [],
        "tracr": [],
        "repeat": [],
        "spacer": [],
        "cas": [],
    }

    for record in records:
        record_tables = load_data2sql.extract_record_info(record)
        
        for table_name, df in record_tables.items():
            if df is not None and not df.empty:
                chunk_tables[table_name].append(df)

    # concatenate each table into one dataframe for the chunk
    result = {}
    for table_name, df_list in chunk_tables.items():
        if df_list:
            result[table_name] = pd.concat(df_list, ignore_index=True)
        else:
            result[table_name] = pd.DataFrame()  # empty dataframe if no data

    return result


engine = load_data2sql.get_connection()


logging.basicConfig(
    filename="failed_operons.log",
    level=logging.ERROR,
    format="%(asctime)s | %(levelname)s | %(message)s"
)

logger = logging.getLogger(__name__)



def main():

    INPUT_FILE = "/spacers_db/crispr-cas-atlas/crispr-cas-atlas-v1.0.jsonl"
    CHUNK_SIZE = 500
    MAX_WORKERS = max(1, os.cpu_count() - 1)
    MAX_INFLIGHT_TASKS = MAX_WORKERS * 2  # limit the number of tasks in flight to avoid overwhelming the system



    json_file_path = INPUT_FILE
    record_iter = load_data2sql.iter_jsonl(json_file_path)
    chunk_iter = chunked_iterator(record_iter, CHUNK_SIZE)


    sql_tables = load_data2sql.get_reflected_tables(engine)

    # keep track of the futures (running tasks) so we can wait for them to finish
    futures = set()

    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        # Fill the initial batch of tasks
        for _ in range(MAX_INFLIGHT_TASKS):
            try:
                chunk = next(chunk_iter)
                future = executor.submit(process_record_chunk, chunk)
                futures.add(future)
            except StopIteration:
                break  # no more chunks to process
        
        while futures:
            done, futures = wait(futures, return_when=FIRST_COMPLETED)

            for future in done:

                try:
                    tables = future.result()
                except Exception as e:
                    logger.error(
                        f"Worker failed while processing chunk | error={str(e)}",
                        exc_info=True
                    )
                    continue  # skip to the next completed future


                try:
                    load_data2sql.upload_chunk_to_sql(engine, sql_tables, tables)

                except Exception as e:
                    # extract operon_ids from the summary table
                    try:
                        operon_ids = tables.get("summary", pd.DataFrame())["operon_id"].tolist()
                    except Exception as e:
                        operon_ids = ["UNKNOWN"]
                    
                    logger.error(
                        f"operon_ids={operon_ids} | error={str(e)}",
                        exc_info=True
                    )

                # Refill the task queue with new chunks as they complete
                try:
                    chunk = next(chunk_iter)
                    future = executor.submit(process_record_chunk, chunk)
                    futures.add(future)
                except StopIteration:
                    pass  # no more chunks to process


    # for chunk in chunk_iter:
    #     tables = process_record_chunk(chunk)
    #     cas_table = tables['cas']
    #     print(cas_table[cas_table.duplicated(subset=['operon_id', 'gene_name', 'hmm_name'], keep=False)])
    #     exit(0)
    #     # load_data2sql.upload_tables(tables)





if __name__ == "__main__":
    main()