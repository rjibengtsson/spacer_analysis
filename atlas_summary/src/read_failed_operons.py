from turtle import pd
import pandas as pd
import typing as t
import load_data2sql
import upload_parallel
import sys, os
import ast


def get_failed_operons(file_path):
    operons_list = []

    with open(file_path, "r") as f:
        failed_operons = [line.strip() for line in f][1]
        op_str = failed_operons.lstrip("operon_ids=")
        data = ast.literal_eval(op_str)
        operons_list.extend(data)
    return operons_list


def process_record_chunk(records: list) -> dict[str, pd.DataFrame]:
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

    result = {}
    for table_name, df_list in chunk_tables.items():
        result[table_name] = (
            pd.concat(df_list, ignore_index=True)
            if df_list else pd.DataFrame()
        )
    return result



def filter_records_by_operon_id(input_file: str, target_operon_ids: set[str]):
    for record in load_data2sql.iter_jsonl(input_file):
        if record.get("operon_id") in target_operon_ids:
            yield record



engine = load_data2sql.get_connection()

def main():
    INPUT_FILE = "/spacers_db/crispr-cas-atlas/crispr-cas-atlas-v1.0.jsonl"
    failed_operons_file = "failed_operons.list"
    failed_operons = get_failed_operons(failed_operons_file)


    json_file_path = INPUT_FILE

    filtered_records = filter_records_by_operon_id(INPUT_FILE, failed_operons)
    chunk_iter = upload_parallel.chunked_iterator(filtered_records, 500)

    for chunk in chunk_iter:
        print(len(chunk))
        # tables = process_record_chunk(chunk)
        # print(tables['cas'])
        # tables['cas'].to_csv("cas_table.tsv", sep ="\t", index=False, header=True, mode='a', encoding='utf-8')
    # record_iter = load_data2sql.iter_jsonl(json_file_path)

    # matching_df = process_record_chunk(record_iter, set(failed_operons))
    # print(matching_df)




if __name__ == "__main__":
    main()