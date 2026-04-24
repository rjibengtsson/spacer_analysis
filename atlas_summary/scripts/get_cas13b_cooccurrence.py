from Bio import SeqIO
from Bio.Seq import Seq
import ijson, json
import orjson
import os, sys
import re
import typing as t
from typing import Optional, Dict, Any, List
from dataclasses import dataclass
import pandas as pd
import numpy as np
from decimal import Decimal
from IPython.display import display
from plotnine import (ggplot, aes, geom_histogram, geom_density, theme_bw, 
                      labs, theme, element_text, scale_x_continuous, geom_vline, 
                      geom_boxplot, geom_jitter, scale_x_discrete, stat_summary,
                      scale_y_continuous, scale_fill_manual, scale_color_manual, geom_hline,
                      element_blank)
import matplotlib.pyplot as plt


def convert(obj):
    if isinstance(obj, Decimal):
        return float(obj)      # or: return str(obj)
    raise TypeError


def generate_genome_subsets():

    atlas_json_file = "/spacers_db/crispr-cas-atlas/crispr-cas-atlas-v1.0.json"
    out_file =  "/spacers_db/crispr-cas-atlas/crispr-cas-atlas_genomes.jsonl"
    with open(atlas_json_file, 'r') as fin, open(out_file, "wb") as fout:
        for record in ijson.items(fin, "item"):
            assembly_type = record["metadata"]["assembly_type"]
            if assembly_type == "Genome":
                fout.write(orjson.dumps(record, default=convert))
                fout.write(b"\n")


@dataclass
class Operon:

    operon_id: Optional[str] = None
    summary: Optional[Dict[str, Any]] = None
    crispr: Optional[List[Dict[str, Any]]] = None
    cas: Optional[List] = None


    def retrieve_operon_info(biosample_id):
        operons_dict = []
        atlas_json_file = "/spacers_db/crispr-cas-atlas/crispr-cas-atlas_genomes.jsonl"
        with open(atlas_json_file, 'r') as f:
            for line in f:
                record = orjson.loads(line)
                if record["metadata"]["biosample_id"] == biosample_id:
                    operon = Operon(
                        operon_id=record["operon_id"],
                        summary=record["summary"],
                        crispr=record["crispr"],
                        cas=record["cas"]
                    )
                    operons_dict.append(operon.__dict__)
        return operons_dict


def retrieve_biosample_id():
    atlas_json_file = "/spacers_db/crispr-cas-atlas/crispr-cas-atlas_genomes.jsonl"
    out_file = "/spacers_db/crispr-cas-atlas/testing/test_genomes_operons_info.jsonl"

    genome_biosample_ids = []

    with open(atlas_json_file, 'r') as f, open(out_file, "w") as fout:
        for i, line in enumerate(f):  # Just process first 10 records for testing
            if i == 50:
                break
            record = orjson.loads(line)
            biosample_id = record["metadata"]["biosample_id"]
            if biosample_id not in genome_biosample_ids:
                operons_list = Operon.retrieve_operon_info(biosample_id)
                data = {"biosample_id": biosample_id,
                        "metadata": record["metadata"],
                        "operons": operons_list}
                json.dump(data, fout, indent=4)
                genome_biosample_ids.append(biosample_id)

    
    # return genome_biosample_ids


def cas_types() -> List[str]:
    return ["Cas3", "Cas9", "Cas10", "Cas12", "Cas13"]


def iter_json_objects(path):
    decoder = json.JSONDecoder()
    buffer = ""

    with open(path, "r") as f:
        for line in f:
            buffer += line
            buffer = buffer.lstrip()
            while buffer:
                try:
                    obj, idx = decoder.raw_decode(buffer)
                    yield obj
                    buffer = buffer[idx:].lstrip()
                except json.JSONDecodeError:
                    # Not enough data yet to decode a full object; read more lines
                    break


def operon_length_distribution():
    atlas_json_file = "/spacers_db/crispr-cas-atlas/testing/test_genomes_operons_info.jsonl"
    out_csv_file = "atlas_summary/testing_atlas_operon_length_distribution.csv"

    df = pd.DataFrame(columns=["operon_id", "type", "subtype", "length"])

    for record in iter_json_objects(atlas_json_file):
        biosample_id = record.get("biosample_id", "")
        all_genes = []
        operons = record.get("operons", [])
        for operon in operons:
            if operon.get("summary", {}).get("subtype", "") is None:
                operon_type = "Unknown"
            else:
                operon_type = operon.get("summary", {}).get("subtype", "").split("-")[0]
            operon_length = operon.get("summary", {}).get("operon_length", 0)
            new_row = {"operon_id": operon.get("operon_id", ""),
                       "type": operon_type,
                       "subtype": operon.get("summary", {}).get("subtype", ""),
                       "length": operon_length}
            df = pd.concat([df, pd.DataFrame([new_row])], ignore_index=True)
    
    df.to_csv(out_csv_file, index=False)



def filter_operons(cas_type: str, operon_lengths: int):
    length_distribution_file = "atlas_summary/atlas_operon_length_distribution.csv"

    # get length distribution info
    df_length = pd.read_csv(length_distribution_file, header=0)
    stats = df_length.groupby('type')['length'].describe()

    stats_index = stats.loc[cas_type]
    mean_length = stats_index['mean']
    std_length = stats_index['std']
    lower_bound = mean_length - 2 * std_length
    upper_bound = mean_length + 2 * std_length

    if operon_lengths >= lower_bound and operon_lengths <= upper_bound:
        return 1
    else:
        return 0




def get_cas_types():

        atlas_json_file = "/spacers_db/crispr-cas-atlas/genomes_operons_info.jsonl"

    
        result_dict = {}
        cas_types_list = cas_types()

        for record in iter_json_objects(atlas_json_file):
            biosample_id = record.get("biosample_id", "")
            all_genes = []
            operons = record.get("operons", [])

            for operon in operons:
                if not operon.get('summary', {}).get('subtype'):
                    continue
                else:
                    cas_type = operon.get('summary', {}).get('subtype').split('-')[0]
                    operon_length = operon.get('summary', {}).get('operon_length')

                    # if filter_operons(cas_type, operon_length) == 1:
                    operon_genes = [cas_gene['gene_name'] for cas_gene in operon.get('cas', [])]
                    all_genes.extend(operon_genes)
        
                    matches = [
                        item for item in all_genes 
                        if any(key in item for key in cas_types_list)
                    ]
                    
                    if len(matches) > 0:
                        # unique_matches = set(matches)
                        result_dict[biosample_id] = matches
                        
        return result_dict


def generate_upset_input(cas_cooccurrence_dict: Dict[str, List[str]]):

    upset_input = "atlas_summary/upset_input.csv"

    cas_types_list = list(cas_cooccurrence_dict.values())
    flat = [item for sublist in cas_types_list for item in sublist]
    cas_list = list(set(flat))   
    df_header = "biosample_id," + ",".join(cas_list) 

    df = pd.DataFrame(columns=df_header.split(","))

    for biosample_id, cas_genes in cas_cooccurrence_dict.items():
        row = {"biosample_id": biosample_id, **{cas_type: (1 if cas_type in cas_genes else 0) for cas_type in cas_list}}
        df = pd.concat([df, pd.DataFrame([row])], ignore_index=True)
    
    df.to_csv(upset_input, index=False)
    


def main():
    cas_list = get_cas_types()
    generate_upset_input(cas_list)
    
    # read_operon_stats()
    # retrieve_biosample_id()
    # operon_length_distribution()



if __name__ == "__main__":
    main()