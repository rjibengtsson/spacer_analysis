from Bio import SeqIO
from Bio.Seq import Seq
import ijson, json
import orjson
import os, sys
import re
import typing as t
import pandas as pd
from IPython.display import display
from plotnine import (ggplot, aes, geom_histogram, geom_density, theme_bw, 
                      labs, theme, element_text, scale_x_continuous, geom_vline, 
                      geom_boxplot, geom_jitter, scale_x_discrete, stat_summary,
                      scale_y_continuous)
import matplotlib.pyplot as plt
import numpy as np
from Bio.SeqRecord import SeqRecord


def get_cluster_dict() -> dict:
    cluster_file = "abacbs2025/data/Cas13b_protein_cd-hit_80i_90c.clstr"
    
    cluster_dict = {}

    with open(cluster_file, 'r') as inf:
        for line in inf:
            line = line.strip()
            if line.startswith('>'):
                cluster_num = line.split()[-1]
            elif (line.endswith('*')) or (line.endswith('%')):
                accn_list = cluster_dict.get(cluster_num, [])
                fields = re.split(r'[,\t ]+', line)
                accn = fields[2].lstrip(' >').replace('...', '')
                accn_list.append(accn)
                cluster_dict[cluster_num] = accn_list

    return cluster_dict    



def extract_cas13b_seqs(operon_ids: t.List[str], json_file: str):
    results = []
    with open(json_file, 'r') as f:
        for line in f:
            record = orjson.loads(line)
            if record['operon_id'] in operon_ids:
                for cas in record.get('cas', []):
                    if cas['gene_name'] == 'Cas13b':
                        cas_protein = cas['protein']
                        results.append(SeqRecord(Seq(cas_protein), id=record['operon_id'], description=""))

    return results


def return_operon_info(operon_id: str) -> dict:
    atlas_json_file = "/spacers_db/crispr-cas-atlas/crispr-cas-atlas-v1.0-VI.jsonl"
    with open(atlas_json_file, 'r') as f:
        for line in f:
            record = orjson.loads(line)
            if record['operon_id'] == operon_id:
                return record
    return {}



def get_spacer_count(operon_id: str, json_file: str) -> int:
    with open(json_file, 'r') as f:
        for line in f:
            record = orjson.loads(line)
            if record['operon_id'] == operon_id:
                spacers = record.get("crispr", {})[0].get("crispr_spacers", [])
                filt_spacers = [s for s in spacers if len(s) == 30]
                return len(filt_spacers)
    return 0



def main():

    cluster_num = "51"

    cluster_dict = get_cluster_dict()    
    cluster_operon_id = cluster_dict.get(cluster_num)

    for operon_id in cluster_operon_id:
        operon_info = return_operon_info(operon_id)
        cas_list = operon_info.get('cas')
        # print(len(cas_list))
        for cas in cas_list:
            print(cas['gene_name'])


    # operon_ids_file = "/home/unimelb.edu.au/rbengtsson/work/spacer_analysis/abacbs2025/data/tmp.list"
    # atlas_json_file = "/spacers_db/crispr-cas-atlas/crispr-cas-atlas-v1.0-VI.jsonl"
    
    # with open(operon_ids_file, 'r') as f:
    #     operon_ids = [line.strip() for line in f.readlines()]
    
    # with open("/home/unimelb.edu.au/rbengtsson/work/spacer_analysis/abacbs2025/data/tmp_30bp.list", 'w') as out_f:
    #     for operon_id in operon_ids:
    #         spacer_count = get_spacer_count(operon_id, atlas_json_file)
    #         out_f.write(f"{operon_id}\t{spacer_count}\n")




    # cas13b_summary_file = "/home/unimelb.edu.au/rbengtsson/work/spacer_analysis/abacbs2025/data/crispr_cas_atlas_Cas13b_summary_edit.csv"

    # out_dir = "/spacers_db/tmp/"

    # species = "Prevotella intermedia"

    # df = pd.read_csv(cas13b_summary_file, header=0)
    # df_filtered = df[df['species'] == species]
    # operon_ids = df_filtered['operon_id'].tolist()
    # cas13b_seqs = extract_cas13b_seqs(operon_ids, atlas_json_file)

    # fasta_outfile = os.path.join(out_dir, f"{species.replace(' ', '_')}_cas13b_seqs.fasta")

    # with open(fasta_outfile, 'w') as f:
    #     SeqIO.write(cas13b_seqs, f, "fasta")


if __name__ == "__main__":
    main()
