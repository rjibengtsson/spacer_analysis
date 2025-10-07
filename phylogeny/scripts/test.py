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



def assembly_type_counter(infile: str):
    """
    Reads a JSON Lines file and prints each JSON object.

    Args:
        infile (str): Path to the input JSON Lines file.
    """

    counter_dict = {}

    with open(infile, 'rb') as inf:
        for line in inf:
            if not line.strip():
                continue
            obj = orjson.loads(line)
            at = obj.get("metadata", {}).get("assembly_type")
            if at:
                counter_dict[at] = counter_dict.get(at, 0) + 1

        for at_type, count in counter_dict.items():
            print(f"{at_type}: {count}")
            



def reverse_complement(seq: str) -> str:
    """
    Returns the reverse complement of a DNA sequence.

    Args:
        seq (str): Input DNA sequence.

    Returns:
        str: Reverse complement of the input sequence.
    """
    complement = str.maketrans('ATCGatcg', 'TAGCtagc')
    return seq.translate(complement)[::-1]



def get_gRNA_seq(spacer_seq: str) -> str:
    """
    Transcribe a DNA spacer sequence into gRNA sequence.
    
    Args:
        spacer_seq (str): DNA spacer sequence (e.g., "ATGCTTACG")
        
    Returns:
        str: Transcribed gRNA sequence (e.g., "AUGCUUACG") in 5' to  3' direction.
    """

    dna = Seq(spacer_seq.upper())
    reverse_complement = dna.reverse_complement()
    gRNA = reverse_complement.transcribe()

    return str(gRNA)



def find_matches(consensus_list: t.List[str], long_seq: str) -> bool:
    rev_comp = reverse_complement(long_seq)

    for spacer in consensus_list:
        if spacer in long_seq or spacer in rev_comp:
            return spacer


def generate_spacer_consensus(spacer_list: t.List[str]) -> t.List[str]:
    u = {s for s in spacer_list if s and len(s) <= 36 }
    return list({min(s, reverse_complement(s)) for s in u})


def correct_spacer_sequence(spacer_list: t.List[str], consensus_list: t.List[str]) -> t.List[str]:
    corrected_spacers = []
    
    for spacer in spacer_list:
        if spacer in consensus_list:
            corrected_spacers.append(spacer)
        elif reverse_complement(spacer) in consensus_list:
            corrected_spacers.append(reverse_complement(spacer))
        elif match := find_matches(consensus_list, spacer):
            corrected_spacers.append(match)
        else:
            corrected_spacers.append(spacer)

    return [s for s in corrected_spacers if s != '']


def extract_from_cdhit_cluster(infile: str) -> t.List[str]:

    # format clstr file from cd-hit into dictionary
    cluster_dict = {}

    with open(infile, 'r') as inf:
        for line in inf:
            line = line.strip()
            if line.startswith('>'):
                cluster_num = line.split()[-1]
            elif (line.endswith('*')) or (line.endswith('%')):
                accn_list = cluster_dict.get(cluster_num, [])
                fields = re.split(r'[,\t ]+', line)
                aa_len = fields[1].replace('aa,', '')
                accn = fields[2].lstrip(' >').replace('...', '')
                accn_list.append(accn)
                cluster_dict[cluster_num] = accn_list

    return cluster_dict


def retrieve_guides(accn_list: t.List[str], json_file:str) -> t.List[str]:

    operon_ids = []
    
    for accn in accn_list:
        if accn == "GCA_000833995.1_SAMN03284264":
            continue
        else:
            assembly_id = accn.split('@')[0]
            operonnum = accn.split('@')[1].split('_')[0]
            operon_id = f"{assembly_id}@{operonnum}"
            operon_ids.append(operon_id)

    spacer_list = []
    spacer_dict = {}
    
    with open(json_file, 'rb') as inf:
        for line in inf:
            line = line.strip()
            if not line:
                continue
            obj = json.loads(line)
            if obj.get("operon_id") in operon_ids:
                spacers = obj.get("crispr", {})[0].get("crispr_spacers")
                spacer_list.extend(spacers)
                # generate crispr array dictionary
                crispr_dict = {}
                crispr_dict[obj.get("metadata", {}).get("biosample_id")] = spacers
                spacer_dict[obj.get("operon_id")] = crispr_dict

    # correct for length variation and orientation of the spacers
    consensus_list = generate_spacer_consensus(spacer_list=spacer_list)

    corrected_spacers_list = []
    for operon_id, n_spacers in spacer_dict.items():
        for biosample_id, spacers in n_spacers.items():
            corrected_spacers = correct_spacer_sequence(spacers, consensus_list)
            corrected_spacers_list.extend(corrected_spacers)

    filtered_spacer_list = [spacer for spacer in corrected_spacers_list if len(spacer) == 29]
    gRNA_seqs = list(set([get_gRNA_seq(seq) for seq in filtered_spacer_list]))

    return gRNA_seqs


def count_nucleotide(spacers: list, position: int) -> t.Dict[str, int]:
    """
    Calculate the nucleotide counts at a specific position across a list of spacer sequences.
    
    Args:
        spacers (list): List of spacer sequences.
        position (int): The position to check (0-based index).
    
    Returns:
        Dict[str, int]: A dictionary with nucleotide counts at the specified position.
    """
    
    counts = {'A': 0, 'C': 0, 'G': 0, 'U': 0}
    
    for seq in spacers:
        if len(seq) > position-1:
            nucleotide = seq[position-1]
            if nucleotide in counts:
                counts[nucleotide] += 1
    
    return counts


def count_nucleotide_freq(spacers: list, position: int) -> t.Dict[str, float]:
    """
    Calculate the nucleotide frequencies at a specific position across a list of spacer sequences.
    
    Args:
        spacers (list): List of spacer sequences.
        position (int): The position to check (0-based index).
    
    Returns:
        Dict[str, float]: A dictionary with nucleotide frequencies at the specified position.
    """
    
    counts = count_nucleotide(spacers, position)
    total = sum(counts.values())
    
    if total == 0:
        return {nt: 0.0 for nt in counts}
    
    freqs = {nt: round(count / total, 4) for nt, count in counts.items()}
    
    return freqs



def position_nucleotide_freq(gRNA_seqs: t.List[str], nucleotide: str) -> t.Dict[int, float]:

    position_nucleotide_freq = {}
    for i in range(1, 31):
        freqs = count_nucleotide_freq(gRNA_seqs, position=i).get(nucleotide)
        position_nucleotide_freq[i] = freqs
    
    return position_nucleotide_freq



def main():
    cluster_file = "/spacers_db/crispr-cas-atlas/cd-hit_all/db_90.clstr"
    json_file = "/spacers_db/crispr-cas-atlas/crispr-cas-atlas-v1.0-VI.jsonl"
    cluster_dict = extract_from_cdhit_cluster(cluster_file)
    
    nucleotide = 'U'

    colour_dict = {"A": "#008000",
                "U": "#FF0000",
                "G": "#ffa100",
                "C": "#0000ff"}

    # keep only clusters with >= 10 members
    new_dict = {k: v for k, v in cluster_dict.items() if len(v) >= 10}
    
    results_dict = {}

    # # collect rows first (faster than repeated concat)
    rows = []
    for key, val in new_dict.items():
        gRNA_seqs = retrieve_guides(val, json_file)
        if not gRNA_seqs:
            continue
        new_row = position_nucleotide_freq(gRNA_seqs, nucleotide=nucleotide)  # should return 30 values
        print(new_row)
        exit(0)
        # rows.append(new_row)


    # # build DataFrame with columns 1..30 (as strings → discrete)
    # pos_order = [i for i in range(1, 31)]
    # df = pd.DataFrame(columns=pos_order)

    # for row in rows:
    #     new_row = pd.DataFrame([row])
    #     df = pd.concat([df, new_row], ignore_index=True)

    # # melt to long form
    # df_melt = df.melt(var_name='position', value_name='frequency')
    # # enforce ordered categorical for correct x ordering
    # df_melt['position'] = pd.Categorical(df_melt['position'], categories=pos_order, ordered=True)


    # # plot
    # p_box = (
    #     ggplot(df_melt, aes(x='position', y='frequency')) +
    #     geom_boxplot(fill=colour_dict.get(nucleotide), alpha=0.4) +
    #     geom_jitter(alpha=0.3, width=0.15, size=1) +
    #     scale_x_discrete(limits=pos_order) +
    #     scale_y_continuous(limits=(0, 1), breaks=np.arange(0, 1.1, 0.1)) +
    #     # geom_hline(yintercept=mean_dict.get(nucleotide), linetype='dashed', color='red', size=1) +
    #     theme_bw() +
    #     labs(
    #         title='',
    #         x='Guide position',
    #         y=f'{nucleotide} frequency'
    #     ) +
    #     theme(axis_text_x=element_text(rotation=90, ha='right'))
    #     )  

    
    # fig = p_box.draw()
    # plt.show()
    # p_box.save(f"{nucleotide}_frequency_plot.png", width=8, height=5, dpi=300)


    # for cluster_id, pos_freq in results_dict.items():
    #     print(cluster_id, pos_freq.get(1))


if __name__ == "__main__":
    main()