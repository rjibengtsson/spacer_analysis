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
from pathlib import Path

# Add the directory containing your other modules
target_dir = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(target_dir))

import modules.atlas_handle as ah



def main():
    fasta_file = "/spacers_db/crispr-cas-atlas/cas13b_all.faa.filtered.faa"
    json_file = "/spacers_db/crispr-cas-atlas/crispr-cas-atlas-v1.0-VI.jsonl"
    acession_list = "phylogeny/pspCas13b_cluster.list"

    accessions = [line.strip() for line in open(acession_list)]

    spacer_list, spacer_dict = ah.retrieve_spacers(accn_list=accessions, json_file=json_file)
        
    consensus_list = ah.generate_spacer_consensus(spacer_list=spacer_list)
    final_spacers_list = []
    for operon_id, n_spacers in spacer_dict.items():
        for biosample_id, spacers in n_spacers.items():
            corrected_spacers = ah.correct_spacer_sequence(spacers, consensus_list)
            final_spacers_list.extend(corrected_spacers)

    

if __name__ == "__main__":
    main()