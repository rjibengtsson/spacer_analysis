import sys
from pathlib import Path
import os
import uuid

# Add the directory containing your other modules
target_dir = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(target_dir))


import modules.database_utils as db_utils
from modules.database_utils import Allthebacteria
import pandas as pd
import modules.genome_utils as genome_utils
from modules.genome_utils import CasContig
import modules.annotation as annotation
from modules.crispr_detecion import CrisprDetection
from modules.analysis_utils import ArrayCandidate


storage = Path("/spacers_db/genomes")
df = db_utils.read_table_from_db(f"{storage}/cas13b_contigs.csv")


for index, row in df.iterrows():
    # check if the cas_gene is not NaN or empty
    if pd.isna(row['cas_gene']) or row['cas_gene'] == "":
        pass
    else:
        # Run data extraction
        bioaccession = row['biosample_accession']
        gbff_file = storage / f"{bioaccession}.gbff"
        gff_file = storage / f"{bioaccession}.gff"
        result_folder = Path(f"/{storage}/{bioaccession}_crispridentify/{bioaccession}")

        candidates = ArrayCandidate.complie_crispr_arrays(result_folder, gbff_file, gff_file, storage)
        filter_candidates = ArrayCandidate.filter_candidates(candidates, avg_dr_len=36, avg_spacer_len=30)

        complete_candidates = []

        for c in filter_candidates:
            print(c.biosample_accn)
            candidate = ArrayCandidate.get_spacer_dr_seq(c, result_folder)
            # print(candidate)
            if candidate is not None:
                complete_candidates.append(candidate)


        # Get the CasContig class instance for the given gbff and gff files
        # This will extract the contig_id and other relevant information
        cls = CasContig.get_cas_info(gbff_file, gff_file)


        array_df = db_utils.generate_array_table(complete_candidates, cls)
        spacer_df = db_utils.generate_spacer_table(complete_candidates)


        db_utils.upload_arraytable_to_sql(array_df, "cas13_bacterial_db", "cas13b_crisprs")
        db_utils.upload_spacertable_to_sql(spacer_df, "cas13_bacterial_db", "spacer_table")

