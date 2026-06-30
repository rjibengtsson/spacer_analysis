import os, sys
import uuid
import re
import subprocess
import typing as t
from typing import Optional
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import sys, os
import pandas as pd
import json
from pathlib import Path
from dataclasses import dataclass
from concurrent.futures import ProcessPoolExecutor


# Add the directory containing your other modules
target_dir = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(target_dir))

import modules.analysis_utils as analysis_utils
import modules.database_utils as db_utils
from modules.blast_utils import PhageBlast
from modules.analysis_utils import PhageElement

seq1 = "TTATGAGTAAAAATAATCCTACTAAAGCTG"
seq2 = "TATAATGGAGGACATTATTATGAAAATCCT"


print(seq1)
print(seq2)
