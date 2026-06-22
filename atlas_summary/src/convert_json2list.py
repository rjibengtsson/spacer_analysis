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



def convert(obj):
    if isinstance(obj, Decimal):
        return float(obj)      # or: return str(obj)
    raise TypeError

def convert_json_to_list():

    atlas_json_file = "/spacers_db/crispr-cas-atlas/crispr-cas-atlas-v1.0.json"
    out_file =  "/spacers_db/crispr-cas-atlas/crispr-cas-atlas-v1.0.jsonl"
    with open(atlas_json_file, 'r') as fin, open(out_file, "wb") as fout:
        for record in ijson.items(fin, "item"):
            fout.write(orjson.dumps(record, default=convert) + b"\n")


if __name__ == "__main__":
    convert_json_to_list()