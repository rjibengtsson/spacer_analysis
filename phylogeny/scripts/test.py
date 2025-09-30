from Bio import SeqIO
import ijson, json
import orjson
import os, sys



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




# assembly_type_counter(infile=sys.argv[1])
print(reverse_complement("TAGAATTAGTTTTAATGTCAACGAGAATA".upper()))





