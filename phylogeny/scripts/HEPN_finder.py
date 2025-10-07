from Bio import SeqIO
import re, sys


def find_hepn_domains(fasta_file):
    
    pattern = r"R.{4}H"
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq)
        matches = [(m.start() + 1, m.end(), m.group()) for m in re.finditer(pattern, seq)]
        print(f"{matches}")


find_hepn_domains(fasta_file=sys.argv[1])