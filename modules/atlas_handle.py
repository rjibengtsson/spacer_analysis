import ijson, json
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import orjson
import os, sys
import re
import typing as t
from pathlib import Path



########################### Functions ####################################################


def subtype_counter(infile: str):
    """
    Checks and prints unique subtype values from a large JSON file.

    Args:
        infile (str): Path to the input JSON file.
    """
    counter_dict = {}

    with open(infile, 'rb') as inf:
        # Stream each top-level object in the array
        for obj in ijson.items(inf, 'item', use_float=True):
            # Extract the nested subtype key
            subtype = obj.get("summary", {}).get("subtype")
            counter_dict[subtype] = counter_dict.get(subtype, 0) + 1
    

    for subtype, count in counter_dict.items():
        print(f"{subtype}: {count}")

    return counter_dict



def filter_json_by_subtype(infile: str, outfile: str, subtype_value="VI-B"):
    """
    Filters a large JSON file to include only entries with a specific subtype.

    Args:
        infile (str): Path to the input JSON file.
        outfile (str): Path to the output JSON file.
        subtype_value (str): The subtype value to filter by. Default is "VI-B".
    """
    with open(infile, 'rb') as inf, open(outfile, 'wb') as outf:
        # Stream each top-level object in the array
        for obj in ijson.items(inf, 'item', use_float=True):
            # Check the nested subtype key
            if obj.get("summary", {}).get("subtype") == subtype_value:
                outf.write(orjson.dumps(obj) + b'\n')




def search_subtypes(infile: str, target_subtypes: str):
    """
    Searches for and prints JSON objects with specific subtype values.

    Args:
        infile (str): Path to the input JSON file.
    """

    with open(infile, 'rb') as inf:
        for obj in ijson.items(inf, 'item', use_float=True):
            subtype = obj.get("metadata", {}).get("biosample_id")
            if subtype == target_subtypes:
                print(orjson.dumps(obj, option=orjson.OPT_INDENT_2).decode("utf-8"))






def gene_counter(infile: str):
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
            for cas in obj.get('cas', []):
                # if cas['evalue'] < 1e-7:
                counter_dict[cas['gene_name']] = counter_dict.get(cas['gene_name'], 0) + 1

        for cas_type, count in counter_dict.items():
            print(f"{cas_type}: {count}")


def find_object(infile: str, search_value):
    """
    Finds and prints JSON objects in a file that contain a specific key-value pair.

    Args:
        infile (str): Path to the input JSON file.
        search_key (str): The key to search for.
        search_value: The value to match.
    """


    with open(infile, 'rb') as inf:
        for line in inf:
            line = line.strip()
            if not line:
                continue
            obj = json.loads(line)
            if obj.get("metadata", {}).get("biosample_id") == search_value:
                print(orjson.dumps(obj, option=orjson.OPT_INDENT_2).decode("utf-8"))
                # print(obj.get("metadata", {}).get("sample_name"))


def get_protein_sequences(infile: str, outfile: str, assembly_type="Genome", gene_name="Cas13b"):
    """
    Extracts protein sequences from a JSON Lines file and writes them to a FASTA file.

    Args:
        infile (str): Path to the input JSON Lines file.
        outfile (str): Path to the output FASTA file.
    """
    
    selection_criteria = {
        "gene_name": gene_name,
        "truncated": "00"
    }
    
    with open(infile, 'rb') as inf, open(outfile, 'w') as outf:
        for line in inf:
            line = line.strip()
            if not line:
                continue
            obj = json.loads(line)
            # if obj.get("metadata", {}).get("assembly_type") == assembly_type:
            for cas in obj.get('cas', []):
                if cas['gene_name'] == selection_criteria['gene_name'] and cas['truncated'] == selection_criteria['truncated']:
                    header = f">{obj['operon_id']}_{obj['metadata']['biosample_id']}"
                    sequence = cas['protein']
                    outf.write(f"{header}\n{sequence}\n")



def P5_125_spacer_seqs():
    return ["TCGTTTTGCCCGAGCCAATCAGCAGAAGC",
            "CCGTTCCGTCGTGGTATTTATAGTTATTCT",
            "AAGTACTAGCACAGAAGAGAATGTCATTTT",
            "GTATGCCATTTGCGATAAATCTCTCGGTG",
            "TTCTAATAACAATCATTAATCCTTAAAGCA",
            "TCACAATGGTTGTGGTATAGTGCAAGTCT",
            "TGTGGCAGACATCCTAAAAAAGCATGGTAT",
            "AAATAGGATAACATCTCCTTAAAAGAAAGT",
            "TATATCTTCATCGAAGTTAACGGTAATGCT",
            "CCTTGAGCAAGGCGTACCGATGAGATAGC",
            "CTCATACTTTTATGAAAAATATTTCAACTG",
            "GGGTATTCTTGAGCGTGAAAAAGTAAGGGA"]



def extract_from_cdhit_cluster(infile: str, fasta_file: str,ref_accn="GCA_000833995.1_SAMN03284264") -> t.List[str]:

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

    # extract all sequences from the cluster containing ref_accn
    cluster_key = next((k for k, v in cluster_dict.items() if ref_accn in v), None)
    target_accns = cluster_dict.get(cluster_key, [])

    return target_accns


def write_fasta_sequences(fasta_file: str, target_accns: t.List[str], outfile: str):

    with open(outfile, 'w') as outf:
        for accn in target_accns:
            for record in SeqIO.parse(fasta_file, "fasta"):
                if record.id == accn:
                    SeqIO.write(record, outf, "fasta")



def retrieve_spacers(accn_list: t.List[str], json_file:str) -> t.List[str]:

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
                

    return spacer_list, spacer_dict


def retrieve_repeat_sequences(accn_list: t.List[str], json_file:str) -> dict:

    operon_ids = []
    
    for accn in accn_list:
        if accn == "GCA_000833995.1_SAMN03284264":
            continue
        else:
            assembly_id = accn.split('@')[0]
            operonnum = accn.split('@')[1].split('_')[0]
            operon_id = f"{assembly_id}@{operonnum}"
            operon_ids.append(operon_id)

    repeat_dict = {}
    
    with open(json_file, 'rb') as inf:
        for line in inf:
            line = line.strip()
            if not line:
                continue
            obj = json.loads(line)
            if obj.get("operon_id") in operon_ids:
                repeats = obj.get("crispr", {})[0].get("crispr_repeat")
                repeat_dict[obj.get("operon_id")] = repeats

    return repeat_dict



def reverse_complement(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


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

############# Main script execution ####################################################



def _retrieve_spacers():
    cluster_file = "/spacers_db/crispr-cas-atlas/cd-hit_all/db_98.clstr"
    fasta_file = "/spacers_db/crispr-cas-atlas/cas13b_all.faa.filtered.faa"
    json_file = "/spacers_db/crispr-cas-atlas/crispr-cas-atlas-v1.0-VI.jsonl"
    accession_list = "/spacers_db/crispr-cas-atlas/cluster3.list"
    # ref_accession = "ERZ4561846@1306_SAMEA4765562"
    
    # cluster_accession_list = extract_from_cdhit_cluster(infile=cluster_file,
    #                         fasta_file=fasta_file,
    #                         ref_accn=ref_accession)


    accessions = [line.strip() for line in open(accession_list)]

    spacers_list, spacer_dict = retrieve_spacers(accn_list=accessions, json_file=json_file)

    # consensus_list = generate_spacer_consensus(spacer_list=spacers_list)

    # final_spacers_list = []

    for operon_id, n_spacers in spacer_dict.items():
        for biosample_id, spacers in n_spacers.items():
            print(operon_id)
            print(spacers)
            # corrected_spacers = correct_spacer_sequence(spacers, consensus_list)
            # for s in spacers:
            #     print(s)
            # final_spacers_list.extend(corrected_spacers)

    # for spacer in set(final_spacers_list):
    #     print(spacer)


    # p5_125_spacers = [reverse_complement(s) for s in P5_125_spacer_seqs()]

    # final_spacers_list.extend(p5_125_spacers)

    # for s in set(final_spacers_list):
    #     # print(get_gRNA_seq(s))
    #     print(s)



def _retrieve_repeats():
    cluster_file = "/spacers_db/crispr-cas-atlas/cd-hit_all/db_98.clstr"
    fasta_file = "/spacers_db/crispr-cas-atlas/cas13b_all.faa.filtered.faa"
    json_file = "/spacers_db/crispr-cas-atlas/crispr-cas-atlas-v1.0-VI.jsonl"
    ref_accession = "ERZ4560302@96_SAMEA4766213"

    cluster_accession_list = extract_from_cdhit_cluster(infile=cluster_file,
                            fasta_file=fasta_file,
                            ref_accn=ref_accession)
    
    repeat_dict = retrieve_repeat_sequences(accn_list=cluster_accession_list, json_file=json_file)
    repeat_list = list(repeat_dict.values())
    consensus_list = generate_spacer_consensus(spacer_list=repeat_list)
    
    for repeat in repeat_list:
        if repeat not in consensus_list:
            rev = reverse_complement(repeat)
            print(rev)
        else:
            print(repeat)



def _get_spacer_seqs():
    fasta_file = "/spacers_db/crispr-cas-atlas/cas13b_all.faa.filtered.faa"
    json_file = "/spacers_db/crispr-cas-atlas/crispr-cas-atlas-v1.0-VI.jsonl"
    accession_list = "/home/unimelb.edu.au/rbengtsson/work/spacer_analysis/tmp.list"

    accessions = [line.strip() for line in open(accession_list)]

    spacer_list, spacer_dict = retrieve_spacers(accn_list=accessions, json_file=json_file)

    for operon_id, n_spacers in spacer_dict.items():
        for biosample_id, spacers in n_spacers.items():
            print(operon_id)
            for spacer in spacers:
                print(spacer)
            print("\n")
            # corrected_spacers = ah.correct_spacer_sequence(spacers, consensus_list)
            # final_spacers_list.extend(corrected_spacers)





if __name__ == "__main__":
    _get_spacer_seqs()
    # _retrieve_spacers()
    # _retrieve_repeats()