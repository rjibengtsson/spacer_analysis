import sys
from pathlib import Path
import unittest
from collections import Counter
from Bio.Seq import Seq

# Add the directory containing your other modules
target_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(target_dir))

import modules.atlas_handle as ah

class TestAtlasHandle(unittest.TestCase):

    def setUp(self):
        self.accn_list = ["ERZ929346@124_SAMEA17979418",
                          "CosteaPI_2017@14178_None",
                          "GCA_000833975.1@4_SAMN03284263",
                          "GCA_000833995.1_SAMN03284264"]
        
        spacer_test_file = "/home/unimelb.edu.au/rbengtsson/work/spacer_analysis/modules/test_files/concensus_spacer_test.list"
        with open(spacer_test_file, 'r') as f:
            self.spacer_list = [line.strip() for line in f]

        self.json_file = "/spacers_db/crispr-cas-atlas/crispr-cas-atlas-v1.0-VI.jsonl"

        self.consensus_list = ["TAGAATTAGTTTTAATGTCAACGAGAATA",
                    "ACTACACCACAAGTATATATTTATCTTTTT",
                    "TTTACTTTATAGTGATAAACTTTAAGAACT",
                    "CTGTAATTCTTACTAGAGCATACTGCTCTT",
                    "TAAACTCAAAAGTGTTTCCTTTGGTGTCTG",
                    "CACACAACGTAAGTAGACCATTGAACATAG",
                    "GTTCGTTGTAGTCTAACGGTGCTTCTTTCT",
                    "AACATTCCATTTAGAGATATCTCCATTGAA",
                    "GTTGAAATTGTTCATAATTTTACACGTTTT"]


    def tearDown(self):
        pass

    def test_retrieve_spacers(self):
        spacer_list_result = ["ACAATTGGTAACACGTTCATAATGTCACGA","CGTGCGATAATTCCAGAAGACATGCAGGAT","CTTCTTCCATCAAATCAATTCTCGTTTCCA",
                              "GTTGAAATTGTTCATAATTTTACACGTTTT", "AAAACGTGTAAAATTATGAACAATTTCAAC","TGGAAACGAGAATTGATTTGATGGAAGAAG",
                              "CAATAATAATAAACCTTATACTCTTAACTA","GTGACAAAAACATCTTAGCTATCAAGTCCA","ATCCTGCATGTCTTCTGGAATTATCGCACG",
                              "TCGTGACATTATGAACGTGTTACCAATTGT","TTATAATCTTGTTAGTAGAGAAAATTCTAC", "TGAAGATGTAAATCTTATAGGTGAAAAAAT",
                              "GCGTTATGCCTATAAAAACATGTTCAACCT","GTTACAATCATGGAAAAGAAATCATTGGTG","GAACAGCAAAGAGAGGTATGTTCAAAGAGG",
                              "AAAACGTACAAATGGTATACGCACTCTGTG","GATTCAGAGTATGAGAATCTGCAGAAAGCA","AATTTATCATAAACTTATTAATAACAATTT",
                              "AAAAACAAGCAGAGTTAATAATCGATGTAT","TATTAGTGGTATTGCATTCTTCGCATTTTA","GATAACCGAAGTCGGCAAGGTGTTCATGAC",
                              "TTGAAACATTGTAGGTGTATCTCTTGCATC"]
        spacer_dict_result = {"ERZ929346@124": {"SAMEA17979418" : ["ACAATTGGTAACACGTTCATAATGTCACGA","CGTGCGATAATTCCAGAAGACATGCAGGAT","CTTCTTCCATCAAATCAATTCTCGTTTCCA","GTTGAAATTGTTCATAATTTTACACGTTTT"]},
                              "CosteaPI_2017@14178": {None: ["AAAACGTGTAAAATTATGAACAATTTCAAC","TGGAAACGAGAATTGATTTGATGGAAGAAG","CAATAATAATAAACCTTATACTCTTAACTA","GTGACAAAAACATCTTAGCTATCAAGTCCA","ATCCTGCATGTCTTCTGGAATTATCGCACG","TCGTGACATTATGAACGTGTTACCAATTGT","TTATAATCTTGTTAGTAGAGAAAATTCTAC"]},
                              "GCA_000833975.1@4": {"SAMN03284263": ["TGAAGATGTAAATCTTATAGGTGAAAAAAT","GCGTTATGCCTATAAAAACATGTTCAACCT","GTTACAATCATGGAAAAGAAATCATTGGTG","GAACAGCAAAGAGAGGTATGTTCAAAGAGG","AAAACGTACAAATGGTATACGCACTCTGTG","GATTCAGAGTATGAGAATCTGCAGAAAGCA","AATTTATCATAAACTTATTAATAACAATTT","AAAAACAAGCAGAGTTAATAATCGATGTAT","TATTAGTGGTATTGCATTCTTCGCATTTTA","GATAACCGAAGTCGGCAAGGTGTTCATGAC","TTGAAACATTGTAGGTGTATCTCTTGCATC"]}}
        
        spacer_list, spacer_dict = ah.retrieve_spacers(accn_list=self.accn_list, json_file=self.json_file)

        self.assertEqual(Counter(spacer_list), Counter(spacer_list_result))
        for accn in spacer_dict_result:
            for biosample in spacer_dict_result[accn]:
                self.assertEqual(
                    Counter(spacer_dict[accn][biosample]),
                    Counter(spacer_dict_result[accn][biosample])
                )        


    def test_generate_spacer_consensus(self):
        generate_spacer_consensus_result = ah.generate_spacer_consensus(self.spacer_list)
        self.assertEqual(Counter(generate_spacer_consensus_result), Counter(self.consensus_list))


    def test_correct_spacer_sequence(self):

        result = ah.correct_spacer_sequence(self.spacer_list, self.consensus_list)
        expected_result_file = "/home/unimelb.edu.au/rbengtsson/work/spacer_analysis/modules/test_files/corrected_spacer.list"
        with open(expected_result_file, 'r') as f:
            expected_result = [line.strip() for line in f]
        self.assertEqual(Counter(result), Counter(expected_result))
        


def test():
    accn_list = ["ERZ4566898@99_SAMEA4765831",
    "ERZ4561431@44_SAMEA4765810",
    "ERZ4561097@157_SAMEA4765866",
    "ERZ4561846@1306_SAMEA4765562",
    "ERZ4561767@437_SAMEA4765564",
    "ERZ4561110@1029_SAMEA4765863",
    "ERZ4560346@597_SAMEA4766029",
    "ERZ4566895@107_SAMEA4765850"]

    json_file = "/spacers_db/crispr-cas-atlas/crispr-cas-atlas-v1.0-VI.jsonl"
    spacer_list, spacer_dict = ah.retrieve_spacers(accn_list=accn_list, json_file=json_file)

    
    spacer_file = "/home/unimelb.edu.au/rbengtsson/work/spacer_analysis/modules/test_files/concensus_spacer_test.list"
    with open(spacer_file, 'r') as f:
        spacer_list = [line.strip() for line in f]
    
    consensus_list = ah.generate_spacer_consensus(spacer_list)
    
    corrected_spacers = ah.correct_spacer_sequence(spacer_list, consensus_list)
    
    long_spacers = ["AAAACGTGTAAAATTATGAACAATTTCAAC",
                    "TTCAATGGAGATATCTCTAAATGGAATGTT",
                    "AGAAAGAAGCACCGTTAGACTACAACGAAC",
                    "CTATGTTCAATGGTCTACTTACGTTGTGTG",
                    "CAGACACCAAAGGAAACACTTTTGAGTTTA",
                    "AAGAGCAGTATGCTCTAGTAAGAATTACAG",
                    "AGTTCTTAAAGTTTATCACTATAAAGTAAA",
                    "AAAAAGATAAATATATACTTGTGGTGTAGT",
                    "TATTCTCGTTGACATTAAAACTAATTCTA"]
    for spacer in long_spacers:
        print(Seq(spacer).reverse_complement())





if __name__ == '__main__':
    # test()
    unittest.main()