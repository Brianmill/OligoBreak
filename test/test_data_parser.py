from pathlib import Path
import string
from src.data_parser import *
import tempfile
import unittest


TEST_FILES = Path(__file__).parent / 'test_files'


class TestDataParser(unittest.TestCase):

    def test_protein_sequence_extraction(self):
        protein_sequence_file = str(TEST_FILES / 'ferredoxin_oxidoreductase.fasta')

        ## getting the original data from the two files
        with open(protein_sequence_file, 'r') as file:
            original_protein_sequence_data = file.read().replace('\n', '')

        ##running the two files through data extraction and reencoding
        extracted_protein_sequence_data = extract_protein_sequence(protein_sequence_file)
        reencoded_protein_sequence_data = store_protein_sequence(extracted_protein_sequence_data)

        self.assertEqual(
            original_protein_sequence_data,
            reencoded_protein_sequence_data,
            "Original protein sequence data does not match re-encoded data."
        )

    def test_multi_amino_data_extraction(self):
        multiple_amino_file = str(TEST_FILES / 'ferredoxin_oxidoreductase_multi_aminos.csv')

        with open(multiple_amino_file, 'r') as file:
            original_multiamino_data = file.read()

        extracted_multiamino_data = extract_multi_amino_data(multiple_amino_file)
        reencoded_multiamino_data = store_multi_amino_data(extracted_multiamino_data)

        self.assertEqual(
            original_multiamino_data.strip(),
            reencoded_multiamino_data.strip(),
            "Original multi-amino data does not match re-encoded data."
        )

    def test_extract_degenerate_codons_key(self):
        degenerate_codon_file = 'C:\\Users\\yehud\\Desktop\\py_protein_varient_library\\data\\orgainsim_AA_decodon_sets\\Human.txt'

        with open(degenerate_codon_file, 'r') as file:
            original_degenerate_codon_data = file.read().translate(str.maketrans("", "", string.whitespace))

        extracted_degenerate_codon_data = extract_degenerate_codons_key(degenerate_codon_file)
        reencoded_degenerate_codon_data = store_degenerate_codons_key(extracted_degenerate_codon_data)

        self.assertEqual(
            original_degenerate_codon_data.strip(),
            reencoded_degenerate_codon_data.strip(),
            "Original degenerate codon data does not match re-encoded data."
        )

    def test_extract_organism_codon_frequencies(self):
        codon_frequencies_file = 'C:\\Users\\yehud\\Desktop\\py_protein_varient_library\\data\\organisim_codon_frequencies\\human_codon_frequencies.csv'

        with open(codon_frequencies_file, 'r') as file:
            original_codon_frequencies_data = file.read()

        extracted_codon_frequencies_data = extract_organism_codon_frequencies(codon_frequencies_file)
        reencoded_codon_frequencies_data = store_organism_codon_frequencies(extracted_codon_frequencies_data)

        self.assertEqual(
            extracted_codon_frequencies_data.strip(),
            reencoded_codon_frequencies_data.strip(),
            "Original codon frequency data does not match re-encoded data."
        )



if __name__ == '__main__':
    unittest.main()
