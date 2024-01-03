from pathlib import Path
import tempfile
import unittest
from src.data_parser import *

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

    def test_multiamino_data_extraction(self):
        multiple_amino_file = str(TEST_FILES / 'ferredoxin_oxidoreductase_degen_condons.txt')

        with open(multiple_amino_file, 'r') as file:
            original_multiamino_data = file.read()

        extracted_multiamino_data = extract_protein_sequence(multiple_amino_file)
        reencoded_multiamino_data = store_protein_sequence(extracted_multiamino_data)

        self.assertEqual(
            original_multiamino_data,
            reencoded_multiamino_data,
            "Original multi-amino data does not match re-encoded data."
        )




if __name__ == '__main__':
    unittest.main()

