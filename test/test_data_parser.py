import sys

sys.path.append("..")
import cProfile
import itertools
from pathlib import Path
import string
from src.data_parser import *
import unittest

TEST_FILES = Path(__file__).parent / 'test_files'


class TestDataParser(unittest.TestCase):
    amino_acid_to_codon_dict = None

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
        # storing the dictionary for use in future tests
        self.amino_acid_to_codon_dict = extracted_degenerate_codon_data
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
            original_codon_frequencies_data.strip(),
            reencoded_codon_frequencies_data.strip(),
            "Original codon frequency data does not match re-encoded data."
        )

    def test_amino_acid_to_codon_roundtrip(self):
        amino_acid_set = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W',
                          'Y'}
        amino_acid_combinations = set()
        for n in range(0, 21):
            amino_acid_combinations.update(set(itertools.combinations(amino_acid_set, n)))

        # test to ensure that all combinations of amino acids are present
        expected_length = 2 ** 20
        self.assertEqual(expected_length, len(amino_acid_combinations))

        amino_acid_to_codon_dict = self.amino_acid_to_codon_dict
        if amino_acid_to_codon_dict is None:
            amino_acid_to_codon_file = 'C:\\Users\\yehud\\Desktop\\py_protein_varient_library\\data\\orgainsim_AA_decodon_sets\\Human.txt'
            amino_acid_to_codon_dict = extract_degenerate_codons_key(amino_acid_to_codon_file)

        for index, expected_amino_set in enumerate(amino_acid_combinations):
            codon_set = convert_from_amino_set_to_codon_set(expected_amino_set, amino_acid_to_codon_dict)
            actual_amino_set = convert_from_codon_set_to_amino_set(codon_set)
            self.assertCountEqual(
                expected_amino_set,
                actual_amino_set,
            )
            if index % 10485 == 0:
                percent = index / 10485
                print(percent)


if __name__ == '__main__':
    cProfile.run("unittest.main()", sort='ncalls')
