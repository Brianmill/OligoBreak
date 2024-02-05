from collections import namedtuple
import cProfile
import itertools
import pandas as pd
from pathlib import Path
import string
from src.data_parser import *
import sys
import unittest

sys.path.append("..")

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

        extracted_degenerate_codon_data = extract_amino_codon_dict(degenerate_codon_file)
        # storing the dictionary for use in future tests
        self.amino_acid_to_codon_dict = extracted_degenerate_codon_data
        reencoded_degenerate_codon_data = store_amino_codon_dict(extracted_degenerate_codon_data)

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
        amino_acid_set = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                          'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'}
        amino_acid_combinations = set()
        for n in range(0, 21):
            amino_acid_combinations.update(set(itertools.combinations(amino_acid_set, n)))

        amino_acid_to_codon_file = 'C:\\Users\\yehud\\Desktop\\py_protein_varient_library\\data\\orgainsim_AA_decodon_sets\\Human.txt'
        amino_acid_to_codon_dict = extract_amino_codon_dict(amino_acid_to_codon_file)

        for index, expected_amino_set in enumerate(amino_acid_combinations):
            codon_set = convert_amino_set_to_codon_set(expected_amino_set, amino_acid_to_codon_dict)
            actual_amino_set = convert_codon_set_to_amino_set(codon_set)
            self.assertCountEqual(
                expected_amino_set,
                actual_amino_set,
            )

    # TODO write a test that checks if DNA produced through reverse translation has right codon frequencies.
    def test_protein_to_dna_roundtrip(self):
        protein_sequence_file = str(TEST_FILES / 'ferredoxin_oxidoreductase.fasta')
        protein_data = extract_protein_sequence(protein_sequence_file)
        original_protein_sequence = protein_data.sequence

        codon_frequencies_file = 'C:\\Users\\yehud\\Desktop\\py_protein_varient_library\\data\\organisim_codon_frequencies\\human_codon_frequencies.csv'
        codon_frequencies = extract_organism_codon_frequencies(codon_frequencies_file)

        dna_sequence = convert_protein_sequence_to_dna_sequence(original_protein_sequence, codon_frequencies)
        roundtrip_protein_sequence = convert_dna_sequence_to_protein_sequence(dna_sequence)

        self.assertEqual(
            original_protein_sequence,
            roundtrip_protein_sequence,
        )

    def test_multi_amino_to_multi_codon_roundtrip(self):
        multiple_amino_file = str(TEST_FILES / 'ferredoxin_oxidoreductase_multi_aminos.csv')
        original_multi_amino_data = extract_multi_amino_data(multiple_amino_file)

        amino_codon_dict_file = 'C:\\Users\\yehud\\Desktop\\py_protein_varient_library\\data\\orgainsim_AA_decodon_sets\\Human.txt'
        amino_codon_dict = extract_amino_codon_dict(amino_codon_dict_file)

        multi_codon_data = convert_multi_amino_to_multi_codon_data(original_multi_amino_data, amino_codon_dict)
        roundtrip_multi_codon_data = convert_multi_codon_to_multi_amino_data(multi_codon_data)

        self.assertCountEqual(
            original_multi_amino_data,
            roundtrip_multi_codon_data,
        )

    def test_dna_temperature_calculator(self):

        file_path = str(TEST_FILES / 'primer.dat')
        df = pd.read_csv(file_path, delimiter='\s+', header=0)
        for ind in df.index:
            dna_length = df['Len'][ind]
            cg_percentage = df['%CG'][ind]
            cg_amount = dna_length * cg_percentage
            calculated_melting_temp = calculate_dna_melting_temp(cg_amount, dna_length)
            expected_melting_temp = float(df['Basic'][ind])

            self.assertAlmostEqual(
                calculated_melting_temp,
                expected_melting_temp,
                places=1
            )

    # Todo: test_overlap_calculator is incomplete. Possible tests would be to go in reverse order to check for 'contains'
    def test_overlap_calculator(self):
        protein_sequence_file = str(TEST_FILES / 'ferredoxin_oxidoreductase.fasta')
        extracted_protein_sequence_data = extract_protein_sequence(protein_sequence_file)
        protein_sequence = extracted_protein_sequence_data.sequence
        max_cg_dna = convert_protein_to_max_cg_dna_binary_sequence(protein_sequence)
        min_cg_dna = convert_protein_to_min_cg_dna_binary_sequence(protein_sequence)

    def test_compute_cg_content_over_a_dna_sequence(self):
        protein_sequence_file = str(TEST_FILES / 'ferredoxin_oxidoreductase.fasta')
        extracted_protein_sequence_data = extract_protein_sequence(protein_sequence_file)
        protein_sequence = extracted_protein_sequence_data.sequence
        max_cg_dna_binary_sequence = convert_protein_to_max_cg_dna_binary_sequence(protein_sequence)
        max_cg_amount = compute_cg_amount_over_a_cg_dna_binary_sequence(max_cg_dna_binary_sequence)

        expected_sequence_length = len(max_cg_dna_binary_sequence)
        actual_sequence_length = len(max_cg_amount)
        self.assertEqual(
            expected_sequence_length,
            actual_sequence_length
        )

        expected_number_of_cg_in_dna_sequence = max_cg_amount.count('1')
        actual_number_of_cg_in_dna_sequence = max_cg_amount[0]
        self.assertEqual(
            expected_number_of_cg_in_dna_sequence,
            actual_number_of_cg_in_dna_sequence
        )

    def test_compute_dna_melting_temp_dict(self):
        min_overlap = 20
        max_overlap = 40
        min_cg = 0
        max_cg = max_overlap
        melting_temp_dict = compute_dna_melting_temp_dict(min_overlap, max_overlap, min_cg, max_cg)

        number_of_tuples_in_dict = 0

        for temperature in melting_temp_dict:
            cg_length_tuples = melting_temp_dict[temperature]

            for cg_length_tuple in cg_length_tuples:
                cg_amount = cg_length_tuple[0]
                overlap_length = cg_length_tuple[1]

                expected_temperature = int(calculate_dna_melting_temp(cg_amount, overlap_length))
                self.assertEqual(
                    temperature,
                    expected_temperature
                )
                number_of_tuples_in_dict = number_of_tuples_in_dict + len(cg_length_tuples)

        overlap_range = max_overlap - min_overlap
        cg_range = max_cg - min_cg
        number_of_expected_tuples = overlap_range * cg_range
        self.assertEqual(
            number_of_tuples_in_dict,
            number_of_expected_tuples
        )


if __name__ == '__main__':
    cProfile.run("unittest.main()", sort='tottime')
