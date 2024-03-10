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


class TestDataParserBeforeCommit(unittest.TestCase):
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

    # Todo: write test_convert_protein_to_max_cg_dna_binary_sequence function
    #  Possible tests would be to go in reverse order to check for 'contains'

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


class TestDataParserAfterEdit(unittest.TestCase):

    def test_protein_sequence_extraction(self):
        protein_sequence_file = str(TEST_FILES / 'ferredoxin_oxidoreductase.fasta')

        # getting the original data from the two files
        with open(protein_sequence_file, 'r') as file:
            original_protein_sequence_data = file.read().replace('\n', '')

        # running the two files through data extraction and reencoding
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

    def test_compute_overlap_temperature_dict(self):
        protein_sequence_file = str(TEST_FILES / 'ferredoxin_oxidoreductase.fasta')
        extracted_protein_sequence_data = extract_protein_sequence(protein_sequence_file)
        protein_sequence = extracted_protein_sequence_data.sequence
        max_cg_dna = convert_protein_to_max_cg_dna_binary_sequence(protein_sequence)
        min_cg_dna = convert_protein_to_min_cg_dna_binary_sequence(protein_sequence)

        # check that there are the expected number of entries in the dictionary.
        # check that the lengths do not exceed the maximum length
        # check that the lengths are not smaller than the minimum length

    def test_generate_overlap_dict(self):
        # makes sure that the calculated and expected temperatures match
        min_overlap = 20
        max_overlap = 40
        min_cg = 0
        max_cg = max_overlap
        cg_overlap_temperature_dict = generate_global_overlap_dict(min_overlap, max_overlap, min_cg, max_cg)

        for overlap_length in cg_overlap_temperature_dict:
            cg_temp_tuples = cg_overlap_temperature_dict[overlap_length]
            for cg_temp_tuple in cg_temp_tuples:
                cg_amount = cg_temp_tuple.cg_amount
                temperature = cg_temp_tuple.temperature

                expected_temperature = round(calculate_dna_melting_temp(cg_amount, overlap_length))
                self.assertEqual(
                    expected_temperature,
                    temperature
                )

    def test_generate_overlap_dict_entries(self):
        # checks to see if there are the correct number of entries in the temperature dictionary
        min_overlap = 20
        max_overlap = 40
        min_cg = 0
        max_cg = max_overlap
        cg_overlap_temperature_dict = generate_global_overlap_dict(min_overlap, max_overlap, min_cg, max_cg)

        number_of_tuples_in_dict = 0
        for overlap_length, cg_temp_tuples in cg_overlap_temperature_dict.items():
            number_of_tuples_in_dict = number_of_tuples_in_dict + len(cg_temp_tuples)

        overlap_range = max_overlap - min_overlap + 1
        cg_range = max_cg - min_cg + 1
        number_of_expected_tuples = overlap_range * cg_range
        self.assertEqual(
            number_of_expected_tuples,
            number_of_tuples_in_dict
        )

    def test_generate_overlap_dict_from_sequence(self):
        dna_binary_sequence = '0100101101010101101010101010100011000010'
        min_overlap = 20
        max_overlap = 40
        overlap_dict = generate_overlap_dict_from_sequence(dna_binary_sequence, min_overlap, max_overlap)
        expected_number_of_entries = max_overlap - min_overlap + 1
        actual_number_of_entries = len(overlap_dict)
        self.assertEqual(
            expected_number_of_entries,
            actual_number_of_entries
        )

    def test_generate_overlap_dict_from_sequence_with_insufficient_sequence_size(self):
        dna_binary_sequence = '0100101101010101101'
        min_overlap = 20
        max_overlap = 40
        self.assertRaises(ValueError, generate_overlap_dict_from_sequence, dna_binary_sequence, min_overlap,
                          max_overlap)

    def test_generate_overlap_dict_from_sequence_with_max_overlap_greater_than_min_overlap(self):
        dna_binary_sequence = '0100101101010101101010101010100011000010'
        min_overlap = 40
        max_overlap = 20
        self.assertRaises(ValueError, generate_overlap_dict_from_sequence, dna_binary_sequence, min_overlap,
                          max_overlap)

    def test_generate_overlap_dict_from_sequence_ensuring_that_sequence_is_reversed(self):
        dna_binary_sequence = '0100101101010101101010101010100011000010'
        min_overlap = 20
        max_overlap = 40

        reversed_sequence = dna_binary_sequence[::-1]
        expected_cg_amount_array = []
        num_cg = 0
        for dna_index, binary_nucleotide in enumerate(reversed_sequence):
            num_cg = num_cg + int(binary_nucleotide)
            overlap_length = dna_index + 1
            if min_overlap <= overlap_length <= max_overlap:
                expected_cg_amount_array.append(num_cg)

        actual_cg_amount_array = []
        overlap_dict = generate_overlap_dict_from_sequence(dna_binary_sequence, min_overlap, max_overlap)
        for overlap, cg_temp_tuple in overlap_dict.items():
            cg_amount = cg_temp_tuple.cg_amount
            actual_cg_amount_array.append(cg_amount)
        self.assertCountEqual(
            expected_cg_amount_array,
            actual_cg_amount_array
        )

    # ToDo: Write a test to ensure that the sequence is infact reversed

    def test_generate_overlap_dict_array(self):
        min_overlap = 20
        max_overlap = 40
        cg_dna_binary_sequence = "0010101000100000101000001001001000000000100000000000101000001100000100110101000000101001001001000000100101"

        sequence_length = len(cg_dna_binary_sequence)
        overlap_amount = sequence_length - min_overlap
        overlap_range = max_overlap - min_overlap
        adjustment_factor = (min_overlap * (min_overlap - 1)) / 2
        expected_number_of_total_entries = overlap_amount * overlap_range - adjustment_factor

        overlap_dict_array = generate_overlap_dict_array(cg_dna_binary_sequence, min_overlap, max_overlap)
        actual_number_of_total_entries = 0
        filtered_overlap_dict_array = [value for value in overlap_dict_array if value is not None]
        for overlap_dict in filtered_overlap_dict_array:
            actual_number_of_total_entries = actual_number_of_total_entries + len(overlap_dict)

        self.assertEqual(expected_number_of_total_entries,
                         actual_number_of_total_entries,
                         str(f"expected: {expected_number_of_total_entries}, actual: {actual_number_of_total_entries}"))

    def test_generate_filtered_melting_temp_dict_array(self):
        min_overlap = 20
        max_overlap = 40
        min_cg = 0
        max_cg = max_overlap

        max_cg_dna = "0010111010110011111011011011011010011010110010011010111010011110010110110111010010111011011011010011"
        min_cg_dna = "0010101000100000101000001001001000000000100000000000101000001100000100110101000000101001001001000000"
        max_temp_dict_list = generate_overlap_dict_array(max_cg_dna, min_overlap, max_overlap)
        min_temp_dict_list = generate_overlap_dict_array(min_cg_dna, min_overlap, max_overlap)

        global_temp_dict = generate_global_overlap_dict(min_overlap, max_overlap, min_cg, max_cg)

        filtered_melting_temp_dict_array = generate_filtered_melting_temp_dict_array(global_temp_dict,
                                                                                     min_temp_dict_list,
                                                                                     max_temp_dict_list)
        self.assertTrue(False)


if __name__ == '__main__':
    cProfile.run("unittest.main()", sort='tottime')
