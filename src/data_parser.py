from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import namedtuple
import itertools
import numpy as np
from math import floor
import sys
import random

CODON_TABLE = {
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S', 'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*', 'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M', 'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K', 'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E', 'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
}

REVERSE_CODON_TABLE = {
    'A': ['GCA', 'GCC', 'GCG', 'GCT'],
    'C': ['TGC', 'TGT'],
    'D': ['GAC', 'GAT'],
    'E': ['GAA', 'GAG'],
    'F': ['TTC', 'TTT'],
    'G': ['GGA', 'GGC', 'GGG', 'GGT'],
    'H': ['CAC', 'CAT'],
    'I': ['ATA', 'ATC', 'ATT'],
    'K': ['AAA', 'AAG'],
    'L': ['TTA', 'TTG', 'CTA', 'CTC', 'CTG', 'CTT'],
    'M': ['ATG'],
    'N': ['AAC', 'AAT'],
    'P': ['CCA', 'CCC', 'CCG', 'CCT'],
    'Q': ['CAA', 'CAG'],
    'R': ['CGA', 'CGC', 'CGG', 'CGT', 'AGA', 'AGG'],
    'S': ['TCA', 'TCC', 'TCG', 'TCT', 'AGC', 'AGT'],
    'T': ['ACA', 'ACC', 'ACG', 'ACT'],
    'V': ['GTA', 'GTC', 'GTG', 'GTT'],
    'W': ['TGG'],
    'Y': ['TAC', 'TAT'],
    '*': ['TAA', 'TAG', 'TGA'],
}

CODON_WITH_LOWEST_CG_CONTENT = {
    'A': '110', 'C': '010', 'D': '100', 'E': '100',
    'F': '000', 'G': '110', 'H': '100', 'I': '000',
    'K': '000', 'L': '000', 'M': '001', 'N': '000',
    'P': '110', 'Q': '100', 'R': '010', 'S': '010',
    'T': '010', 'V': '100', 'W': '011', 'Y': '000',
    '*': '000',
}

CODON_WITH_HIGHEST_CG_CONTENT = {
    'A': '111', 'C': '011', 'D': '101', 'E': '101',
    'F': '001', 'G': '111', 'H': '101', 'I': '001',
    'K': '001', 'L': '101', 'M': '001', 'N': '001',
    'P': '111', 'Q': '101', 'R': '111', 'S': '011',
    'T': '011', 'V': '101', 'W': '011', 'Y': '001',
    # * can be coded as 010, but this does not matter much, since it is a stop codon after all
    '*': '001',
}

# a dictionary that maps ambiguous DNA nucleotides to unambiguous DNA nucleotides
AMBIGUOUS_UNAMBIGUOUS_NUCLEOTIDE_MAP = {
    'A': {'A'},
    'C': {'C'},
    'G': {'G'},
    'T': {'T'},
    'R': {'A', 'G'},
    'Y': {'C', 'T'},
    'S': {'G', 'C'},
    'W': {'A', 'T'},
    'K': {'G', 'T'},
    'M': {'A', 'C'},
    'B': {'C', 'G', 'T'},
    'D': {'A', 'G', 'T'},
    'H': {'A', 'C', 'T'},
    'V': {'A', 'C', 'G'},
    'N': {'A', 'T', 'C', 'G'}
}


def extract_protein_sequence(protein_seq_file_path: str):
    """Extracts protein and name information from a protein FASTA file.

    Parameters
    ----------
    protein_seq_file_path: str
        The path of a protein sequence FASTA file

    Returns
    -------
    Protein
        namedtuple('name', 'sequence')
    """
    fasta_file = open(protein_seq_file_path, "r")
    protein_name = fasta_file.readline().strip()
    protein_seq = fasta_file.read().strip()
    fasta_file.close()

    Protein = namedtuple('Protein', ['name', 'sequence'])
    protein = Protein(protein_name, protein_seq)
    return protein


def extract_multi_amino_data(multi_amino_file_path: str):
    """Extracts multi_amino data and stores them in a set.

    Multi amino data is locations in the protein sequence where more than one amino
    acid is desired. This is indicated through the use of index, amino list pairs.

    Parameters
    ----------
    multi_amino_file_path : str
        the path of the file containing the multi amino data

    Returns
    -------
    multi_amino_data : dict
        extracted multi amino data. multi_amino_data[position] = (amino_acids)
    """
    multi_amino_file = open(multi_amino_file_path, "r")
    multi_amino_data = {}

    for line in multi_amino_file:
        index, amino_acids = line.strip().split(',')
        multi_amino_data[index] = tuple(amino_acids)

    multi_amino_file.close()
    return multi_amino_data


def extract_organism_codon_frequencies(codon_frequencies_file_path: str):
    """Extracts codon frequencies for an organism and stores them in a dictionary

    Parameters
    ----------
    codon_frequencies_file_path : str
        Path to the file containing the organism codon frequencies.

    Returns
    -------
    dict
        codon frequencies stored as permille (like percent but 1000 instead of 100)
        codon_per_mille_data[codon] = per_mille
    """
    codon_frequency_file = open(codon_frequencies_file_path, "r")
    codon_per_mille_data = {}

    for line in codon_frequency_file:
        codon, per_mille = line.strip().split(',')
        codon_per_mille_data[codon] = per_mille
    codon_frequency_file.close()
    return codon_per_mille_data


def extract_amino_codon_dict(amino_codon_file_path: str):
    """ Extracts mappings between amino acid sets and corresponding codon sets.

        Amino_set Codon_set correspondences are derived using the codon genie
        tool. The tool finds the optimal codon set, using ambiguous DNA, that
        codes for a given amino acid set.


        Parameters
        ----------
        amino_codon_file_path : str
            Path to the file containing amino_set codon_set correspondences.

        Returns
        -------
        dict
            Extracted amino_set codon_set equivalencies.
            amino_codon_dict[amino_acid_key] = {
                'num_codons',
                'codon_set',
                'score'
            }
    """
    amino_codon_file = open(amino_codon_file_path, "r")
    amino_codon_dict = {}

    for line in amino_codon_file:
        num_codons, amino_acid_key, codons, score = line.strip().split('\t')
        # stored as tuple, since order needs to be conserved to ensure roundtrip test is successful
        codons = tuple(codons.split())

        # Amino acid key has the format of a 20 digit binary number.
        # i.e. {'V', 'Y', 'T', 'H', 'A'} -> 11010000000100000001
        # The logic is as follows:
        # VYWTSPFMKLIHGEQCDNRA   List of amino acids in reverse alphabetical order, according to full name
        # ||.|.......|.......|
        # 11010000000100000001

        # TODO convert data structure from a set into a namedtuple()
        amino_codon_dict[amino_acid_key] = {
            'num_codons': num_codons,
            'codon_set': codons,
            'score': score
        }
    amino_codon_file.close()
    return amino_codon_dict


def store_protein_sequence(protein_sequence_data) -> str:
    # Converts protein data back into a String for roundtrip testing purposes

    name = protein_sequence_data.name
    sequence = protein_sequence_data.sequence
    fasta_string = name + sequence
    return fasta_string


def store_multi_amino_data(multi_amino_data) -> str:
    # Converts multi amino data back into a String for roundtrip testing purposes

    output_string = ''

    for index in multi_amino_data:
        multi_amino_string = ''.join(str(amino) for amino in multi_amino_data[index])
        output_string = output_string + index + ',' + multi_amino_string + '\n'

    return output_string


def store_amino_codon_dict(amino_codon_dict) -> str:
    # Converts amino_set codon_set correspondences back into a String for roundtrip testing purposes

    output_list = []
    for amino_acid_key in amino_codon_dict:
        num_decodons = amino_codon_dict[amino_acid_key]['num_codons']
        decodons = ''.join(amino_codon_dict[amino_acid_key]['codon_set'])
        score = amino_codon_dict[amino_acid_key]['score']
        output_list.append(f"{num_decodons}{amino_acid_key}{decodons}{score}")

    return ''.join(output_list)


def store_organism_codon_frequencies(codon_per_mille_data) -> str:
    # Converts codon frequencies back into a String for roundtrip testing purposes

    output_list = []
    for amino_acid in codon_per_mille_data:
        per_mille = codon_per_mille_data[amino_acid]
        output_list.append(f"{amino_acid},{per_mille}")
    return '\n'.join(output_list)


def convert_amino_set_to_codon_set(amino_set, aminos_to_codons_mapping):
    if not amino_set:
        return set()

    # Step 1: Convert amino acid into binary representation

    # Note: amino acids are encoded in reverse alphabetical order according to full name
    amino_list = "VYWTSPFMKLIHGEQCDNRA"
    binary_representation = ""
    for index, amino in enumerate(amino_list):
        if amino in amino_set:
            binary_representation = binary_representation + "1"
        else:
            binary_representation = binary_representation + "0"

    # Step 2: pass the binary representation into the dictionary to get the codon set
    codon_set = aminos_to_codons_mapping[binary_representation]['codon_set']
    return codon_set


def convert_codon_set_to_amino_set(ambiguous_codons):
    amino_set = set()
    for ambiguous_codon in ambiguous_codons:
        # Step 1: converting the three nucleotides from ambiguous to unambiguous
        # {'ANK'} -> {'A', 'ATCG', 'GT'}
        unambiguous_nucleotides = (AMBIGUOUS_UNAMBIGUOUS_NUCLEOTIDE_MAP[nucleotide] for nucleotide in ambiguous_codon)

        # Step 2: forming all possible codons through a cartesian product of nucleotides
        # {'A', 'ATCG', 'GT'} -> {'AAG', 'ATG', 'ACG', 'AGG', 'AAT', 'ATT', 'ACT', 'AGT'}
        unambiguous_codons = (''.join(codon) for codon in itertools.product(*unambiguous_nucleotides))

        # Step 3: translating codons into amino acids
        # {'AAT', 'ATT', 'ACT', 'AGT'} -> {'K', 'M', 'T', 'R', 'N', 'I', 'S'}
        amino_set.update(set(CODON_TABLE[codon] for codon in unambiguous_codons))

    return amino_set


def convert_multi_amino_to_multi_codon_data(multi_amino_data, aminos_to_codons_mapping):
    # converts a list of multi amino positions to a list of multi codon positions

    multi_codon_data = {}
    for position, aminos in multi_amino_data.items():
        codons = convert_amino_set_to_codon_set(aminos, aminos_to_codons_mapping)
        multi_codon_data[position] = codons
    return multi_codon_data


def convert_multi_codon_to_multi_amino_data(multi_codon_data):
    # converts a list of multi codon positions to a list of multi amino positions

    multi_amino_data = {}
    for position, codons in multi_codon_data.items():
        aminos = convert_codon_set_to_amino_set(codons)
        multi_amino_data[position] = aminos
    return multi_amino_data


# TODO refactor function to reduce nesting
def convert_protein_sequence_to_dna_sequence(protein_sequence, codon_frequencies):
    amino_apportionment = apportion_aminos_to_codons(protein_sequence, codon_frequencies)

    codon_list = []
    for amino_acid in protein_sequence:
        current_amino_apportionment = amino_apportionment[amino_acid]
        number_of_aminos_left = sum(current_amino_apportionment.values())
        chosen_codon_index = random.randint(0, number_of_aminos_left - 1)

        codon_to_be_inserted = ''
        cumulative_count = 0
        for codon, count in current_amino_apportionment.items():
            cumulative_count += count
            if chosen_codon_index <= cumulative_count:
                codon_to_be_inserted = codon
                current_amino_apportionment[codon] -= 1
                if current_amino_apportionment[codon] == 0:
                    del current_amino_apportionment[codon]
                break
        codon_list.append(codon_to_be_inserted)

    dna_sequence = ''.join(codon_list)
    return dna_sequence


def apportion_aminos_to_codons(protein_sequence, codon_frequencies):
    # apportions amino acids to codons based on codon frequencies using the Hamiltonian method of apportionment
    # The following video was heavily relied on to write the function: https://www.youtube.com/watch?v=YWfEqWLz9pc

    codon_allocations = {}
    for amino_acid in "VYWTSPFMKLIHGEQCDNRA":
        amino_acid_amount = protein_sequence.count(amino_acid)
        codons = REVERSE_CODON_TABLE[amino_acid]

        codon_amounts = [float(codon_frequencies[codon]) for codon in codons]
        total_codon_amount = sum(codon_amounts)

        divisor = total_codon_amount / amino_acid_amount
        codon_quotas = [(amount / divisor) for amount in codon_amounts]
        codon_apportionments = [floor(codon_quota) for codon_quota in codon_quotas]
        codon_priority_scores = [codon_quotas[i] - codon_apportionments[i] for i in range(len(codon_quotas))]

        remaining_aminos = amino_acid_amount - sum(codon_apportionments)
        while remaining_aminos > 0:
            max_index = codon_priority_scores.index(max(codon_priority_scores))
            codon_apportionments[max_index] += 1
            codon_priority_scores[max_index] = 0
            remaining_aminos -= 1

        amino_codon_amounts = {codons[i]: codon_apportionments[i] for i in range(len(codons))}
        codon_allocations[amino_acid] = amino_codon_amounts

    return codon_allocations


def convert_dna_sequence_to_protein_sequence(dna_sequence):
    print(dna_sequence)
    codon_list = [dna_sequence[i:i + 3] for i in range(0, len(dna_sequence), 3)]
    amino_list = [CODON_TABLE[codon] for codon in codon_list]
    protein_sequence = ''.join(amino_list)
    return protein_sequence


def convert_protein_to_high_melting_point_dna(protein_sequence):
    codon_list = [CODON_WITH_HIGHEST_CG_CONTENT[amino] for amino in protein_sequence]
    print(codon_list)
    dna_sequence = ''.join(codon_list)
    return dna_sequence


def convert_protein_to_low_melting_point_dna(protein_sequence):
    codon_list = [CODON_WITH_LOWEST_CG_CONTENT[amino] for amino in protein_sequence]
    dna_sequence = ''.join(codon_list)
    return dna_sequence


# Todo error handling: Division by 0
def calculate_dna_melting_temp(cg_amount, dna_length):
    # basic melting temperature equation obtained from http://biotools.nubic.northwestern.edu/OligoCalc.html#helpthermo

    if dna_length == 0:
        print("Error: Division by zero is not allowed. Please provide a non-zero denominator.")
        sys.exit(1)

    melting_temperature = 64.9 + 41 * (cg_amount - 16.4) / dna_length
    return melting_temperature


def compute_dna_melting_temp_dict(min_overlap, max_overlap, min_cg, max_cg):
    """
    dictionary keys are melting temperatures and dictionary values are lists containing all the
    (cg_amount, overlap_length) tuples that result in the given melting temperature
    """
    melting_temp_dict = {}
    for overlap_length in range(min_overlap, max_overlap + 1):
        for cg_amount in range(min_overlap, max_overlap + 1):
            melting_temp = int(calculate_dna_melting_temp(cg_amount, overlap_length))
            amount_overlap_tuple = (cg_amount, overlap_length)
            melting_temp_dict.setdefault(melting_temp, []).append(amount_overlap_tuple)
    return melting_temp_dict


def calculate_overlap(cg_dna_sequence, desired_melting_temp, dna_index):
    # values derived from minimum and maximum recommended overlap length for Gibbs Isothermal Assembly
    min_overlap = 20
    max_overlap = 40
    counter: int = 0
    for overlap_length in range(min_overlap, max_overlap):

        if dna_index - overlap_length < 0:
            continue
        cg_amount = cg_dna_sequence[dna_index - overlap_length: dna_index + 1].count('1')
        calculated_melting_temp = calculate_dna_melting_temp(cg_amount, overlap_length)
        if abs(calculated_melting_temp - desired_melting_temp) < 2.5:
            counter = counter + 1

            # return calculated_melting_temp

    return counter


def calculate_overlap_matrix(min_temp, max_temp, low_cg_dna, high_cg_dna):
    length = len(low_cg_dna)
    temp_range = max_temp - min_temp

    minimum_overlap_matrix = np.zeros((temp_range, length)).astype(int)
    maximum_overlap_matrix = np.zeros((temp_range, length)).astype(int)

    for temp_index in range(temp_range):
        desired_temp = min_temp + temp_index
        for dna_index in range(length):
            minimum_overlap = calculate_overlap(low_cg_dna, desired_temp, dna_index)
            maximum_overlap = calculate_overlap(high_cg_dna, desired_temp, dna_index)

            minimum_overlap_matrix[temp_index][dna_index] = minimum_overlap
            maximum_overlap_matrix[temp_index][dna_index] = maximum_overlap

    with np.printoptions(threshold=np.inf, linewidth=10000):
        print(minimum_overlap_matrix)
        print(maximum_overlap_matrix)

    return minimum_overlap_matrix, maximum_overlap_matrix
