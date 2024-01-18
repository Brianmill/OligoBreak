from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import namedtuple
import itertools
import numpy as np
from math import floor
import time

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
    'W': ['TGG'],
    'M': ['ATG'],
    'H': ['CAC', 'CAT'],
    'Q': ['CAA', 'CAG'],
    'N': ['AAC', 'AAT'],
    'K': ['AAA', 'AAG'],
    'D': ['GAC', 'GAT'],
    'E': ['GAA', 'GAG'],
    'Y': ['TAC', 'TAT'],
    'F': ['TTC', 'TTT'],
    'C': ['TGC', 'TGT'],
    'I': ['ATA', 'ATC', 'ATT'],
    '*': ['TAA', 'TAG', 'TGA'],
    'P': ['CCA', 'CCC', 'CCG', 'CCT'],
    'T': ['ACA', 'ACC', 'ACG', 'ACT'],
    'V': ['GTA', 'GTC', 'GTG', 'GTT'],
    'A': ['GCA', 'GCC', 'GCG', 'GCT'],
    'G': ['GGA', 'GGC', 'GGG', 'GGT'],
    'S': ['TCA', 'TCC', 'TCG', 'TCT', 'AGC', 'AGT'],
    'L': ['TTA', 'TTG', 'CTA', 'CTC', 'CTG', 'CTT'],
    'R': ['CGA', 'CGC', 'CGG', 'CGT', 'AGA', 'AGG']
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
    """Converts an amino acid set into a corresponding DNA codon set.

    Parameters
    ----------
    amino_set : set
        Set of amino acids.
    aminos_to_codons_mapping : dict
        Mapping of amino acid set to codon set.

    Returns
    -------
    codon_set : set
        Corresponding DNA codon set.
    """
    # if amino set is empty
    if not amino_set:
        return set()
    # if there is only one amino acid in the set, simply lookup the corresponding codon
    elif len(amino_set) == 1:
        return CODON_TABLE[amino_set]

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
    """Returns the corresponding amino acid set for a set of DNA codons

    Parameters
    ----------
    ambiguous_codons : set
        A set of codons made up of ambiguous DNA nucleotides.

    Returns
    -------
    aminos : set
        Corresponding amino acid set.
    """

    amino_set = set()
    for ambiguous_codon in ambiguous_codons:
        # Step 1: converting the three nucleotides from ambiguous to unambiguous
        # {'ANK'} -> {'A', 'ATCG', 'GT'}
        unambiguous_nucleotides = (AMBIGUOUS_UNAMBIGUOUS_NUCLEOTIDE_MAP[nucleotide]
                                   for nucleotide in ambiguous_codon)

        # Step 2: forming all possible codons through a cartesian product of nucleotides
        # {'A', 'ATCG', 'GT'} -> {'AAG', 'ATG', 'ACG', 'AGG', 'AAT', 'ATT', 'ACT', 'AGT'}
        unambiguous_codons = (''.join(codon) for codon in itertools.product(*unambiguous_nucleotides))

        # Step 3: translating codons into amino acids
        # {'AAT', 'ATT', 'ACT', 'AGT'} -> {'K', 'M', 'T', 'R', 'N', 'I', 'S'}
        amino_set.update(set(CODON_TABLE[codon] for codon in unambiguous_codons))

    return amino_set


def convert_multi_amino_to_multi_codon_data(multi_amino_data, aminos_to_codons_mapping):
    multi_codon_data = {}
    for position, aminos in multi_amino_data.items():
        codons = convert_amino_set_to_codon_set(aminos, aminos_to_codons_mapping)
        multi_codon_data[position] = codons
    return multi_codon_data


def convert_multi_codon_to_multi_amino_data(multi_codon_data):
    multi_amino_data = {}
    for position, codons in multi_codon_data.items():
        aminos = convert_codon_set_to_amino_set(codons)
        multi_amino_data[position] = aminos
    return multi_amino_data


def convert_protein_sequence_to_dna_sequence(protein_sequence, codon_frequencies):

    codon_apportionments = codon_apportionment(protein_sequence, codon_frequencies)
    print(str(codon_apportionment))
    pass


def codon_apportionment(protein_sequence, codon_frequencies):
    # apportions amino acids to codons based on codon frequencies using the Hamiltonian method of apportionment

    # codon amounts for all aminos
    global_codon_amounts = {}
    for amino in "VYWTSPFMKLIHGEQCDNRA":
        # determining the number of occurrences of the given amino
        amino_amount = protein_sequence.count(amino)

        codons = REVERSE_CODON_TABLE[amino]

        codon_amounts = [float(codon_frequencies[codon]) for codon in codons]
        total_codon_amount = sum(codon_amounts)

        # divisor = total_population / number_of_seats
        divisor = total_codon_amount / amino_amount

        # quota = population / divisor
        quotas = [(amount / divisor) for amount in codon_amounts]

        # initial = floor(quota)
        initial_apportionment = [floor(quota) for quota in quotas]

        priority_scores = [quotas[i] - initial_apportionment[i] for i in range(len(quotas))]

        remaining_aminos = amino_amount - sum(initial_apportionment)

        while remaining_aminos > 0:
            max_index = priority_scores.index(max(priority_scores))
            initial_apportionment[max_index] += 1
            priority_scores[max_index] = 0
            remaining_aminos -= 1

        for index, codon in enumerate(codons):
            global_codon_amounts[codon] = initial_apportionment[index]

        return global_codon_amounts


# TODO come up with how to test this function


def convert_dna_sequence_to_protein_sequence(dna_sequence):
    pass
