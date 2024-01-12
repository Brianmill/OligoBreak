import sys

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def extract_protein_sequence(protein_seq_file_path: str):
    fasta_file = open(protein_seq_file_path, "r")
    protein_name = fasta_file.readline().strip()
    protein_seq = fasta_file.read().strip()
    fasta_file.close()

    record = SeqRecord(
        Seq(protein_seq),
        name=protein_name,
    )
    return record


def store_protein_sequence(protein_sequence_data) -> str:
    name = str(protein_sequence_data.name)
    sequence = str(protein_sequence_data.seq)
    fasta_string = name + sequence
    return fasta_string


def extract_multi_amino_data(multi_amino_file_path: str):
    multi_amino_file = open(multi_amino_file_path, "r")
    multi_amino_data = {}

    for line in multi_amino_file:
        index, amino_acids = line.strip().split(',')
        multi_amino_data[index] = amino_acids

    multi_amino_file.close()
    return multi_amino_data


def store_multi_amino_data(multi_amino_data) -> str:
    output_string = ''

    for index in multi_amino_data:
        output_string = output_string + index + ',' + multi_amino_data[index] + '\n'

    return output_string


def extract_degenerate_codons_key(degenerate_codons_file_path: str):
    degenerate_codon_file = open(degenerate_codons_file_path, "r")
    decodon_data = {}
    decodon_data_binary = {}

    for line in degenerate_codon_file:
        num_decodons, amino_acid_key, decodons, score = line.strip().split('\t')
        decodon_set = tuple(decodons.split())
        """ 
        amino_acid_key is a 20 digit binary number represented as a String. i.e.    10100011010111101101
        the twenty digits each represent the presence or absence of 20 amino acids: ARNDCQEGHILKMFPSTWYV
        they are in alphabetical order based on full amino name
        """
        decodon_data[amino_acid_key] = {
            'num_decodons': num_decodons,
            'decodon_set': decodon_set,
            'score': score
        }
    degenerate_codon_file.close()
    return decodon_data


def store_degenerate_codons_key(decodon_data) -> str:
    output_list = []
    for amino_acid_key in decodon_data:
        num_decodons = decodon_data[amino_acid_key]['num_decodons']
        decodons = ''.join(decodon_data[amino_acid_key]['decodon_set'])
        score = decodon_data[amino_acid_key]['score']
        output_list.append(f"{num_decodons}{amino_acid_key}{decodons}{score}")

    return ''.join(output_list)


def extract_organism_codon_frequencies(codon_frequencies_file_path: str):
    codon_frequency_file = open(codon_frequencies_file_path, "r")
    codon_per_mille_data = {}

    for line in codon_frequency_file:
        amino_acid, per_mille = line.strip().split(',')
        codon_per_mille_data[amino_acid] = per_mille
    codon_frequency_file.close()
    return codon_per_mille_data


def store_organism_codon_frequencies(codon_per_mille_data) -> str:
    output_list = []
    for amino_acid in codon_per_mille_data:
        per_mille = codon_per_mille_data[amino_acid]
        output_list.append(f"{amino_acid},{per_mille}")
    return '\n'.join(output_list)


class aminoAcid:
    @staticmethod
    def binary_to_letter(binary_representation):
        # the list of amino acid letters in alphabetical order according to full name, used for decoding and encoding
        const_amino_list = "ARNDCQEGHILKMFPSTWYV"

        letter_representation = set()
        for index, bit in enumerate(binary_representation):
            if bit == '1':
                letter_representation.add(const_amino_list[index])

        return letter_representation

    @staticmethod
    def letter_to_binary(letter_representation):
        # the list of amino acid letters in alphabetical order according to full name, used for decoding and encoding
        const_amino_list = "ARNDCQEGHILKMFPSTWYV"
        binary_representation = ""
        for index, letter in enumerate(const_amino_list):
            if letter in letter_representation:
                binary_representation = binary_representation + "1"
            else:
                binary_representation = binary_representation + "0"
        return binary_representation


def convert_from_multi_amino_to_multi_decodon(amino_acids, amino_decodon_dict):
    raise NotImplementedError("convert_from_multi_amino_to_mutli_decodon is not yet implemented.")


def convert_from_multi_decodon_to_multi_amino(multi_decodons):
    raise NotImplementedError("convert_from_multi_decodon_to_multi_amino is not yet implemented.")
