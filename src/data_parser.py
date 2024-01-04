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


def extract_multiple_amino_data(multi_amino_file_path: str):
    multi_amino_file = open(multi_amino_file_path, "r")
    multi_amino_data = {}

    for line in multi_amino_file:
        index, amino_acids = line.strip().split(',')
        multi_amino_data[index] = amino_acids

    multi_amino_file.close()
    return multi_amino_data


def store_multiple_amino_data(multi_amino_data) -> str:

    output_string = ''

    for index in multi_amino_data:
        output_string = output_string + index + ',' + multi_amino_data[index] + '\n'

    return output_string


def extract_organism_codon_frequencies(codon_frequencies_file: str):

    organism_file = open(codon_frequencies_file, "r")

    # Adjust this part according to your file format
    # This example assumes each line in the file is a key-value pair separated by a comma
    data = {}
    for line in file:
        key, value = line.strip().split(',')
        data[key] = value

    return data
    raise (NotImplementedError('extract_organism_codon_frequencies not yet implemented'))


def store_organism_codon_frequencies() -> str:
    raise (NotImplementedError('store_organism_codon_frequencies not yet implemented'))


def extract_degenerate_codons_key(degenerate_codons_key_file: str):
    raise (NotImplementedError('extract_organism_codon_frequencies not yet implemented'))


def store_degenerate_codons_key() -> str:
    raise (NotImplementedError('store_organism_codon_frequencies not yet implemented'))


