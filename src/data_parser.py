from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def extract_protein_sequence(protein_seq_file: str):
    fasta_file = open(protein_seq_file, "r")
    protein_name = fasta_file.readline().strip()
    protein_seq = fasta_file.read().strip()
    fasta_file.close()

    record = SeqRecord(
        Seq(protein_seq),
        name=protein_name,
    )
    return record

def store_protein_sequence(protein_sequence_data) -> str:
    name = protein_sequence_data.name
    print(name)
    sequence = protein_sequence_data.seq
    print(sequence)

    print(str(name) + '\n' + str(sequence))
    fasta_string = str(name) + str(sequence)
    return fasta_string

def extract_multiple_amino_data(multiple_amino_file: str):
    raise (NotImplementedError('extract_multiple_amino_data not yet implemented'))


def store_multiple_amino_data(multiple_amino_data) -> str:
    raise (NotImplementedError('store_multiple_amino_data not yet implemented'))


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


