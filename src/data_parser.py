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
        output_string = output_string + index + ' ,' + multi_amino_data[index] + '\n'

    return output_string


def extract_degenerate_codons_key(degenerate_codons_file_path: str):

    degenerate_codon_file = open(degenerate_codons_file_path, "r")
    decodon_data = {}

    for line in degenerate_codon_file:
        num_decodons, amino_acid_key, decodons, score = line.strip().split('\t')
        decodon_set = tuple(decodons.split())
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
        codon_per_mille_data[amino_acid] = {
            per_mille
        }
    codon_frequency_file.close()
    return codon_per_mille_data


def store_organism_codon_frequencies() -> str:
    raise (NotImplementedError('store_organism_codon_frequencies not yet implemented'))
