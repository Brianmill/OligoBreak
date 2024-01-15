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
            'num_codons': num_decodons,
            'codon_set': decodon_set,
            'score': score
        }
    degenerate_codon_file.close()
    return decodon_data


def store_degenerate_codons_key(decodon_data) -> str:
    output_list = []
    for amino_acid_key in decodon_data:
        num_decodons = decodon_data[amino_acid_key]['num_codons']
        decodons = ''.join(decodon_data[amino_acid_key]['codon_set'])
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

    # converts a binary representation of an amino acid set to its letter representation.
    # i.e. 10000000100000001011 -> {'A', 'H', 'T', 'Y', 'V'}
    # The logic is as follows:
    # 10000000100000001011
    # |.......|.......|.||
    # ARNDCQEGHILKMFPSTWYV

    # converts a letter representation of an amino acid set to its binary representation.
    # i.e. {'A', 'H', 'T', 'Y', 'V'} -> 10000000100000001011
    # The logic is as follows:
    # ARNDCQEGHILKMFPSTWYV
    # |.......|.......|.||
    # 10000000100000001011


def convert_from_amino_set_to_codon_set(amino_set, aminos_to_codons):
    if not amino_set:
        return set()

    # the list of amino acid letters in alphabetical order according to full name, used for decoding and encoding
    amino_list = "VYWTSPFMKLIHGEQCDNRA"

    binary_representation = ""
    for index, amino in enumerate(amino_list):
        if amino in amino_set:
            binary_representation = binary_representation + "1"
        else:
            binary_representation = binary_representation + "0"

    codon_set = aminos_to_codons[binary_representation]['codon_set']
    return codon_set


def convert_from_codon_set_to_amino_set(multi_codon_set):
    # converts from a set of codons to the set of amino acids that are coded for by the codon set

    full_codon_set = set()  # all the codons derived from every multi-codon
    for multi_codon in multi_codon_set:
        partial_codon_set = convert_from_multi_codon_to_single_codon_set(
            multi_codon)  # codons derived from individual multi-codon
        full_codon_set.update(partial_codon_set)

    amino_set = convert_from_codon_set_to_amino_acid_set(full_codon_set)
    return amino_set


def convert_from_multi_codon_to_single_codon_set(multi_codon):
    # converts from multi-codon to corresponding codon sets
    # i.e. {'ANT'} -> {'AAT', 'ATT', 'ACT', 'AGT'}

    iupac_dict = {  # a dictionary that maps IUPAC DNA nucleotide codes to corresponding DNA nucleotide sets
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

    # converting the three nucleotides from IUPAC into the corresponding nucleotide sets
    first_nucleotide_set = iupac_dict[multi_codon[0]]
    second_nucleotide_set = iupac_dict[multi_codon[1]]
    third_nucleotide_set = iupac_dict[multi_codon[2]]

    # stitching the nucleotides together to form all possible codons
    codon_set = set()
    for first_nucleotide in first_nucleotide_set:
        for second_nucleotide in second_nucleotide_set:
            for third_nucleotide in third_nucleotide_set:
                codon = first_nucleotide + second_nucleotide + third_nucleotide
                codon_set.add(codon)
    return codon_set


def convert_from_codon_set_to_amino_acid_set(codon_set):
    # {'AAT', 'ATT', 'ACT', 'AGT'} -> {'N', 'I', 'T', 'S'}

    codon_table = {'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S', 'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
                   'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*', 'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
                   'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
                   'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
                   'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M', 'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
                   'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K', 'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
                   'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
                   'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E', 'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
                   }

    amino_acid_set = set()
    for codon in codon_set:
        amino_acid = codon_table[codon]
        amino_acid_set.add(amino_acid)

    return amino_acid_set
