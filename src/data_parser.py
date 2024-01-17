import itertools
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

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
    seqRecord
        the extracted protein information
    """
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
    # Converts protein data back into a String for roundtrip testing purposes

    name = str(protein_sequence_data.name)
    sequence = str(protein_sequence_data.seq)
    fasta_string = name + sequence
    return fasta_string


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
        extracted multi amino data. multi_amino_data[position] = 'amino_acids'
    """
    multi_amino_file = open(multi_amino_file_path, "r")
    multi_amino_data = {}

    for line in multi_amino_file:
        index, amino_acids = line.strip().split(',')
        multi_amino_data[index] = amino_acids

    multi_amino_file.close()
    return multi_amino_data


def store_multi_amino_data(multi_amino_data) -> str:
    # Converts multi amino data back into a String for roundtrip testing purposes

    output_string = ''

    for index in multi_amino_data:
        output_string = output_string + index + ',' + multi_amino_data[index] + '\n'

    return output_string


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

        amino_codon_dict[amino_acid_key] = {
            'num_codons': num_codons,
            'codon_set': codons,
            'score': score
        }
    amino_codon_file.close()
    return amino_codon_dict


def store_amino_codon_dict(amino_codon_dict) -> str:
    # Converts amino_set codon_set correspondences back into a String for roundtrip testing purposes

    output_list = []
    for amino_acid_key in amino_codon_dict:
        num_decodons = amino_codon_dict[amino_acid_key]['num_codons']
        decodons = ''.join(amino_codon_dict[amino_acid_key]['codon_set'])
        score = amino_codon_dict[amino_acid_key]['score']
        output_list.append(f"{num_decodons}{amino_acid_key}{decodons}{score}")

    return ''.join(output_list)


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
        unambiguous_nucleotides = [AMBIGUOUS_UNAMBIGUOUS_NUCLEOTIDE_MAP[nucleotide] for nucleotide in ambiguous_codon]

        # Step 2: forming all possible codons through a cartesian product of nucleotides
        # {'A', 'ATCG', 'GT'} -> {'AAG', 'ATG', 'ACG', 'AGG', 'AAT', 'ATT', 'ACT', 'AGT'}
        unambiguous_codons = [''.join(codon) for codon in itertools.product(*unambiguous_nucleotides)]

        # Step 3: translating codons into amino acids
        # {'AAT', 'ATT', 'ACT', 'AGT'} -> {'K', 'M', 'T', 'R', 'N', 'I', 'S'}
        amino_set.update(set([CODON_TABLE[codon] for codon in unambiguous_codons]))

    return amino_set
