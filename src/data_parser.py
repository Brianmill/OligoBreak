import random
import itertools
import sys
from math import floor
from typing import NamedTuple

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
    # The stop codon can be coded as 010.
    # This has little impact because it is a stop codon.
    '*': '001',
}

# Mapping from ambiguous DNA nucleotides
# to unambiguous DNA nucleotides.
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

class protein(NamedTuple):
    name: str
    sequence: str

def parse_arguments():
    """
    Parse command-line arguments and load the referenced input data.

    Returns
    -------
    tuple
        A tuple containing the protein record, organism codon-frequency
        dictionary, amino-acid-to-codon dictionary, multi-amino dictionary,
        output path, minimum and maximum oligo lengths, minimum and maximum
        temperatures, maximum temperature difference, and cost per base.
    """
    args = sys.argv

    if len(args) < 6:
        raise ValueError("Not enough arguments provided")

    protein_sequence = args[1]
    organism_codon_frequency = args[2]
    amino_codon_data = args[3]
    multi_amino_data = args[4]
    oligo_output = args[5]

    min_oligo_len = None
    max_oligo_len = None
    min_temp = None
    max_temp = None
    max_temp_diff = None
    cost_per_base = None

    i = 6
    while i < len(args):
        if args[i] == "-nl":
            min_oligo_len = int(args[i + 1])
            i += 2
        elif args[i] == "-xl":
            max_oligo_len = int(args[i + 1])
            i += 2
        elif args[i] == "-nt":
            min_temp = int(args[i + 1])
            i += 2
        elif args[i] == "-xt":
            max_temp = int(args[i + 1])
            i += 2
        elif args[i] == "-d":
            max_temp_diff = int(args[i + 1])
            i += 2
        elif args[i] == "-c":
            cost_per_base = float(args[i + 1])
            i += 2
        else:
            raise ValueError(f"Unknown argument: {args[i]}")

    protein_sequence_file = extract_protein_sequence(protein_sequence)
    organism_codon_frequency_file = extract_organism_codon_frequencies(organism_codon_frequency)
    amino_codon_data_file = extract_amino_codon_dict(amino_codon_data)
    multi_amino_data_file = extract_multi_amino_data(multi_amino_data)
    oligo_output_file = oligo_output

    return (
        protein_sequence_file, 
        organism_codon_frequency_file, 
        amino_codon_data_file, 
        multi_amino_data_file, 
        oligo_output_file, 
        min_oligo_len, 
        max_oligo_len, 
        min_temp, 
        max_temp, 
        max_temp_diff, 
        cost_per_base
    )

def convert_dna_sequence_to_binary_sequence(dna_seq):
    """
    Convert a DNA sequence to binary, using 1 for C/G and 0 for A/T.
    """
    binary_seq = ""

    for nucleotide in dna_seq:
        if nucleotide in ['C', 'G']:
            binary_seq = binary_seq + "1"
        else:
            binary_seq = binary_seq + "0"

    return binary_seq

def calculate_final_order_cost(cost_per_base, oligos, position_degeneracy_dict):
    """
    Calculate the total synthesis cost of a segmented DNA sequence.

    Each oligo's base cost is multiplied by the product of degeneracies for
    amino-acid positions that overlap that oligo.
    """
    total_cost = 0

    for oligo in oligos.keys():
        start, end = oligo
        oligo_len = end - start + 1
        oligo_cost = oligo_len * cost_per_base

        oligo_position_degeneracy_set = set()

        for position in position_degeneracy_dict.keys():
            position = int(position)

            nt_positions = [(position * 3) - 2, (position * 3) - 1, position * 3]

            if any(start <= nt_position <= end and nt_position >= 0 for nt_position in nt_positions):
                oligo_position_degeneracy_set.add(position_degeneracy_dict[str(position)])

        oligo_degeneracy = 1

        for position_degeneracy in oligo_position_degeneracy_set:
            oligo_degeneracy = oligo_degeneracy * position_degeneracy

        total_oligo_cost = oligo_cost * oligo_degeneracy

        total_cost += total_oligo_cost

    return total_cost

def write_oligos_to_fasta(dna_sequence, optimal_path_dict, multi_codon_dict, multi_fasta_file):
    """
    Write the segmented sequence and its degenerate variants to a FASTA file.

    Returns
    -------
    dict
        Mapping of ``(start, end)`` nucleotide-bound tuples to lists of oligo
        sequence variants. Bounds use the same coordinates used to write the
        FASTA headers.
    """
    
    oligos = {}

    for breakpoint, oligo_len in optimal_path_dict.items():

        oligo_start = breakpoint - oligo_len
        oligo_end = breakpoint
        base_seq = dna_sequence[oligo_start:oligo_end]
        oligo_variants = [base_seq]

        for position in multi_codon_dict.keys():
            position = int(position)

            nt_positions = [(position * 3) - 2, (position * 3) - 1, position * 3]

            if any(oligo_start <= nt_position <= oligo_end for nt_position in nt_positions):
                multi_codons = multi_codon_dict[str(position)]
                new_variants = []
                for seq in oligo_variants:

                    start_idx = nt_positions[0] - oligo_start
                    end_idx = nt_positions[2] - oligo_start + 1

                    for degen_codon in multi_codons:
                        new_seq = (seq[:start_idx] + degen_codon + seq[end_idx:])

                        new_variants.append(new_seq)

                oligo_variants = new_variants

        oligos[(oligo_start, oligo_end)] = oligo_variants

    with open(multi_fasta_file, "w") as f:
        for i, seq in enumerate(oligos.keys(), start=1):
            start, end = seq
            oligo = oligos[seq]
            for oligo_varient in oligo:
                header = f">oligo_{i}, positions {start}-{end}"
                f.write(header + "\n")
                f.write(oligo_varient + "\n")

    return oligos

def get_optimal_path(backtracking_arr, oligo_cost_arr, seq_len):
    """
    Reconstruct the lowest-cost oligo path from the backtracking arrays.

    Each backtracking entry stores the previous oligo's index and overlap.

    Returns
    -------
    dict
        Mapping of each oligo's ending sequence position to its length.
    """

    optimal_path = {}

    optimal_temp = None
    for temp in backtracking_arr:
        if backtracking_arr[temp][seq_len - 1] is not None:
            current_cost = oligo_cost_arr[temp][seq_len - 1]
            if optimal_temp is None or current_cost < oligo_cost_arr[optimal_temp][seq_len - 1]:
                optimal_temp = temp

    current_index = len(backtracking_arr[temp]) - 1

    while current_index is not None:
        if backtracking_arr[optimal_temp][current_index].parent is not None:
            parent = backtracking_arr[optimal_temp][current_index].parent + 1
            overlap = backtracking_arr[optimal_temp][current_index].overlap
        if (parent != current_index + 1):
            optimal_path[current_index + 1] = current_index - parent + overlap + 1
        else:
            optimal_path[current_index + 1] = current_index + 1

        current_index = backtracking_arr[optimal_temp][current_index].parent

    optimal_path = dict(reversed(optimal_path.items()))
    return optimal_path

def calculate_position_degeneracy(multi_codon_dict):
    """
    Create a dictionary that maps each degenerate position to the number
    of degenerate codons required at the position.

    Returns
    -------
    dict
        Mapping of amino-acid position keys to the number of codons in the
        corresponding codon set.
    """
    position_degeneracy_dict = {}

    position_degeneracy_dict = {
        position: len(codon_set)
        for position, codon_set in multi_codon_dict.items()
    }

    return position_degeneracy_dict

def calculate_position_overlap(dna_binary_sequence, min_temp, max_temp, min_overlap, max_overlap):
    """
    Find an overlap length for each sequence position and target temperature.

    Returns
    -------
    dict
        Mapping of integer temperatures to lists indexed by nucleotide
        position. Each value is the matching overlap length, or None when no
        permitted overlap reaches that temperature at that position.
    """
    sequence_len = len(dna_binary_sequence)

    overlap_dict = {temp: [None] * sequence_len for temp in range(min_temp, max_temp + 1)}

    for temp in range(min_temp, max_temp + 1):

        for dna_index, binary_nucleotide in enumerate(dna_binary_sequence):

            cg_amount = 0

            for overlap_length in range(min_overlap, max_overlap + 1):

                overlap_index = dna_index - overlap_length + 1
                
                if overlap_index < 0 or overlap_index >= sequence_len:
                    break
            
                cg_amount += int(dna_binary_sequence[overlap_index])

                melting_temperature = round(
                    calculate_dna_melting_temp(cg_amount, overlap_length)
                )

                if melting_temperature == temp:
                    overlap_dict[temp][dna_index] = overlap_length
                    break

    return overlap_dict

def calculate_cost(start, end, single_nucleotide_cost, position_degeneracy_dict):
    """
    Calculate synthesis cost for an oligo covering nucleotide indices start
    through end, inclusive.

    Returns None if any degenerate amino-acid position (3-nt codon) would be
    partially included by the oligo (i.e. some of the 3 nt are inside and some
    outside). Degenerate positions fully inside the oligo contribute their
    degeneracy multiplicatively to the oligo cost. Malformed degeneracy values
    are treated as 1.
    """
    if end <= start:
        return None

    oligo_length = end - start + 1
    total_cost = oligo_length * single_nucleotide_cost

    position_degeneracy_set = []

    if not position_degeneracy_dict:
        return total_cost

    for pos_key, degeneracy_val in position_degeneracy_dict.items():
        try:
            aa_pos = int(pos_key)
        except Exception:
            continue

        nt0 = (aa_pos - 1) * 3
        nt_positions = (nt0, nt0 + 1, nt0 + 2)

        inside_flags = [start <= nt <= end for nt in nt_positions]

        if not any(inside_flags):
            continue

        if not all(inside_flags):
            return None

        try:
            position_degeneracy_set.append(int(degeneracy_val))
        except Exception:
            position_degeneracy_set.append(1)

    total_oligo_degeneracy = 1
    for position_degeneracy in position_degeneracy_set:
        total_oligo_degeneracy *= position_degeneracy

    total_cost *= total_oligo_degeneracy

    return total_cost

def extract_protein_sequence(protein_seq_file_path: str):
    """
    Extract the protein name and sequence from a FASTA file.

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

    this_protein = protein(name=protein_name, sequence=protein_seq)
    return this_protein


def extract_multi_amino_data(multi_amino_file_path: str):
    """
    Extract multi-amino positions and their amino-acid alternatives.

    Multi amino data is locations in the protein sequence where more than one amino
    acid is desired. This is indicated through the use of index, amino list pairs.

    Parameters
    ----------
    multi_amino_file_path : str
        the path of the file containing the multi amino data

    Returns
    -------
    multi_amino_data : dict
        Mapping of position strings to tuples of allowed amino-acid symbols:
        ``multi_amino_data[position] = (amino_acid, ...)``.
    """
    multi_amino_file = open(multi_amino_file_path, "r")
    multi_amino_data = {}

    for line in multi_amino_file:
        index, amino_acids = line.strip().split(',')
        multi_amino_data[index] = tuple(amino_acids)

    multi_amino_file.close()
    return multi_amino_data


def extract_organism_codon_frequencies(codon_frequencies_file_path: str):
    """
    Extract organism codon frequencies into a dictionary.

    Parameters
    ----------
    codon_frequencies_file_path : str
        Path to the file containing the organism codon frequencies.

    Returns
    -------
    dict
        Mapping of codon strings to their frequency values as read from the
        file, stored in permille (per thousand):
        ``codon_per_mille_data[codon] = per_mille``.
    """
    codon_frequency_file = open(codon_frequencies_file_path, "r")
    codon_per_mille_data = {}

    for line in codon_frequency_file:
        codon, per_mille = line.strip().split(',')
        codon_per_mille_data[codon] = per_mille
    codon_frequency_file.close()
    return codon_per_mille_data


def extract_amino_codon_dict(amino_codon_file_path: str):
    """
    Extract mappings between amino-acid sets and corresponding codon sets.

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
            Mapping of 20-character amino-acid binary keys to records with
            ``num_codons`` (string), ``codon_set`` (ordered tuple of codons),
            and ``score`` (string):
            ``amino_codon_dict[key] = {'num_codons': ..., 'codon_set': ..., 'score': ...}``.

    Codon sets are stored as tuples to preserve order in round-trip tests.
    """
    amino_codon_file = open(amino_codon_file_path, "r")
    amino_codon_dict = {}

    for line in amino_codon_file:
        num_codons, amino_acid_key, codons, score = line.strip().split('\t')
        codons = tuple(codons.split())

        """
        Amino acid key has the format of a 20 digit binary number.
        i.e. {'V', 'Y', 'T', 'H', 'A'} -> 11010000000100000001
        The logic is as follows:
        VYWTSPFMKLIHGEQCDNRA   List of amino acids in reverse alphabetical order, according to full name
        ||.|.......|.......|
        11010000000100000001
        """

        amino_codon_dict[amino_acid_key] = {
            'num_codons': num_codons,
            'codon_set': codons,
            'score': score
        }
    amino_codon_file.close()
    return amino_codon_dict


def store_protein_sequence(protein_sequence_data) -> str:
    """
    Serialize protein name and sequence data for round-trip testing.
    """
    name = protein_sequence_data.name
    sequence = protein_sequence_data.sequence
    fasta_string = name + sequence
    return fasta_string


def store_multi_amino_data(multi_amino_data) -> str:
    """
    Serialize multi-amino position data for round-trip testing.
    """
    output_string = ''

    for index in multi_amino_data:
        multi_amino_string = ''.join(str(amino) for amino in multi_amino_data[index])
        output_string = output_string + index + ',' + multi_amino_string + '\n'

    return output_string


def store_amino_codon_dict(amino_codon_dict) -> str:
    """
    Serialize amino-acid-to-codon mappings for round-trip testing.
    """
    output_list = []
    for amino_acid_key in amino_codon_dict:
        num_decodons = amino_codon_dict[amino_acid_key]['num_codons']
        decodons = ''.join(amino_codon_dict[amino_acid_key]['codon_set'])
        score = amino_codon_dict[amino_acid_key]['score']
        output_list.append(f"{num_decodons}{amino_acid_key}{decodons}{score}")

    return ''.join(output_list)


def store_organism_codon_frequencies(codon_per_mille_data) -> str:
    """
    Serialize organism codon frequencies for round-trip testing.
    """
    output_list = []
    for amino_acid in codon_per_mille_data:
        per_mille = codon_per_mille_data[amino_acid]
        output_list.append(f"{amino_acid},{per_mille}")
    return '\n'.join(output_list)


def convert_amino_set_to_codon_set(amino_set, aminos_to_codons_mapping):
    """
    Look up the codon set corresponding to a set of amino acids.

    The amino-acid set is encoded as a binary key in the supplied mapping.
    """
    if not amino_set:
        return set()

    amino_list = "VYWTSPFMKLIHGEQCDNRA"
    binary_representation = ""
    for index, amino in enumerate(amino_list):
        if amino in amino_set:
            binary_representation = binary_representation + "1"
        else:
            binary_representation = binary_representation + "0"

    codon_set = aminos_to_codons_mapping[binary_representation]['codon_set']
    return codon_set


def convert_codon_set_to_amino_set(ambiguous_codons):
    """
    Expand ambiguous codons and return the amino acids they can encode.
    """
    amino_set = set()
    for ambiguous_codon in ambiguous_codons:
        unambiguous_nucleotides = (AMBIGUOUS_UNAMBIGUOUS_NUCLEOTIDE_MAP[nucleotide] for nucleotide in ambiguous_codon)

        unambiguous_codons = (''.join(codon) for codon in itertools.product(*unambiguous_nucleotides))

        amino_set.update(set(CODON_TABLE[codon] for codon in unambiguous_codons))

    return amino_set


def convert_multi_amino_to_multi_codon_data(multi_amino_data, aminos_to_codons_mapping):
    """
    Convert amino-acid alternatives at each position to codon alternatives.

    Returns
    -------
    dict
        Mapping of position keys from ``multi_amino_data`` to sets of
        corresponding ambiguous codon strings.
    """
    multi_codon_data = {}
    for position, aminos in multi_amino_data.items():
        codons = convert_amino_set_to_codon_set(aminos, aminos_to_codons_mapping)
        multi_codon_data[position] = codons
    return multi_codon_data


def convert_multi_codon_to_multi_amino_data(multi_codon_data):
    """
    Convert codon alternatives at each position to amino-acid alternatives.

    Returns
    -------
    dict
        Mapping of position keys from ``multi_codon_data`` to sets of amino
        acid symbols encoded by those codons.
    """
    multi_amino_data = {}
    for position, codons in multi_codon_data.items():
        aminos = convert_codon_set_to_amino_set(codons)
        multi_amino_data[position] = aminos
    return multi_amino_data

def convert_protein_sequence_to_dna_sequence(protein_sequence, codon_frequencies):
    """
    Generate a DNA sequence by assigning codons according to organism usage.

    Codons are apportioned by amino-acid frequency, then selected randomly
    without replacement from each amino acid's allocation.
    """
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
    """
    Apportion each amino acid across synonymous codons by codon frequency.

    Uses the Hamilton method of apportionment to assign integer codon counts.

    Returns
    -------
    dict
        Mapping of amino-acid symbols to dictionaries mapping synonymous
        codons to their allocated integer counts:
        ``codon_allocations[amino_acid][codon] = count``.
    """

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
    """
    Translate a DNA sequence into its protein sequence.
    """
    codon_list = [dna_sequence[i:i + 3] for i in range(0, len(dna_sequence), 3)]
    amino_list = [CODON_TABLE[codon] for codon in codon_list]
    protein_sequence = ''.join(amino_list)
    return protein_sequence


def convert_protein_to_max_cg_dna_binary_sequence(protein_sequence):
    """
    Return the binary DNA representation with maximum C/G content.
    """
    codon_list = [CODON_WITH_HIGHEST_CG_CONTENT[amino] for amino in protein_sequence]
    dna_sequence = ''.join(codon_list)
    return dna_sequence


def convert_protein_to_min_cg_dna_binary_sequence(protein_sequence):
    """
    Return the binary DNA representation with minimum C/G content.
    """
    codon_list = [CODON_WITH_LOWEST_CG_CONTENT[amino] for amino in protein_sequence]
    dna_sequence = ''.join(codon_list)
    return dna_sequence


def calculate_dna_melting_temp(cg_amount, dna_length):
    """
    Estimate DNA melting temperature from sequence length and C/G count.

    Uses the short-oligo rule below 14 bases and the long-oligo equation for
    sequences of at least 14 bases.
    """
    if dna_length < 0:
        raise ValueError("Error: dna length cannot be negative")
    if dna_length == 0:
        raise ZeroDivisionError("Error: Division by zero is not allowed. Please provide a non-zero denominator.")
    if dna_length < 14:
        at_amount = dna_length - cg_amount
        melting_temperature = at_amount * 2 + cg_amount * 4
        return melting_temperature
    if dna_length >= 14:
        melting_temperature = 64.9 + 41 * (cg_amount - 16.4) / dna_length
        return melting_temperature
