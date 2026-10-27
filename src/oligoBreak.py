# #!/user/bin/env python 3

import data_parser as dp
import numpy as np
from typing import NamedTuple

def OligoBreak(min_oligo_len, max_oligo_len, min_temp, max_temp, dna_sequence, position_overlap, nt_cost, position_degeneracy_dict):
    
    sequence_len = len(dna_sequence)
    temp_range = max_temp - min_temp + 1

    backtracking_arr = {temp: [None] * sequence_len for temp in range(min_temp, max_temp + 1)}

    Backtrack = NamedTuple('Backtrack', [('parent', int), ('overlap', int), ('parent_temp', int)])

    # initialize the cost array, which will be used to store the minimum cost of synthesizing 
    # the DNA sequence from position 0 to position i at temp t
    # init to infinity
    oligo_cost_arr = {temp: [np.inf] * sequence_len for temp in range(min_temp, max_temp + 1)}

    for temp in range(min_temp, max_temp + 1):

        backtracking_arr[temp][0] = Backtrack(parent=None, overlap=0, parent_temp=None)

        oligo_cost_arr[temp][0] = 0

        for oligo_len in range(min_oligo_len, max_oligo_len + 1):
            oligo_cost = dp.calculate_cost(0, oligo_len, nt_cost, position_degeneracy_dict)
            if oligo_cost is not None:
  
                oligo_cost_arr[temp][oligo_len] = oligo_cost
                backtracking_arr[temp][oligo_len] = Backtrack(parent=None, overlap=0, parent_temp=temp)

        for position in range(min_oligo_len, sequence_len):

            prev_oligo_cost = oligo_cost_arr[temp][position - 1]

            if prev_oligo_cost != np.inf:

                # get overlap distance
                overlap = position_overlap[temp][position]

                #if overlap is None and temp > min_temp:
                    #overlap = position_overlap[temp - 1][position]
                    #parent_temp = temp - 1
                #else:
                    #parent_temp = temp

                if overlap is None:
                    continue

                for oligo_len in range(min_oligo_len, max_oligo_len + 1):

                    # start and end of oligo at a certain position with a certain overlap length
                    oligo_start = position - overlap
                    oligo_end = position + oligo_len - overlap - 1
                    
                    # check if oligo end is within the bounds of the sequence, if not, skip
                    if oligo_end < sequence_len:
                        
                        # calculate the cost of synthesizing the oligo at the current position
                        new_cost = dp.calculate_cost(oligo_start, oligo_end, nt_cost, position_degeneracy_dict)

                        if new_cost is None:
                            continue

                        current_oligo_cost = prev_oligo_cost + new_cost

                        # if current cost is less then previous cost, replace the cost
                        if current_oligo_cost < oligo_cost_arr[temp][oligo_end]:
                            oligo_cost_arr[temp][oligo_end] = current_oligo_cost
                            backtracking_arr[temp][oligo_end] = Backtrack(parent=position - 1, overlap=overlap, parent_temp=temp)

    return oligo_cost_arr, backtracking_arr

if __name__ == "__main__":

    (
        protein, 
        codon_frequencies, 
        amino_codon_dict, 
        multi_amino_dict, 
        oligo_output_file, 
        min_oligo_len, 
        max_oligo_len, 
        min_temp, 
        max_temp, 
        max_temp_diff, 
        nt_cost 
    ) = dp.parse_arguments()

    dna_seq_normal = dp.convert_protein_sequence_to_dna_sequence(protein.sequence, codon_frequencies)
    # TEST
    # print(f"DNA sequence: {dna_seq_normal}")
    min_overlap, max_overlap = 1, 700
    dna_seq_binary = dp.convert_dna_sequence_to_binary_sequence(dna_seq_normal)
    seq_position_overlap = dp.calculate_position_overlap(dna_seq_binary, min_temp, max_temp, min_overlap, max_overlap)
    multi_codon_dict = dp.convert_multi_amino_to_multi_codon_data(multi_amino_dict, amino_codon_dict)
    position_degeneracy = dp.calculate_position_degeneracy(multi_codon_dict)

    oligo_cost_arr, backtracking_arr = OligoBreak(min_oligo_len, 
                                                  max_oligo_len, 
                                                  min_temp, 
                                                  max_temp, 
                                                  dna_seq_binary, 
                                                  seq_position_overlap, 
                                                  nt_cost, 
                                                  position_degeneracy)

    optimal_path_dict = dp.get_optimal_path(backtracking_arr, oligo_cost_arr, len(dna_seq_binary))

    oligos = dp.write_oligos_to_fasta(dna_seq_normal, optimal_path_dict, multi_codon_dict, oligo_output_file)

