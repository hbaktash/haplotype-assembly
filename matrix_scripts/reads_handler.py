import numpy as np
import collections


def remove_empty_reads(sub_reads: list):
    init_size = len(sub_reads)
    non_empty_indexes = []
    i = 0
    for read in sub_reads:
        if not all(e == "X" for e in read):
            non_empty_indexes.append(i)
        i += 1
    pruned_reads = [sub_reads[i].copy() for i in non_empty_indexes]
    print("after {}, before {}".format(len(pruned_reads), init_size))
    return pruned_reads


def remove_repeated_reads(reads: list):
    deletables = []
    for i in range(len(reads)):
        for k in range(0,i):
            if collections.Counter(reads[i]) == collections.Counter(reads[k]):
                deletables.append(i)
                break
    print("{} repeated reads deleted".format(len(deletables)))
    return [reads[i] for i in range(len(reads)) if not deletables.__contains__(i)], deletables


def check_read_inclusion(big_r: list, small_r: list):
    if len(big_r) != len(small_r):
        print("WTF")
        return False
    for i in range(len(big_r)):
        if (big_r[i] != "X") and (small_r[i] != "X"):
            if big_r[i] != small_r[i]:
                return False
            else:
                continue
        elif small_r[i] != "X":
            return False
        else:
            continue
    return True


def remove_submissive_reads(reads: list):
    deletables = []
    for i in range(len(reads)):
        for k in range(len(reads)):
            if i != k:
                if check_read_inclusion(big_r=reads[k], small_r=reads[i]):
                    deletables.append(i)
                    break
    print("{} submissive reads deleted".format(len(deletables)))
    return [reads[i] for i in range(len(reads)) if not deletables.__contains__(i)], deletables


def find_position_in_sorted_indexes(position: int, indexes: list):
    if position < indexes[0]:
        return 0
    for i in range(len(indexes)):
        if position == indexes[i]:
            return i
        elif indexes[i] < position < indexes[i + 1]:
            return i + 1
    return len(indexes)


def extract_separation_positions_in_my_indexes(index_pairs: list, separation_indexes: list):
    separation_positions = []
    indexes = [a[0] for a in index_pairs]
    for sep_index in separation_indexes:
        separation_positions.append(find_position_in_sorted_indexes(sep_index, indexes))
    return separation_positions


def separate_reads(read_arrays: list, index_pairs: list, separation_indexes: list):
    """returns separated read arrays with independent sizes
    each element is a read_arrays list"""

    print("    -extracting separation positions")
    separation_positions = extract_separation_positions_in_my_indexes(index_pairs, separation_indexes)
    indexes = [a[0] for a in index_pairs]
    print("    -separating reads")
    sub_read_arrays = []
    sub_index_pairs = []
    first_position = 0
    separation_positions.append(len(indexes))
    for i in range(len(separation_positions)):
        tmp_sub_read_array = [read[first_position:separation_positions[i]] for read in read_arrays]
        sub_read_arrays.append(tmp_sub_read_array.copy())
        sub_index_pairs.append(index_pairs[first_position:separation_positions[i]].copy())
        first_position = separation_positions[i]
    print("    -removing empty reads")
    independent_read_arrays = [remove_empty_reads(sub_read_array) for sub_read_array in sub_read_arrays]
    return independent_read_arrays, sub_index_pairs


def separate_exons(read_arrays: list, index_pairs: list, exon_intervals: list):
    print("     -extracting interval positions")
    indexes = [a[0] for a in index_pairs]
    interval_positions = [(find_position_in_sorted_indexes(p[0], indexes),
                           find_position_in_sorted_indexes(p[1], indexes))
                          for p in exon_intervals]

    print("     -separating reads")
    sub_read_arrays = []
    sub_index_pairs = []
    for interval_position in interval_positions:
        tmp_sub_read_array = [read[interval_position[0]:interval_position[1]] for read in read_arrays]
        sub_read_arrays.append(tmp_sub_read_array.copy())

        sub_index_pairs.append(index_pairs[interval_position[0]:interval_position[1]].copy())
    print("    -removing empty reads")
    independent_read_arrays = [remove_empty_reads(sub_read_array) for sub_read_array in sub_read_arrays]
    return independent_read_arrays, sub_index_pairs


def build_variations(read_arrays: list):
    variations = []
    for j in range(len(read_arrays[0])):
        col_vars = []
        for read in read_arrays:
            if not col_vars.__contains__(read[j]):
                col_vars.append(read[j])
        variations.append(col_vars.copy())
    return variations


def read_array_to_matrix(read_arrays: list):
    new_variations = build_variations(read_arrays)
    var_map = {"A": 1, "C": -1, "T": 1j, "G": -1j, "X": 0}
    variation_to_value_dicts = []
    for variation in new_variations:
        # print(variation)
        tmp_dict = {}
        i = 1
        for var in variation:
            tmp_dict[var] = var_map[var]
        variation_to_value_dicts.append(tmp_dict.copy())
    rows = len(read_arrays)
    cols = len(read_arrays[0])
    # print(read_arrays[0])
    final_matrix = np.zeros((rows, cols), dtype=np.complex)
    for (j, vars_dict) in zip(range(cols), variation_to_value_dicts):
        for (read, i) in zip(read_arrays, range(rows)):
            final_matrix[i, j] = vars_dict[read[j]]
    # print(final_matrix[0, :])
    return final_matrix, variation_to_value_dicts


def read_mat_to_arrays(read_mat: np.ndarray, dicts: list = None, use_dicts: bool = True):
    if use_dicts and dicts != None:
        if len(dicts) != read_mat.shape[1]:
            print("dict and array size don't match")
            return
    read_arrays = []
    reverse_var_map = {1: "A", -1: "C", 1j: "T", -1j: "G", 0: "X"}
    for i in range(read_mat.shape[0]):
        if use_dicts:
            read_arrays.append([
                [key for (key, value) in dicts[j].items() if value == read_mat[i, j]][0]
                for j in range(len(dicts))])
        else:
            read_arrays.append([reverse_var_map[read_mat[i, j]] for j in range(read_mat.shape[1])])
    return read_arrays


def get_block_frequencies(blocks: list, read_arrays: list):
    freqs = []
    for block in blocks:
        count = 0
        for read in read_arrays:
            if check_read_inclusion(block, read) or check_read_inclusion(read, block):
                count += 1
        freqs.append(count)
    return freqs


def map_blocks_to_source_blocks(source_blocks: list, individual_blocks: list):
    mapping = []
    # print("SS")
    # for b in source_blocks:
    #     print(b)
    # print("indi")
    # for b in individual_blocks:
    #     print(b)
    for ind_block in individual_blocks:
        diffs = [blocks_distance(ind_block, source_block) for source_block in source_blocks]
        min_index = diffs.index(min(diffs))
        mapping.append(min_index)
    return mapping


def blocks_distance(b1: list, b2: list):
    if len(b1) != len(b2):
        return -1
    diff = 0
    for e1, e2 in zip(b1, b2):
        if e1 != e2:
            diff += 1
        if (e1 == "X") or (e2 == "X"):  # if either is X then assume 1/2 chance for equality
            diff -= 0.5
    return diff
