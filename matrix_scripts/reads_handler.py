import numpy as np


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


def separate_reads(read_arrays: list, index_pairs: list, separation_indexes: list):
    """returns separated read arrays with independent sizes
    each element is a read_arrays list"""

    print("extracting separation positions")
    separation_positions = []
    indexes = [a[0] for a in index_pairs]
    for sep_index in separation_indexes:
        for i in range(len(indexes)):
            if sep_index == indexes[i]:
                separation_positions.append(i)
                break
            elif indexes[i] < sep_index < indexes[i + 1]:
                separation_positions.append(i + 1)
                break

    print("separating reads")
    sub_read_arrays = []
    sub_index_pairs = []
    first_position = 0
    separation_positions.append(len(indexes))
    for i in range(len(separation_positions)):
        tmp_sub_read_array = [read[first_position:separation_positions[i]] for read in read_arrays]
        sub_read_arrays.append(tmp_sub_read_array.copy())
        sub_index_pairs.append(index_pairs[first_position:separation_positions[i]].copy())
        first_position = separation_positions[i]
    print("removing empty reads")
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


def read_array_to_matrix(read_arrays: list, index_pairs: list, difference_magnitude=1j):
    new_variations = build_variations(read_arrays)
    var_map = {"A": 1, "C": -1, "T": 1j, "G": -1j, "X": 0}
    variation_to_value_dicts = []
    for variation in new_variations:
        # print(variation)
        tmp_dict = {}
        i = 1
        for var in variation:
            tmp_dict[var] = var_map[var]
            # if var == "X":
            #     tmp_dict[var] = 0
            #     continue
            # tmp_dict[var] = i * difference_magnitude
            # if difference_magnitude.real != 0:  # if complex
            #     print("LINEAR!!!!!!")
            #     i += 1
            # else:
            #     i *= difference_magnitude
        # print(tmp_dict)
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


def read_mat_to_arrays(read_mat: np.ndarray, dicts: list):
    if len(dicts) != read_mat.shape[1]:
        print("dict and array size don't match")
    read_arrays = []
    for i in range(read_mat.shape[0]):
        read_arrays.append([
            [key for (key, value) in dicts[j].items() if value == read_mat[i, j]][0]
            for j in range(len(dicts))])
    return read_arrays
