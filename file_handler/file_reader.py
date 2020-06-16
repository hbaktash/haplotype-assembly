import errno
import os
import numpy as np
from statistic_tools import stat_tools
import pickle
from scipy import io

SAM_PATH = "./"
DATA_PATH = os.path.join(os.getcwd(), "data")
DATA_INSTANCE = "outBamBWAMEM_MOT{}_GAPDH"


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def safe_open_w(path, mode):
    ''' Open "path" for writing, creating any parent directories as needed.
    '''
    mkdir_p(os.path.dirname(path))
    return open(path, mode)


def vcf_to_indices(vcf_file):
    indexes = []
    for line in vcf_file.readlines():
        if line[0] == "#":
            continue
        params = line.split()
        variations = params[4]
        if len(variations) <= 5:
            continue
        # if len(variations) >= 9:
        # print("V V V", variations)
        # print(tuple(variations.split(",")[0:-1]))
        # print(params[1])
        pos = int(params[1])
        variations = variations.split(",")[0:-1]
        indexes.append([pos, tuple(variations)])
    return indexes


def vcf_bunch_to_indices(path: str):
    all_indexes = []
    for file in os.listdir(path):
        if file.endswith(".vcf"):
            with open(os.path.join(path, file), 'r') as vcf_file:
                all_indexes += vcf_to_indices(vcf_file)
    return all_indexes


def sam_to_read_array(sam_file, indexes: list, variations=None, verbose=False):  # ignoring mate indexes are one_based
    reads_array = []
    line_count = -1
    for line in sam_file.readlines():
        line_count += 1
        if line[0] == "@":
            continue
        params = line.split()
        pos_one_based = int(params[3])
        read_str = params[9]
        cigar = params[5]
        if (cigar == "*") or (read_str == "*"):  # no alignment
            continue
        parsed_read = parse_CIGAR_read_str(cigar, read_str)
        if parsed_read is None:
            continue
        parsed_read = parsed_read.replace(" ", "X")
        parsed_read = parsed_read.replace("N", "X")
        # print(len(parsed_read))
        start = pos_one_based
        end = pos_one_based + len(parsed_read)
        read_in_indexes = ["X" for _ in range(len(indexes))]
        for i in range(len(indexes)):
            if (indexes[i] < end) and (indexes[i] >= start):
                read_index = indexes[i] - start
                # print("ri:", read_index)
                if parsed_read[read_index] != "X":
                    if not variations[i].__contains__(parsed_read[read_index]) and verbose:
                        print("********* Anomaly", "expected:", variations[i], "\ngot:", parsed_read[read_index])
                        print("     in location:", indexes[i], "line:", line_count)
                    read_in_indexes[i] = parsed_read[read_index]
            elif indexes[i] > end:
                break
        reads_array.append(read_in_indexes)
    return reads_array


def sam_bunch_to_read_array(path: str, indexes: list, variations: list):
    all_reads = []
    i = 0
    for file in os.listdir(path):
        if file.endswith(".sam"):
            print(file, "\n", i + 1)
            i += 1
            with open(os.path.join(path, file), 'r') as sam_file:
                tmp_reads = sam_to_read_array(sam_file, indexes, variations)
                stats = variation_stats(variations, tmp_reads, indexes)
                stat_tools.plot_ranked_hist(stats)
                all_reads += tmp_reads
    bad_rows = check_reads_sanity(all_reads)
    print("bad rows:", bad_rows)
    print("bad indexes:\n", [indexes[bad_row] for bad_row in bad_rows])
    return all_reads


def parse_CIGAR_read_str(cigar: str, read: str):  # will delete the insertions TODO would we get I or D?
    final_align = ""
    # print(cigar)
    i = -1
    if cigar == "*":
        return None
    while cigar != "":
        # print(cigar)
        i += 1
        if cigar[i] == "S":
            size = int(cigar[:i])
            read = read[size:]
            cigar = cigar[i + 1:]
            i = 0
        elif cigar[i] == "M":
            size = int(cigar[:i])
            final_align += read[:size]
            read = read[size:]
            cigar = cigar[i + 1:]
            i = 0
        elif cigar[i] == "I":
            size = int(cigar[:i])
            read = read[size:]
            cigar = cigar[i + 1:]
            i = 0
        elif cigar[i] == "D":
            size = int(cigar[:i])
            final_align += " " * size
            cigar = cigar[i + 1:]
            i = 0
        elif cigar[i] == "H":  # ignore cigar till now
            cigar = cigar[i + 1:]
            i = 0
    return final_align


def check_reads_sanity(reads_list: list, file_name: str = ""):
    bad_cols = []
    for j in range(len(reads_list[0])):
        # print(j)
        flag = True
        for read in reads_list:
            if read[j] != "X":
                flag = False
        if flag:
            bad_cols.append(j)
    if len(bad_cols) != 0:
        print("*****************BAD COLS!!!!!\n", bad_cols)
    return bad_cols


# def single_test(file_num: int):
#     with open(DATA_PATH + "\\vcf\\outBamBWAMEM_MOT{}_GAPDH.vcf".format(file_num), 'r') as vcf_file:
#         index_pairs = vcf_to_indices(vcf_file)
#     print("indexes**********", len(index_pairs))
#     for pair in index_pairs:
#         print(pair)
#
#     indexes = [a[0] for a in index_pairs]
#     with open(DATA_PATH + "\\sam_sorted\\outBamBWAMEM_MOT{}_GAPDH.sam".format(file_num), 'r') as sam_file:
#         reads_list = sam_to_read_array(sam_file, indexes)


def variation_stats(variations, all_read_arrays, indexes: list):
    index_variation_stats = []
    for j in range(len(all_read_arrays[0])):
        varz = {}
        variation = variations[j]
        for read_array in all_read_arrays:
            if read_array[j] != "X":
                if varz.__contains__(read_array[j]):
                    varz[read_array[j]] += 1
                else:
                    varz[read_array[j]] = 1
        # print(varz)
        # print("original:", variation, "\n", indexes[j])
        # print("", end="\n")
        index_variation_stats.append(varz)
    return index_variation_stats


def read_merged():
    merged_vcf_path = os.path.join(DATA_PATH, "merges")
    merged_sam_path = os.path.join(DATA_PATH, "merges")
    with open(os.path.join(merged_vcf_path, "merged-vars.vcf"), 'r') as vcf_file, \
            open(os.path.join(merged_sam_path, "merged-ndup.sam"), 'r') as sam_file:
        index_pairs = vcf_to_indices(vcf_file)
        indexes = [pair[0] for pair in index_pairs]
        variations = [pair[1] for pair in index_pairs]
        read_arrays = sam_to_read_array(sam_file, indexes, variations, verbose=False)
    check_reads_sanity(read_arrays)
    print("num of reads:", len(read_arrays))
    print("variant columns:", len(indexes))
    return read_arrays, index_pairs
    # stats = variation_stats(variations, read_arrays, indexes)
    # stat_tools.plot_ranked_hist(stats)


def test():
    vcf_path = os.path.join(DATA_PATH, "vcf")
    sam_path = os.path.join(DATA_PATH, "sam_sorted")
    all_variant_index_pairs = vcf_bunch_to_indices(vcf_path)
    indexes = [pair[0] for pair in all_variant_index_pairs]
    variations = [pair[1] for pair in all_variant_index_pairs]
    all_read_arrays = sam_bunch_to_read_array(sam_path, indexes, variations)
    print("rows", len(all_read_arrays), "cols", len(indexes))
    cols = len(indexes)
    # print("cols,", cols)
    count = 0
    # variation_stats(variations, all_read_arrays)

    # for read_array in all_read_arrays:
    #     for r in read_array:
    #         if r == "X":
    #             count += 1
    # print(count, " out of ", rows*len(all_read_arrays))


def save_indexes(indexes: list, separated: bool = False):
    if not separated:
        save_as_file(indexes, "all-indexes")
    else:
        i = 0
        for index_list in indexes:
            i += 1
            save_as_file(index_list, "indexes-part-{}".format(i))


def save_as_file(obj, name: str, version: str = "", part_number: int = -1):
    if version == "":
        with safe_open_w(os.path.join(DATA_PATH, "results", name + ".pkl"), "wb+") as file:
            pickle.dump(obj, file)
        return
    if part_number != -1:
        version_part = os.path.join(version, "part_{}".format(str(part_number)))
    else:
        version_part = version
    with safe_open_w(os.path.join(DATA_PATH, "results", version_part, name + ".pkl"), "wb+") as file:
        pickle.dump(obj, file)


def save_parts_as_mat(version: str = "", part_number: int = -1):
    version_part = os.path.join(version, "part_{}".format(str(part_number)))
    with safe_open_w(os.path.join(DATA_PATH, "results", version_part, "completed-matrix.pkl"), 'rb') as file:
        complete_mat = pickle.load(file)
    io.savemat(os.path.join(DATA_PATH, "results", "exon-mat{}.mat".format(part_number)), {'mat{}'.format(part_number): complete_mat})

def save_result(completed_mat, reads_arr, mapping_dicts, version: str, part_number: int = -1):
    if part_number != -1:
        version_part = os.path.join(version, "part_{}".format(str(part_number)))
    else:
        version_part = version
    with safe_open_w(os.path.join(DATA_PATH, "results", version_part, "completed-matrix" + ".pkl"), "wb+") as completed_mat_file:
        pickle.dump(completed_mat, completed_mat_file)
    with safe_open_w(os.path.join(DATA_PATH, "results", version_part, "read-arrays" + ".pkl"), "wb+") as read_arr_file:
        pickle.dump(reads_arr, read_arr_file)
    with safe_open_w(os.path.join(DATA_PATH, "results", version_part, "mapping-dict" + ".pkl"), "wb+") as mapping_dicts_file:
        pickle.dump(mapping_dicts, mapping_dicts_file)


def load_results(version: str, part_number: int = -1):
    if part_number != -1:
        version_part = os.path.join(version, "part_{}".format(str(part_number)))
    else:
        version_part = version
    with open(os.path.join(DATA_PATH, "results", version_part, "completed-matrix.pkl"), 'rb') as completed_mat_file:
        completed_matrix = pickle.load(completed_mat_file)
    with open(os.path.join(DATA_PATH, "results", version_part, "read-arrays.pkl"), 'rb') as read_arr_file:
        read_arrays = pickle.load(read_arr_file)
    with open(os.path.join(DATA_PATH, "results", version_part, "mapping-dict.pkl"), 'rb') as mapping_dict_file:
        mapping_dicts = pickle.load(mapping_dict_file)
    return completed_matrix, read_arrays, mapping_dicts


if __name__ == '__main__':
    a = [1, 2, 3, 4]
    save_as_file(a, "test")
