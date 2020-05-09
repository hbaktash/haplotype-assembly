import os
import numpy as np

SAM_PATH = "./"
DATA_PATH = 'C:\\Users\\hosse\\Desktop\\Internship\\code\\data\\'
DATA_INSTANCE = "outBamBWAMEM_MOT{}_GAPDH.bam"


def vcf_to_indices(vcf_file):
    indexes = []
    for line in vcf_file.readlines():
        if line[0] == "#":
            continue
        params = line.split()
        alleles = params[4]
        if len(alleles) <= 5:
            continue
        pos = int(params[1])
        x_1, x_2 = alleles.split(",")[0], alleles.split(",")[1]
        indexes.append([pos, (x_1, x_2)])
    return indexes


def sam_to_read_array(sam_file, indexes: list):  # ignoring mate indexes are one_based
    reads_array = []
    for line in sam_file.readlines():
        if line[0] == "@":
            continue
        params = line.split()
        pos_one_based = int(params[3])
        read_str = params[9]
        cigar = params[5]
        parsed_read = parse_CIGAR_read_str(cigar, read_str)

        start = pos_one_based
        end = pos_one_based + len(parsed_read)
        read_in_indexes = ["X" for i in range(len(indexes))]
        for i in range(len(indexes)):
            if (indexes[i] <= end) and (indexes[i] >= start):
                read_index = indexes[i] - start
                read_in_indexes[i] = parsed_read[read_index]
            elif indexes[i] > end:
                break
        reads_array.append(read_in_indexes)
    return reads_array


def parse_CIGAR_read_str(cigar: str, read: str):  # will delete the insertions TODO would we get I or D?
    final_align = ""
    for i in range(len(cigar)):
        if cigar[i] == "S":
            size = int(cigar[:i])
            read = read[size:]
            cigar = cigar[i + 1:]
        elif cigar[i] == "M":
            size = int(cigar[:i])
            final_align += read[:size]
            read = read[size:]
            cigar = cigar[i + 1:]
    return final_align


def test(file_num):
    with open(DATA_PATH+"\\vcf\\outBamBWAMEM_MOT{}_GAPDH.vcf".format(file_num), 'r') as vcf_file:
        indice_pairs = vcf_to_indices(vcf_file)
        print(indice_pairs[:10])


if __name__ == '__main__':
    test(101)
