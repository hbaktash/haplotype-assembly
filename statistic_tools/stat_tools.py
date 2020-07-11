from matplotlib import pyplot
from file_handler import file_reader as fr
from matrix_scripts import matrix_tools as mt
import numpy as np
import time
import reads_handler as rh
import collections

def plot_ranked_hist(domination_stats: list):
    first_ranks = []
    second_ranks = []
    for stat_dict in domination_stats:
        stat = sorted(list(stat_dict.values()), reverse=True)
        if len(stat) <= 1:
            print("What??????????")
            continue
        count = sum(stat)
        first_ranks.append(stat[0] / count)
        second_ranks.append(stat[1] / count)
    bins = np.linspace(0, 1, 30)
    pyplot.figure()
    pyplot.hist(first_ranks, bins, alpha=0.5, label='rank_1s', color='red')
    pyplot.hist(second_ranks, bins, alpha=0.5, label='rank_2s', color='blue')
    pyplot.legend(loc='center')
    pyplot.show()
    time.sleep(2)


def analyse_result(version: str, part_number: int = -1):
    mat, reads, dicts = fr.load_results(version, part_number=part_number)
    rounded_mat = mt.complex_matrix_projection(mat)
    zero_reads = 0
    for read in reads:
        for e in read:
            if e == "X":
                zero_reads += 1
    print("number of completed elements: {} out of {} zeros".format(zero_reads - np.sum(rounded_mat == 0), zero_reads))
    return mat, reads, dicts


def find_pair_block_frequencies(all_exon_blocks: list, individuals_exon_blocks: list):
    freq_tables = []
    for exon_blocks, per_exon_individual_blocks in zip(all_exon_blocks, individuals_exon_blocks):
        freq_table = find_pair_block_frequencies_for_single_exon(exon_blocks, per_exon_individual_blocks)
        freq_tables.append(freq_table)
    return freq_tables


def find_pair_block_frequencies_for_single_exon(exon_blocks: list, individuals_exon_blocks: list):
    print("############################")
    m = len(exon_blocks)
    if m == 1:
        return np.array([len(individuals_exon_blocks)])

    double_freq_table = np.zeros((m, m))

    read_less_individuals = 0
    for ind_blocks in individuals_exon_blocks:
        ind_blocks_mapping = rh.map_blocks_to_source_blocks(source_blocks=exon_blocks, individual_blocks=ind_blocks)
        detected_exon_blocks = list(collections.Counter(ind_blocks_mapping).keys())
        print(detected_exon_blocks)
        if len(detected_exon_blocks) == 0:
            read_less_individuals += 1
        if len(detected_exon_blocks) == 1:
            index = detected_exon_blocks[0]
            double_freq_table[index, index] += 1
        else:
            for i1 in detected_exon_blocks:
                for i2 in detected_exon_blocks:
                    if i1 != i2:
                        double_freq_table[i1, i2] += 1
    double_freq_table = double_freq_table/(len(individuals_exon_blocks) - read_less_individuals)
    return double_freq_table


def get_frequency(l: list):
    d = {}
    for e in l:
        if list(d.keys()).__contains__(e):
            d[e] += 1
        else:
            d[e] = 1
    return d


def main_double_freqs():
    exons_count = 9
    individuals_count = 75
    all_individual_blocks = fr.load_all_individual_blocks("svt", individuals_count, 9)
    all_exon_blocks = fr.load_all_exon_blocks("exons-svt", max_parts=exons_count)
    per_exon_individual_blocks = [[ind_blocks[i] for ind_blocks in all_individual_blocks] for i in range(exons_count)]
    freq_tables = find_pair_block_frequencies(all_exon_blocks, per_exon_individual_blocks)
    c = 0
    for freq_table, exon_blocks in zip(freq_tables, all_exon_blocks):
        c += 1
        print("***************** exon {} ******************".format(c))
        print("exon blocks:")
        for block in exon_blocks:
            print(block)
        print("frequencies:")
        print("entry (i,j) indicates freq of b_i and b_j,\n (i,i) indicates frequency of b_i\n")
        freqs = np.floor(freq_table*100) / 100
        print(freqs)


if __name__ == '__main__':
    main_double_freqs()
