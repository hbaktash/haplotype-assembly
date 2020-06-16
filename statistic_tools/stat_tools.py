from matplotlib import pyplot
from file_handler import file_reader as fr
from matrix_scripts import matrix_tools as mt
import numpy as np
import time


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


def get_frequency(l: list):
    d = {}
    for e in l:
        if list(d.keys()).__contains__(e):
            d[e] += 1
        else:
            d[e] = 1
    return d
