from matplotlib import pyplot
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
