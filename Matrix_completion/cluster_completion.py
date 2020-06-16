import numpy as np
from statistic_tools import stat_tools as st


def correlation_matrix(R: np.ndarray):
    return np.matmul(R, R.transpose().conjugate())


def fix_diagonals(corr_matrix: np.ndarray):
    corr_matrix[np.diag_indices(corr_matrix.shape[0], 2)] = 1


def min_cut_clusters(corr_matrix: np.ndarray):
    pass


def get_clusters(R: np.ndarray):
    correlation_mat = correlation_matrix(R)
    fix_diagonals(correlation_mat)
    clusters = min_cut_clusters(correlation_mat)
    return clusters


def get_block_from_clustered_reads(R_c: np.ndarray):
    freq_dicts = []
    for j in range(R_c.shape[1]):
        freq_dict = st.get_frequency(list(R_c[:, j]))
        del freq_dict[0]
        freq_dicts.append(freq_dict)
    # TODO
    pass


def complete_with_clusters(R: np.ndarray, clusters: list):
    # TODO
    pass
