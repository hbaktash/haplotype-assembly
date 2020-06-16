from matrix_scripts import matrix_tools as mt
import numpy as np


def complete_matrix(R: np.ndarray, delta: float = 0.1, shrinkage_threshold: float = 0.5, epsilon: float = 0.02, verbose: bool = True):
    selection_matrix = np.ones(R.shape)
    selection_matrix[R == 0] = 0
    tmp_Y = R
    _, s, _ = mt.svd(R)
    print("initial s:\n", list(s))
    norm_R = mt.frob_norm(R)
    k = 0
    i = 1
    tmp_X = np.zeros(R.shape)
    pre_X = np.zeros(R.shape)
    print("norm of R:\n", mt.frob_norm(R))
    while True:
        tmp_X = mt.singular_value_shrinkage(tmp_Y, threshold=shrinkage_threshold*delta)

        diff = mt.frob_norm(mt.selective_matrix(tmp_X - R, selection_matrix))
        if verbose:
            print("known entries difference: ", diff)
        diff = mt.frob_norm(pre_X - tmp_X)
        if diff <= epsilon*norm_R:
            break
        if verbose:
            print("X difference: ", diff)
        tmp_Y = tmp_X - delta*mt.selective_matrix(tmp_X - R, selection_matrix)
        pre_X = tmp_X
    return tmp_X


def hapt_SVT(R: np. ndarray):
    H_hat = complete_matrix(R)
    h_p, h_m = mt.extract_distinctive_blocks(H_hat)
    return h_p, h_m


if __name__ == '__main__':
    r = np.array([
        [1, 0, 0, -1],
        [0, 0, -1, 1],
        [1, 1, 0, 0],
        [1, 1, 0, -1],
        [0, 1, 0, -1],
        [-1, 1, 0, 0],
        [-1, 0, 0, 1],
        [0, 0, 1, -1],
         ])
    p, m = hapt_SVT(r)
    print(p)
    print("**************************\n", m)
