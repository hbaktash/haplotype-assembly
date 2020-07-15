from matrix_scripts import matrix_tools as mt
import numpy as np


def complete_matrix(R: np.ndarray,
                    delta: float = 0.8,
                    shrinkage_threshold: float = 0.5,
                    epsilon: float = 0.02,
                    min_iters: int = 500,
                    verbose: bool = True):
    selection_matrix = np.ones(R.shape)
    selection_matrix[R == 0] = 0
    tmp_Y = R
    _, s, _ = mt.svd(R)
    print("initial s:\n", list(s))
    norm_R = mt.frob_norm(R)
    print("norm of R:", norm_R)
    k = 0
    i = 1
    tmp_X = np.zeros(R.shape)
    pre_X = np.zeros(R.shape)
    print("norm of R:\n", mt.frob_norm(R))
    number_of_iters = 0
    while True:
        tmp_X = mt.singular_value_shrinkage(tmp_Y, threshold=shrinkage_threshold)

        if verbose:
            diff = mt.frob_norm(mt.selective_matrix(tmp_X - R, selection_matrix))
            print("known entries difference: ", diff)
        diff = mt.frob_norm(pre_X - tmp_X)
        if (diff <= epsilon*norm_R) and (number_of_iters >= min_iters):  # diff <= epsilon*norm_R and
            break
        if verbose:
            print("X difference: ", diff)
        tmp_Y = tmp_X - delta*mt.selective_matrix(tmp_X - R, selection_matrix)
        pre_X = tmp_X
        number_of_iters += 1
    print("     ***number of iterations:", number_of_iters)
    print("     ***difference with init:", mt.frob_norm(mt.selective_matrix(R - tmp_X, selection_matrix)))
    return tmp_X


def hapt_SVT(R: np. ndarray):
    H_hat = complete_matrix(R)
    h_p, h_m = mt.extract_distinctive_blocks(H_hat)
    return h_p, h_m


# if __name__ == '__main__':
    # r = np.array([
    #     [1, 0, 0, -1],
    #     [0, 0, -1, 1],
    #     [1, 1, 0, 0],
    #     [1, 1, 0, -1],
    #     [0, 1, 0, -1],
    #     [-1, 1, 0, 0],
    #     [-1, 0, 0, 1],
    #     [0, 0, 1, -1],
    #      ])
    # p, m = hapt_SVT(r)
    # print(p)
    # print("**************************\n", m)
