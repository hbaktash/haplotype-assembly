from matrix_scripts import matrix_tools as mt
import numpy as np


def complete_matrix(R: np.ndarray, delta: float = 0.1, shrinkage_threshold: float = 0.5, epsilon: float = 0.02):
    selection_matrix = np.ones(R.shape)
    selection_matrix[R == 0] = 0
    tmp_Y = R
    k = 0
    i = 1
    tmp_X = np.zeros(R.shape)
    print("norm of R:\n", mt.frob_norm(R))
    while True:
        tmp_X = mt.singular_value_shrinkage(tmp_Y, threshold=shrinkage_threshold)
        tmp_Y = tmp_Y + delta*mt.selective_matrix(R - tmp_X, selection_matrix)

        if mt.frob_norm(mt.selective_matrix(tmp_X - R, selection_matrix)) <= epsilon*mt.frob_norm(R):
            break
        print(mt.frob_norm(mt.selective_matrix(tmp_X - R, selection_matrix)) - epsilon*mt.frob_norm(R))
    return tmp_X


def hapt_SVT(R: np. ndarray):
    H_hat = complete_matrix(R)
    h_p, h_m = mt.extract_haplotypes(H_hat)
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
