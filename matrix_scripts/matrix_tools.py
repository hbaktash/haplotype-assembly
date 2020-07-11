import math

import numpy as np
import sympy
from scipy.linalg import lu
from scipy import io
import file_reader as fr
import os
import matlab
import reads_handler as rh
from Matrix_completion import HapSVT
import file_reader as fr


def svd(A):
    return np.linalg.svd(A, full_matrices=False)


def singular_value_shrinkage(A, threshold):
    u, s, vh = svd(A)
    # print("*********************************\n", list(s))
    # print(s.shape)
    s = s - threshold
    s[s < 0] = 0
    shrinked = np.matmul(u, np.matmul(np.diag(s), vh))
    return shrinked


def selective_matrix(A, flags):
    return np.multiply(A, flags)


def frob_norm(A):
    return np.linalg.norm(A, ord='fro')


def find_pivots(upper_tri: np.ndarray, epsilon=0.00001):
    pivots = []
    for i in range(upper_tri.shape[0]):
        for j in range(upper_tri.shape[1]):
            if complex_norm(upper_tri[i, j]) > epsilon:
                pivots.append(j)
                break
    return pivots


def rref(matrix, with_lu=False, cheat: bool = False, part_number=-1):
    if cheat:  # :D
        p = io.loadmat(os.path.join(fr.DATA_PATH,
                                    "results",
                                    "exons-large-svt",
                                    "exons-pivots",
                                    "p{}.mat".format(part_number)))['p{}'.format(part_number)]
        return None, p[0]
    if not with_lu:
        a = sympy.Matrix(matrix)
        rref_mat, pivots = a.rref(normalize_last=True)
        return sympy.matrix2numpy(rref_mat), pivots
    else:
        pl, u = lu(matrix, permute_l=True)
        pivots = find_pivots(upper_tri=u)
        return u, pivots


def get_block_frequencies(blocks: list, rounded_H:np.ndarray):
    freqs = []
    for block in blocks:
        freq = 0
        for numeric_read in list(rounded_H):
            if np.allclose(np.array(block), numeric_read):
                freq += 1
        freqs.append(freq)
    return freqs


def extract_distinctive_blocks(H: np.ndarray, single_individual: bool = True, part_number=-1, cheat: bool = False):
    if single_individual:
        rref_H, pivots = rref(H.transpose())
        binary_H = np.ones(H.shape)
        binary_H[H > 0] = 2
        binary_H[H < 0] = 0
        binary_H = binary_H - 1
        h_p = binary_H[pivots[0], :]
        h_m = binary_H[pivots[1], :]
        return np.array(h_p, h_m)
    else:
        # print("    -estimating rank")
        # estimated_r = estimated_rank(H)
        print("    -rref")
        _, pivots = rref(H.transpose(), with_lu=True, part_number=part_number, cheat=cheat)
        print("    -projecting")
        # print(H.shape)
        # print(list(H))
        rounded_H = complex_matrix_projection(H)
        print("    -pivots: ", len(pivots), list(pivots))
        haplo_blocks = []
        print("    -saving blocks")
        for i in range(len(pivots)):
            haplo_blocks.append(list(rounded_H[pivots[i], :]))
        # block_frequencies = get_block_frequencies(haplo_blocks, rounded_H)
        return np.array(haplo_blocks)


def estimated_rank(H: np.ndarray, singular_threshold: float = 0.02):
    _, singular_values, _ = svd(H)
    valid_count = np.sum(singular_values >= singular_threshold)
    return int(valid_count)


def complex_norm(c: complex):
    return math.sqrt(c.real ** 2 + c.imag ** 2)


def complex_value_projection(c: complex):
    if complex_norm(c) < 0.01:
        return 0
    if complex_norm(c) > 5:
        print("TOO BIG!!!")
    else:
        if c.real >= c.imag:
            if c.imag >= -c.real:
                return 1
            else:
                return -1j
        else:
            if c.imag >= -c.real:
                return 1j
            else:
                return -1


def complex_matrix_projection(r_hat: np.ndarray):
    # complex_projector_vectorized = np.vectorize(complex_value_projection)
    new_mat = np.zeros(r_hat.shape, dtype=np.complex)
    for i in range(new_mat.shape[0]):
        for j in range(new_mat.shape[1]):
            # print(i, j)
            # print(r_hat[i, j])
            new_mat[i, j] = complex_value_projection(r_hat[i, j])
    return new_mat
    # return complex_projector_vectorized(r_hat)


def single_individual_blocks(reads: list, indexes: list):
    pass


def complete_read_arrays(test_title, read_arrays, part: int):
    if len(read_arrays) == 0:
        completed_mat = [0]
        dicts = {"X": 0}
    else:
        print("     {} building matrix and dicts".format(part))
        read_mat, dicts = rh.read_array_to_matrix(read_arrays)
        print("     matrix shape:", read_mat.shape)
        print("     non-zero entries:", np.sum(read_mat != 0))
        print("     {} completing matrix".format(part))
        completed_mat = HapSVT.complete_matrix(read_mat, delta=0.99, shrinkage_threshold=0.8, epsilon=0.0001,
                                               verbose=False,
                                               min_iters=600)
        u, s, vh = svd(completed_mat)
        print("     {} completed s:\n".format(part), list(s))
    print("     {} saving to file..".format(part))
    fr.save_result(completed_mat, read_arrays, dicts, test_title, part_number=part)
    return completed_mat


# if __name__ == '__main__':
#     A = np.random.rand(8, 4)
#     u, s, vh = svd(A)
#     print(u.shape, s.shape, vh.shape)
