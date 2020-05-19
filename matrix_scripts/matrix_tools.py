import math

import numpy as np
import sympy


def svd(A):
    return np.linalg.svd(A, full_matrices=False)


def singular_value_shrinkage(A, threshold):
    u, s, vh = svd(A)
    print("*********************************\n", list(s))
    print(s.shape)
    s = s - threshold
    s[s < 0] = 0
    shrinked = np.matmul(u, np.matmul(np.diag(s), vh))
    return shrinked


def selective_matrix(A, flags):
    return np.multiply(A, flags)


def frob_norm(A):
    return np.linalg.norm(A, ord='fro')


def rref(A):
    a = sympy.Matrix(A)
    rref_mat, pivots = a.rref()
    return sympy.matrix2numpy(rref_mat), pivots


def extract_haplotypes(H):
    rref_H, p = rref(H.transpose())
    binary_H = np.ones(H.shape)
    binary_H[H > 0] = 2
    binary_H[H < 0] = 0
    binary_H = binary_H - 1
    h_p = binary_H[p[0], :]
    h_m = binary_H[p[1], :]
    return h_p, h_m


def complex_matrix_projection(r_hat: np.ndarray):
    ans = np.zeros(r_hat.shape, dtype=np.complex)
    for i in range(r_hat.shape[0]):
        print(i)
        for j in range(r_hat.shape[1]):
            ans[i, j] = complex_value_projection(r_hat[i, j])
    return ans


def complex_value_projection(c: complex):
    if complex_norm(c) < 1 / 2:
        return 0
    elif complex_norm(c) > 2:
        print("TOO BIG")
        return 0
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


def complex_norm(c: complex):
    return math.sqrt(c.real ** 2 + c.imag ** 2)


if __name__ == '__main__':
    A = np.random.rand(8, 4)
    u, s, vh = svd(A)
    print(u.shape, s.shape, vh.shape)
