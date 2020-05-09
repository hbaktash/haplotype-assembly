import numpy as np
import sympy

def svd(A):
    return np.linalg.svd(A, full_matrices=False)


def singular_value_shrinkage(A, threshold):
    u, s, vh = svd(A)
    s = s - threshold
    s[s < 0] = 0
    return np.matmul(u, np.matmul(np.diag(s), vh))


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

if __name__ == '__main__':
    A = np.random.rand(8, 4)
    u, s, vh = svd(A)
    print(u.shape, s.shape, vh.shape)
