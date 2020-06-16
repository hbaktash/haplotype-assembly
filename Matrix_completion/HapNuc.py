import numpy as np
import cvxpy as cp
from matrix_scripts import matrix_tools as mt


def matrix_completion(R, epsilon=0.1):
    selection_matrix = (R == 0)
    known_value_indices = tuple(zip(*[[x[0], x[1]] for x in list(zip(np.where(R != 0)[0], np.where(R != 0)[1]))]))
    known_values = (R[R != 0])
    print(R)
    X = cp.Variable(R.shape)
    objective_fn = cp.norm(X, p='nuc')
    constraints = [
        cp.norm(X[known_value_indices] - known_values) <= epsilon
    ]
    problem = cp.Problem(objective=cp.Minimize(objective_fn), constraints=constraints)
    print("CONS:", problem.status)
    problem.solve(gp=False, verbose=True)
    return X.value


def hap_nuc(R):
    H_hat = matrix_completion(R)

    h_p, h_m = mt.extract_distinctive_blocks(H_hat)
    return h_p, h_m


def test():
    r = np.array([[1, 1, 1, 0],
                  [-1, 0.8, -1, 0],
                  [0, 1, 1, -1],
                  [1, -1, 0, -1],
                  [0, -1, 1, -1]
                  ])
    r_hat = matrix_completion(r, epsilon=1)
    print(np.round(r_hat))
    print(type(r_hat))
    _, s, _ = mt.svd(r_hat)
    print(np.round(s))


if __name__ == '__main__':
    test()

