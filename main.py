import numpy as np
from file_handler import file_reader
from Matrix_completion import HapNuc, HapSVT
from matrix_scripts import matrix_tools as mt

if __name__ == '__main__':
    test_title = (input("enter test title\n"))
    print("reading files")
    read_arrs, index_pars = file_reader.read_merged()
    print("building matrix and dicts")
    read_mat, dicts = file_reader.read_array_to_matrix(read_arrs, index_pars, difference_magnitude=1j)
    print("completing matrix")
    completed_mat = HapSVT.complete_matrix(read_mat, delta=0.9, shrinkage_threshold=150, epsilon=0.001)
    print("saving to file..")
    file_reader.save_result(completed_mat, read_arrs, dicts)

    u, s, vh = mt.svd(completed_mat)
    print("svd:\nshape:")
    print(np.array(s).shape)
    print("")
    _, s2, _ = mt.svd(completed_mat)
    print("completed s:\n", s2)
