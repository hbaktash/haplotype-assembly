import numpy as np
from file_handler import file_reader
from Matrix_completion import HapNuc, HapSVT
from matrix_scripts import matrix_tools as mt
from matrix_scripts import reads_handler as rh
import os

SEPERATION_POSITIONS = []


def final_block_results(completed_matrix: np.ndarray, idx_pairs: list, version, part_number, cheat=False):
    print("finding blocks")
    haplo_blocks = mt.extract_distinctive_blocks(completed_matrix,
                                                 single_individual=False,
                                                 part_number=part_number,
                                                 cheat=cheat)

    haplo_blocks_as_string = rh.read_mat_to_arrays(haplo_blocks, None, use_dicts=False)
    haplo_blocks_as_string, _ = rh.remove_repeated_reads(haplo_blocks_as_string)
    haplo_blocks_as_string, _ = rh.remove_submissive_reads(haplo_blocks_as_string)
    print("final results:")
    print("variant indexes:")
    print([a[0] for a in idx_pairs])
    with open(os.path.join(file_reader.DATA_PATH,"results","freqs.txt"),'a+') as file:
        file.write("distinct blocks:\n")
        print("distinct blocks:")
        block_frequencies = rh.get_block_frequencies(haplo_blocks_as_string,
                                                     rh.read_mat_to_arrays(mt.complex_matrix_projection(completed_matrix),
                                                                           None,
                                                                           use_dicts=False))
        for block, freq in zip(haplo_blocks_as_string, block_frequencies):
            print(block)
            file.write(str(block))
            file.write("\n")
            print("freq: {}% out of {} reads in this part".format(int(1000*(freq/sum(block_frequencies)))/10,
                                                                  completed_matrix.shape[0]))
            file.write("freq: {}% out of {} reads in this part\n".format(int(1000*(freq/sum(block_frequencies)))/10,
                                                                  completed_matrix.shape[0]))
    print("saving final result to file")
    file_reader.save_as_file(haplo_blocks_as_string, "extracted-blocks", version=version, part_number=part_number)


if __name__ == '__main__':
    test_title = (input("enter test title\n"))
    separation_positions = [6644200, 6645000, 6645826, 6646025, 6646235, 6646439, 6646661, 6647241]
    exon_intervals = [(6643683, 6644027), (6644430, 6644635), (6645660, 6645759), (6645850, 6645956),
                      (6646083, 6646176), (6646267, 6646382), (6646475, 6646556), (6646750, 6647162),
                      (6647267, 6647537)]
    # separation_positions = [int(a) for a in input("enter separations indexes:\n").split(" ")]
    print("reading files")
    all_read_arrs, all_index_pars = file_reader.read_merged()
    print("separating reads")
    # independent_read_arrays, sub_index_pairs = rh.separate_reads(all_read_arrs, all_index_pars,
    #                                                              separation_positions)
    independent_read_arrays, sub_index_pairs = rh.separate_exons(all_read_arrs, all_index_pars,
                                                                 exon_intervals=exon_intervals)
    print("saving indexes")
    file_reader.save_indexes([a[0] for a in all_index_pars], False)
    file_reader.save_indexes([[a[0] for a in index_par_list] for index_par_list in sub_index_pairs], True)
    part = 0
    print("working on parts:")
    for read_arrays, index_pairs in zip(independent_read_arrays, sub_index_pairs):
        part += 1
        print("            PART {}           ".format(part))
        # print("     {} building matrix and dicts".format(part))
        # read_mat, dicts = rh.read_array_to_matrix(read_arrays, index_pairs, difference_magnitude=1j)
        # print("     matrix shape:", read_mat.shape)
        # print("     {} completing matrix".format(part))
        # completed_mat = HapSVT.complete_matrix(read_mat, delta=0.8, shrinkage_threshold=1.9, epsilon=0.003,
        #                                        verbose=False)
        # print("     {} saving to file..".format(part))
        # file_reader.save_result(completed_mat, read_arrays, dicts, test_title, part_number=part)
        # u, s, vh = mt.svd(completed_mat)
        # print("     {} svd:\nshape:".format(part))
        # print(np.array(s).shape)
        # print("")
        # _, s2, _ = mt.svd(completed_mat)
        # print("     {} completed s:\n".format(part), list(s2))

        completed_mat, _, _ = file_reader.load_results(test_title, part_number=part)
        print("         reads size {}".format(completed_mat.shape))
        # file_reader.save_parts_as_mat(version=test_title, part_number=part)
        final_block_results(completed_mat, index_pairs, test_title, part, cheat=True)
