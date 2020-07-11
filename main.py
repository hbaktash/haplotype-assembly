import numpy as np
from file_handler import file_reader as fr
from Matrix_completion import HapNuc, HapSVT
from matrix_scripts import matrix_tools as mt
from matrix_scripts import reads_handler as rh
import os
from statistic_tools import stat_tools

SEPERATION_POSITIONS = []


def remove_all_individual_insane_blocks(version: str, number_of_individuals: int = 1492, parts: int = 9):
    individual_blocks = fr.load_all_individual_blocks(version, number_of_individuals)
    individual_read_arrays = fr.load_all_individual_read_arrays(version, number_of_individuals)
    ind_count = 0
    for individual_reads_by_exon, individual_blocks_by_exon in zip(individual_read_arrays, individual_blocks):
        print("----------------------------------------------")
        ind_count += 1
        part_count = 0
        for reads, blocks in zip(individual_reads_by_exon, individual_blocks_by_exon):
            part_count += 1
            if (reads == None) or (len(reads) < 1):
                continue
            sane_blocks = rh.remove_insane_blocks(blocks, reads)
            if len(sane_blocks) > 2:
                print("still insane ind {} part {} : {} -> {}".format(ind_count,
                                                                      part_count,
                                                                      len(blocks),
                                                                      len(sane_blocks)))
            fr.save_as_file(sane_blocks, "extracted-blocks", version, part_count, ind_count)


def final_block_results(completed_numeric_matrix: np.ndarray, idx_pairs: list, version, part_number, cheat=False,
                        use_rref: bool = False):
    print("finding blocks")
    low_value = 0.5
    print("low value enteries:", np.sum(np.abs(completed_numeric_matrix) < low_value))
    if not use_rref:
        haplo_blocks_as_string = []
        projected = mt.complex_matrix_projection(completed_numeric_matrix)
        projected_as_strings = rh.read_mat_to_arrays(projected, None, use_dicts=False)
        for semi_block in projected_as_strings:
            if haplo_blocks_as_string.__contains__(semi_block):
                continue
            else:
                haplo_blocks_as_string.append(semi_block)
        print("found {} blocks".format(len(haplo_blocks_as_string)))

    else:
        haplo_blocks = mt.extract_distinctive_blocks(completed_numeric_matrix,
                                                     single_individual=False,
                                                     part_number=part_number,
                                                     cheat=cheat)
        haplo_blocks_as_string = rh.read_mat_to_arrays(haplo_blocks, None, use_dicts=False)
    haplo_blocks_as_string, _ = rh.remove_repeated_reads(haplo_blocks_as_string)
    haplo_blocks_as_string, _ = rh.remove_submissive_reads(haplo_blocks_as_string)
    print("final results:")
    print("variant indexes:")
    print([a[0] for a in idx_pairs])
    with open(os.path.join(fr.DATA_PATH, "results", version, "freqs.txt"), 'a+') as file:
        # file.write("distinct blocks:\n")
        print("distinct blocks:")
        block_frequencies = rh.get_block_frequencies(haplo_blocks_as_string,
                                                     rh.read_mat_to_arrays(
                                                         mt.complex_matrix_projection(completed_numeric_matrix),
                                                         None,
                                                         use_dicts=False))
        print("sums:", sum(block_frequencies))
        print("matrix shape:", completed_numeric_matrix.shape)
        for block, freq in zip(haplo_blocks_as_string, block_frequencies):
            print(block)
            print("freq: {}% out of {} reads in this part".format(freq,  # int(1000*(freq/sum(block_frequencies)))/10
                                                                  completed_numeric_matrix.shape[0]))
    print("saving final result to file")
    fr.save_as_file(haplo_blocks_as_string, "extracted-blocks", version=version, part_number=part_number)


def main(for_individuals: bool):
    if for_individuals:
        print("This is for all individuals, stop otherwise.")
    test_title = (input("enter test title\n"))
    # separation_positions = [6644200, 6645000, 6645826, 6646025, 6646235, 6646439, 6646661, 6647241]
    exon_intervals = [(6643683, 6644027), (6644430, 6644635), (6645660, 6645759), (6645850, 6645956),
                      (6646083, 6646176), (6646267, 6646382), (6646475, 6646556), (6646750, 6647162),
                      (6647267, 6647537)]
    print("reading files")
    all_index_pairs = fr.read_merged_indexes()
    all_indexes = [a[0] for a in all_index_pairs]

    if for_individuals:
        all_read_arrs = fr.load_individuals_reads(all_indexes)
    else:
        all_read_arrs = [fr.read_merged(all_index_pairs)]

    indi_counter = 0
    for read_arrs in all_read_arrs:
        indi_counter += 1
        print("#######################################################################################################")
        print("#######################################################################################################")
        print("#######################################################################################################")
        print("individual number {}/{}".format(indi_counter, len(all_read_arrs)))

        print("separating reads")
        independent_read_arrays, sub_index_pairs = rh.separate_exons(read_arrs, all_index_pairs,
                                                                     exon_intervals=exon_intervals)
        print("saving indexes")
        fr.save_indexes([a[0] for a in all_index_pairs], False)
        fr.save_indexes([[a[0] for a in index_par_list] for index_par_list in sub_index_pairs], True)
        part = 0
        print("working on parts:")
        for read_arrays, index_pairs in zip(independent_read_arrays, sub_index_pairs):
            print("***********************************************************************************")
            part += 1
            print("            PART {}           ".format(part))
            title = test_title
            if for_individuals:
                title = "individual-{}-{}".format(indi_counter, test_title)
            completed_numeric_mat = mt.complete_read_arrays(title, read_arrays, part)
            # completed_numeric_mat, _, _ = fr.load_results(title, part_number=part)
            # print("         reads size {}".format(len(completed_numeric_mat)))
            if len(completed_numeric_mat) == 1:
                continue
            # # # fr.save_part_as_mat(version=test_title, part_number=part)
            final_block_results(completed_numeric_mat, index_pairs, title, part, cheat=True, use_rref=False)


if __name__ == '__main__':
    # main(for_individuals=True)
    stat_tools.main_double_freqs("large-svt")
    # remove_all_individual_insane_blocks("large-svt")
