#! /usr/bin/python3

import pytest

import two_break_sorting

def init_test_from_file(path):

    cycles_set = two_break_sorting.organize_inputs(path)
    genomes, synteny_count = two_break_sorting.convert_cycles_strings_to_int_lists(cycles_set)

    process, tb_distances = two_break_sorting.transform_second_genome_to_first(genomes,synteny_count)
    results = two_break_sorting.transform_graph_lists_to_synteny_strings(process)
    return results, tb_distances

# note:
# an error in this test is that it takes tb_distance from the testable code, rather than re-calculating it itself. too trusting

def file_based_test(in_file, assert_file):
    results, tb_distances = init_test_from_file(in_file)

    assert_values = []
    with open(assert_file) as f:
        for line in f.readlines():
            assert_values.append(line.rstrip())

    print(results)
    print(assert_values)
    assert len(results) == len(assert_values)
    assert results[0] == assert_values[0]
    assert results[len(results) -1] == assert_values[len(assert_values) - 1]
    if len(results) > 1:
        prev = tb_distances[0]
        for next in tb_distances[1:]:
            assert prev > next
            prev = next


def test_sample_0():
    rt = "test/two_break_sorting/"
    file_based_test(rt + "input/sample_0.txt", rt + "output/sample_0.txt")

def test_sample_1():
    rt = "test/two_break_sorting/"
    file_based_test(rt + "input/sample_1.txt", rt + "output/sample_1.txt")

def test_sample_2():
    rt = "test/two_break_sorting/"
    file_based_test(rt + "input/sample_2.txt", rt + "output/sample_2.txt")

def test_sample_3():
    rt = "test/two_break_sorting/"
    file_based_test(rt + "input/sample_3.txt", rt + "output/sample_3.txt")

def test_sample_4():
    rt = "test/two_break_sorting/"
    file_based_test(rt + "input/sample_4.txt", rt + "output/sample_4.txt")

def test_sample_5():
    rt = "test/two_break_sorting/"
    file_based_test(rt + "input/sample_5.txt", rt + "output/sample_5.txt")
