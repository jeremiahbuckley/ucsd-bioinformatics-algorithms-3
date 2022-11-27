#! /usr/bin/python3

import pytest

import two_break_distance

def init_test_from_file(path):

    cycles_set = two_break_distance.organize_inputs(path)
    genomes, synteny_count = two_break_distance.convert_cycles_strings_to_int_lists(cycles_set)

    results = two_break_distance.calc_two_break_distance(genomes, synteny_count)
    return results

def file_based_test(in_file, assert_file):
    results = init_test_from_file(in_file)
    result_idx = 0
    with open(assert_file) as f:
        assert results == int(f.readline().rstrip())


def test_sample_0():
    rt = "test/two_break_distance/"
    file_based_test(rt + "input/sample_0.txt", rt + "output/sample_0.txt")

def test_sample_1():
    rt = "test/two_break_distance/"
    file_based_test(rt + "input/sample_1.txt", rt + "output/sample_1.txt")

def test_sample_2():
    rt = "test/two_break_distance/"
    file_based_test(rt + "input/sample_2.txt", rt + "output/sample_2.txt")

