#! /usr/bin/python3

import pytest

import greedy_reversal_sorting

def init_test_from_file(path):

    with open(path) as f:
        permstr = f.readline().rstrip()

    permutation = []
    for x in permstr.split():
        permutation.append(int(x))

    results = greedy_reversal_sorting.greedy_sort(permutation)
    return results

def file_based_test(in_file, assert_file):
    results = init_test_from_file(in_file)
    result_idx = 0
    with open(assert_file) as f:
        assert_val = f.readline().rstrip()
        assert results[result_idx] == assert_val
        result_idx += 1


def test_sample_0():
    rt = "test/greedy_reversal_sorting/"
    file_based_test(rt + "input/sample0.txt", rt + "output/sample0.txt")


def test_sample_1():
    rt = "test/greedy_reversal_sorting/"
    file_based_test(rt + "input/sample1.txt", rt + "output/sample1.txt")

def test_test_1():
    rt = "test/greedy_reversal_sorting/"
    file_based_test(rt + "input/test_1.txt", rt + "output/test_1.txt")

def test_test_2():
    rt = "test/greedy_reversal_sorting/"
    file_based_test(rt + "input/test_2.txt", rt + "output/test_2.txt")

def test_test_3():
    rt = "test/greedy_reversal_sorting/"
    file_based_test(rt + "input/test_3.txt", rt + "output/test_3.txt")

def test_test_4():
    rt = "test/greedy_reversal_sorting/"
    file_based_test(rt + "input/test_4.txt", rt + "output/test_4.txt")

