#! /usr/bin/python3

import pytest

import number_of_breakpoints

def init_test_from_file(path):

    with open(path) as f:
        permstr = f.readline().rstrip()

    with open(path) as f:
        permstr = f.readline().rstrip()

    permutation = [0]
    for x in permstr.split():
        permutation.append(int(x))
    permutation.append(len(permutation))

    results = number_of_breakpoints.num_of_breakpoints(permutation)
    return results

def file_based_test(in_file, assert_file):
    results = init_test_from_file(in_file)
    result_idx = 0
    with open(assert_file) as f:
        assert results == int(f.readline().rstrip())


def test_sample_0():
    rt = "test/number_of_breakpoints/"
    file_based_test(rt + "input/sample0.txt", rt + "output/sample0.txt")

