#! /usr/bin/python3

import pytest

import mlcs

def init_test_from_file(path):

    with open(path) as f:
        nuc_1 = f.readline().rstrip()
        nuc_2 = f.readline().rstrip()
        nuc_3 = f.readline().rstrip()

    match_reward = 1
    mismatch_penalty = 0
    indel_penalty = 0
    scoring = mlcs.Scoring(match_reward, mismatch_penalty, indel_penalty)

    results = mlcs.find_common_subseq([nuc_1, nuc_2, nuc_3],  scoring)
    return results

def file_based_test(in_file, assert_file):
    results = init_test_from_file(in_file)
    with open(assert_file) as f:
        assert_val0 = int(f.readline().rstrip())
        assert_val1 = f.readline().rstrip()
        assert_val2 = f.readline().rstrip()
        assert_val3 = f.readline().rstrip()
    assert results[0] == assert_val0
    found_winner = False
    for set in results[1]:
        if (set[0] == assert_val1 and set[1] == assert_val2 and set[2] == assert_val3) or \
           (set[0] == assert_val1 and set[2] == assert_val2 and set[1] == assert_val3) or \
           (set[1] == assert_val1 and set[0] == assert_val2 and set[2] == assert_val3) or \
           (set[1] == assert_val1 and set[2] == assert_val2 and set[0] == assert_val3) or \
           (set[2] == assert_val1 and set[0] == assert_val2 and set[1] == assert_val3) or \
           (set[2] == assert_val1 and set[1] == assert_val2 and set[0] == assert_val3):
            found_winner = True
    assert found_winner

def test_sample_0():
    rt = "test/13_MultipleSequenceAlignment/"
    file_based_test(rt + "input/sample.txt", rt + "output/sample.txt")

# def test_sample_2():
#     rt = "test/13_MultipleSequenceAlignment/"
#     file_based_test(rt + "input/sample2.txt", rt + "output/sample2.txt")

def test_input_1():
    rt = "test/13_MultipleSequenceAlignment/"
    file_based_test(rt + "input/input_1.txt", rt + "output/output_1.txt")

def test_input_2():
    rt = "test/13_MultipleSequenceAlignment/"
    file_based_test(rt + "input/input_2.txt", rt + "output/output_2.txt")

def test_input_3():
    rt = "test/13_MultipleSequenceAlignment/"
    file_based_test(rt + "input/input_3.txt", rt + "output/output_3.txt")

def test_input_4():
    rt = "test/13_MultipleSequenceAlignment/"
    file_based_test(rt + "input/input_4.txt", rt + "output/output_4.txt")

def test_input_5():
    rt = "test/13_MultipleSequenceAlignment/"
    file_based_test(rt + "input/input_5.txt", rt + "output/output_5.txt")

def test_input_6():
    rt = "test/13_MultipleSequenceAlignment/"
    file_based_test(rt + "input/input_6.txt", rt + "output/output_6.txt")

# def test_input_7():
#     rt = "test/13_MultipleSequenceAlignment/"
#     file_based_test(rt + "/input/input_7.txt", rt + "output/output_7.txt")




