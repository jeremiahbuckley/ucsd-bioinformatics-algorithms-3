#! /usr/bin/python3

import sys
import time

import globalalignment
import localalignment
import editdistance
import fittingalignment
import overlapalignment

def globalalignment_test(testfile):
    with open(testfile) as f:
        nucleotide_h = f.readline().rstrip()
        nucleotide_w = f.readline().rstrip()

    scoring_strs = []
    with open("./BLOSUM62.txt") as f:
        scoring_strs = f.readlines()

    scoring = {}
    virt_keys = 'ACDEFGHIKLMNPQRSTVWY'
    for i in range(1, len(scoring_strs)):
        row = scoring_strs[i].split()
        key = row[0]
        row_dict = {}
        for j in range(1, len(row)):
            virt_key = virt_keys[j-1]
            row_dict[virt_key] = int(row[j])
        scoring[key] = row_dict
    print(scoring)

    return globalalignment.lcsbacktrack(nucleotide_h, nucleotide_w, scoring)
    print(f_score)
    print(out_h)
    print(out_w)

def globalalignment_test_1():
    results = globalalignment_test("test/globalalignment.txt")

    with open("test/globalalignment_expected_result.txt") as f:
        score = int(f.readline())
        val_1 = f.readline().rstrip()
        val_2 = f.readline().rstrip()

    for i in range(3):
        print(results[i])

    assert results[0] == score
    assert results[1] == val_1
    assert results[2] == val_2

def globalalignment_test_2():
    results = globalalignment_test("test/global_alignment_test_dataset_2.txt")

    with open("test/global_alignment_test_dataset_2_expected_result.txt") as f:
        score = int(f.readline())
        val_1 = f.readline().rstrip()
        val_2 = f.readline().rstrip()

    for i in range(3):
        print(results[i])

    assert results[0] == score
    assert results[1] == val_1
    assert results[2] == val_2

def globalalignment_test_3():
    results = globalalignment_test("test/dataset_247_3.txt")

    with open("test/dataset_247_3_expected_result.txt") as f:
        score = int(f.readline())
        val_1 = f.readline().rstrip()
        val_2 = f.readline().rstrip()

    for i in range(3):
        print(results[i])

    assert results[0] == score
    assert results[1] == val_1
    assert results[2] == val_2

def run_all_tests():
    globalalignment_test_1()
    globalalignment_test_2()
    globalalignment_test_3()
    return

if __name__ == '__main__':
    start = time.process_time()
    run_all_tests()
    end = time.process_time()
    print("Time: {0}".format(end-start))