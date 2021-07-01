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

def globalalignment_verify(results, verifyfile):

    with open(verifyfile) as f:
        score = int(f.readline())
        val_1 = f.readline().rstrip()
        val_2 = f.readline().rstrip()

    for i in range(3):
        print(results[i])

    assert results[0] == score
    assert results[1] == val_1
    assert results[2] == val_2

def localalignment_test(testfile):
    with open(testfile) as f:
        nucleotide_h = f.readline().rstrip()
        nucleotide_w = f.readline().rstrip()

    scoring_strs = []
    with open("./PAM250_scoring.txt") as f:
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

    return localalignment.lcsbacktrack(nucleotide_h, nucleotide_w, scoring)

def localalignment_verify(results, verifyfile):

    with open(verifyfile) as f:
        score = int(f.readline())
        val_1 = f.readline().rstrip()
        val_2 = f.readline().rstrip()

    for i in range(3):
        print(results[i])
        print("aaa")

    assert results[0] == score
    assert results[1] == val_1
    assert results[2] == val_2


def fittingalignment_test(testfile):
    with open(testfile) as f:
        nucleotide_h = f.readline().rstrip()
        nucleotide_w = f.readline().rstrip()

    scoring_strs = []
    with open("./PAM250_scoring.txt") as f:
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

    return fittingalignment.lcsbacktrack(nucleotide_h, nucleotide_w, scoring)

def fittingalignment_verify(results, verifyfile):

    with open(verifyfile) as f:
        score = int(f.readline())
        val_1 = f.readline().rstrip()
        val_2 = f.readline().rstrip()

    for i in range(3):
        print(results[i])
        print("aaa")

    assert results[0] == score
    assert results[1] == val_1
    assert results[2] == val_2


def overlapalignment_test(testfile):
    with open(testfile) as f:
        nucleotide_h = f.readline().rstrip()
        nucleotide_w = f.readline().rstrip()

    scoring_strs = []
    with open("./PAM250_scoring.txt") as f:
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

    return overlapalignment.lcsbacktrack(nucleotide_h, nucleotide_w, scoring)

def overlapalignment_verify(results, verifyfile):

    with open(verifyfile) as f:
        score = int(f.readline())
        val_1 = f.readline().rstrip()
        val_2 = f.readline().rstrip()

    print("aaa")
    print(verifyfile)
    for i in range(3):
        print(results[i])

    print(score)
    print(val_1)
    print(val_2)

    assert results[0] == score
    assert results[1] == val_1
    assert results[2] == val_2

def globalalignment_test_1():
    results = globalalignment_test("test/globalalignment.txt")
    globalalignment_verify(results, "test/globalalignment_expected_result.txt")

def globalalignment_test_2():
    results = globalalignment_test("test/global_alignment_test_dataset_2.txt")
    globalalignment_verify(results, "test/global_alignment_test_dataset_2_expected_result.txt")

def globalalignment_test_3():
    results = globalalignment_test("test/dataset_247_3.txt")
    globalalignment_verify(results, "test/dataset_247_3_expected_result.txt")

def localalignment_test_1():
    results = localalignment_test("test/localalignment.txt")
    localalignment_verify(results, "test/localalignment_expected_result.txt")

def localalignment_test_2():
    results = localalignment_test("test/local_alignment_test_2.txt")
    localalignment_verify(results, "test/local_alignment_test_2_expected_result.txt")

def localalignment_test_3():
    results = localalignment_test("test/local_alignment_test_3.txt")
    localalignment_verify(results, "test/local_alignment_test_3_expected_result.txt")

def localalignment_test_4():
    results = localalignment_test("test/local_alignment_test_4.txt")
    localalignment_verify(results, "test/local_alignment_test_4_expected_result.txt")

def localalignment_test_5():
    results = localalignment_test("test/local_alignment_test_5.txt")
    localalignment_verify(results, "test/local_alignment_test_5_expected_result.txt")

def fittingalignment_test_1():
    results = fittingalignment_test("test/08_FittingAlignment/inputs/test1.txt")
    fittingalignment_verify(results, "test/08_FittingAlignment/outputs/test1.txt")

def fittingalignment_test_2():
    results = fittingalignment_test("test/08_FittingAlignment/inputs/test2.txt")
    fittingalignment_verify(results, "test/08_FittingAlignment/outputs/test2.txt")

def fittingalignment_test_3():
    results = fittingalignment_test("test/08_FittingAlignment/inputs/test3.txt")
    fittingalignment_verify(results, "test/08_FittingAlignment/outputs/test3.txt")

def fittingalignment_test_4():
    results = fittingalignment_test("test/08_FittingAlignment/inputs/test4.txt")
    fittingalignment_verify(results, "test/08_FittingAlignment/outputs/test4.txt")

def fittingalignment_test_5():
    results = fittingalignment_test("test/08_FittingAlignment/inputs/test5.txt")
    fittingalignment_verify(results, "test/08_FittingAlignment/outputs/tes51.txt")

def overlapalignment_test_1():
    results = overlapalignment_test("test/overlapalignment.txt")
    overlapalignment_verify(results, "test/overlapalignment_expected_result.txt")

def overlapalignment_test_2():
    results = overlapalignment_test("test/09_OverlapAlignment/inputs/test1.txt")
    overlapalignment_verify(results, "test/09_OverlapAlignment/outputs/test1.txt")

def overlapalignment_test_3():
    results = overlapalignment_test("test/09_OverlapAlignment/inputs/test2.txt")
    overlapalignment_verify(results, "test/09_OverlapAlignment/outputs/test2.txt")

def overlapalignment_test_4():
    results = overlapalignment_test("test/09_OverlapAlignment/inputs/test3.txt")
    overlapalignment_verify(results, "test/09_OverlapAlignment/outputs/test3.txt")

def overlapalignment_test_5():
    results = overlapalignment_test("test/09_OverlapAlignment/inputs/test4.txt")
    overlapalignment_verify(results, "test/09_OverlapAlignment/outputs/test4.txt")

def overlapalignment_test_6():
    results = overlapalignment_test("test/09_OverlapAlignment/inputs/test5.txt")
    overlapalignment_verify(results, "test/09_OverlapAlignment/outputs/test5.txt")

def overlapalignment_test_7():
    results = overlapalignment_test("test/09_OverlapAlignment/inputs/test6.txt")
    overlapalignment_verify(results, "test/09_OverlapAlignment/outputs/test6.txt")

def overlapalignment_test_8():
    results = overlapalignment_test("test/dataset_248_7.txt")
    overlapalignment_verify(results, "test/dataset_247_10_expected_result.txt")


def run_all_tests():
    #globalalignment_test_1()
    #globalalignment_test_2()
    #globalalignment_test_3()
    
    localalignment_test_1()
    localalignment_test_2()
    localalignment_test_3()
    localalignment_test_4()
    localalignment_test_5()
    #localalignment_test_6()

    fittingalignment_test_1()
    fittingalignment_test_2()
    fittingalignment_test_3()
    fittingalignment_test_4()

    overlapalignment_test_1()
    overlapalignment_test_2()
    overlapalignment_test_3()
    overlapalignment_test_4()
    overlapalignment_test_5()
    #overlapalignment_test_6()
    #overlapalignment_test_7()
    overlapalignment_test_8()

    return

if __name__ == '__main__':
    start = time.process_time()
    run_all_tests()
    end = time.process_time()
    print("Time: {0}".format(end-start))