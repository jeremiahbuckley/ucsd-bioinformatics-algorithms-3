#! /usr/bin/python3

import sys
import time

import any_alignment

def print_scoring(scoring, keys):
    print("      " + "    ".join([k for k in keys]))
    for i in range(len(keys)):
        val_str = ""
        for j in range(len(keys)):
            val_str += " " + str(scoring[keys[i]][keys[j]]).rjust(4)
        print("{0} {1}".format(keys[i], val_str))

def load_scoring(scoring_file_loc):
    scoring_strs = []
    with open(scoring_file_loc) as f:
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

    print_scoring(scoring, virt_keys)

    return scoring

def simple_scoring(match, mismatch):

    scoring = {}
    virt_keys = 'ACDEFGHIKLMNPQRSTVWY'
    for i in range(len(virt_keys)):
        key = virt_keys[i]
        row_dict = {}
        for j in range(len(virt_keys)):
            virt_key = virt_keys[j]
            row_dict[virt_key] = match if i == j else mismatch
        scoring[key] = row_dict

    print_scoring(scoring, virt_keys)

    return scoring

def alignment_test_2(testfile, alignment_function, scoring, alignment_type):
    with open(testfile) as f:
        nucleotide_h = f.readline().rstrip()
        nucleotide_w = f.readline().rstrip()

    return alignment_function(nucleotide_h, nucleotide_w, scoring, alignment_type)

def test_verify(results, verifyfile):

    with open(verifyfile) as f:
        score = int(f.readline())
        val_1 = f.readline().rstrip()
        val_2 = f.readline().rstrip()

    print("calcualted result:")
    for i in range(3):
        print(results[i])

    print("expected result:")
    print(score)
    print(val_1)
    print(val_2)

    assert results[0] == score
    assert results[1] == val_1
    assert results[2] == val_2

def single_test(input_file, expected_results_file, scoring_table, alignment_type):
    print(input_file)
    print(expected_results_file)
    results = alignment_test_2(input_file, any_alignment.lcsbacktrack, scoring_table, alignment_type)
    test_verify(results, expected_results_file)

def global_tests():
    scoring = load_scoring("./BLOSUM62.txt")
    for in_f, out_f in [ \
                        ["test/globalalignment.txt", "test/globalalignment_expected_result.txt"], \
                        ["test/global_alignment_test_dataset_2.txt","test/global_alignment_test_dataset_2_expected_result.txt"], \
                        ["test/dataset_247_3.txt","test/dataset_247_3_expected_result.txt"] \
                        ]:
        single_test(in_f, out_f, scoring,"global")

def local_tests():
    scoring = load_scoring("./PAM250_scoring.txt")
    for in_f, out_f in [ \
                        ["test/localalignment.txt", "test/localalignment_expected_result.txt"], \
                        ["test/local_alignment_test_2.txt","test/local_alignment_test_2_expected_result.txt"], \
                        ["test/local_alignment_test_3.txt","test/local_alignment_test_3_expected_result.txt"], \
                        ["test/local_alignment_test_4.txt","test/local_alignment_test_4_expected_result.txt"], \
                        ["test/local_alignment_test_5.txt","test/local_alignment_test_5_expected_result.txt"], \
                        ]:
        single_test(in_f, out_f, scoring,"local")

def fitting_tests():
    scoring = simple_scoring(1, -1)
    for in_f, out_f in [ \
                        ["test/08_FittingAlignment/inputs/test1.txt", "test/08_FittingAlignment/outputs/test1.txt"], \
                        ["test/08_FittingAlignment/inputs/test2.txt", "test/08_FittingAlignment/outputs/test2.txt"], \
                        ["test/08_FittingAlignment/inputs/test3.txt", "test/08_FittingAlignment/outputs/test3.txt"], \
                        ["test/08_FittingAlignment/inputs/test4.txt", "test/08_FittingAlignment/outputs/test4.txt"], \
#                        ["test/08_FittingAlignment/inputs/test5.txt", "test/08_FittingAlignment/outputs/test5.txt"], \
                        ]:
        single_test(in_f, out_f, scoring,"fitting")

def overlap_tests():
    scoring = simple_scoring(1, -2)
    for in_f, out_f in [ \
                        ["test/overlapalignment.txt", "test/overlapalignment_expected_result.txt"], \
                        ["test/09_OverlapAlignment/inputs/test1.txt","test/09_OverlapAlignment/outputs/test1.txt"], \
                        ["test/09_OverlapAlignment/inputs/test2.txt","test/09_OverlapAlignment/outputs/test2.txt"], \
                        ["test/09_OverlapAlignment/inputs/test3.txt","test/09_OverlapAlignment/outputs/test3.txt"], \
                        ["test/09_OverlapAlignment/inputs/test4.txt","test/09_OverlapAlignment/outputs/test4.txt"], \
#                        ["test/09_OverlapAlignment/inputs/test5.txt","test/09_OverlapAlignment/outputs/test5.txt"], \
#                        ["test/09_OverlapAlignment/inputs/test6.txt","test/09_OverlapAlignment/outputs/test6.txt"], \
                        ["test/dataset_248_7.txt","test/dataset_248_7_expected_result.txt"], \
                        ]:
        single_test(in_f, out_f, scoring,"overlap")

def run_all_tests():

    global_tests()
    
    local_tests()

    fitting_tests()

    overlap_tests()

    return

if __name__ == '__main__':
    start = time.process_time()
    run_all_tests()
    end = time.process_time()
    print("Time: {0}".format(end-start))