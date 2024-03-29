#! /usr/bin/python3

import sys
import time
#from week2.any_alignment import AlignmentStrategyGlobal

import any_alignment

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

def single_test(input_file, expected_results_file, alignment_type, alignment_strategy):
    print(input_file)
    print(expected_results_file)


    with open(input_file) as f:
        nucleotide_h = f.readline().rstrip()
        nucleotide_w = f.readline().rstrip()

    results = any_alignment.lcsbacktrack(nucleotide_h, nucleotide_w, alignment_strategy)
    test_verify(results, expected_results_file)

def global_tests():
    for in_f, out_f in [ \
                        ["test/globalalignment.txt", "test/globalalignment_expected_result.txt"], \
                        ["test/global_alignment_test_dataset_2.txt","test/global_alignment_test_dataset_2_expected_result.txt"], \
                        ["test/dataset_247_3.txt","test/dataset_247_3_expected_result.txt"] \
                        ]:
        single_test(in_f, out_f,"global", any_alignment.AlignmentStrategyGlobal(5, True, scoring_filename = "./BLOSUM62.txt"))

def local_tests():
    for in_f, out_f in [ \
                        ["test/localalignment.txt", "test/localalignment_expected_result.txt"], \
                        ["test/local_alignment_test_2.txt","test/local_alignment_test_2_expected_result.txt"], \
                        ["test/local_alignment_test_3.txt","test/local_alignment_test_3_expected_result.txt"], \
                        ["test/local_alignment_test_4.txt","test/local_alignment_test_4_expected_result.txt"], \
                        ["test/local_alignment_test_5.txt","test/local_alignment_test_5_expected_result.txt"], \
                        ]:
        single_test(in_f, out_f,"local", any_alignment.AlignmentStrategyLocal(5, True, scoring_filename = "./PAM250_scoring.txt"))

def fitting_tests():
    for in_f, out_f in [ \
                        ["test/08_FittingAlignment/inputs/test1.txt", "test/08_FittingAlignment/outputs/test1.txt"], \
                        ["test/08_FittingAlignment/inputs/test2.txt", "test/08_FittingAlignment/outputs/test2.txt"], \
                        ["test/08_FittingAlignment/inputs/test3.txt", "test/08_FittingAlignment/outputs/test3.txt"], \
                        ["test/08_FittingAlignment/inputs/test4.txt", "test/08_FittingAlignment/outputs/test4.txt"], \
#                        ["test/08_FittingAlignment/inputs/test5.txt", "test/08_FittingAlignment/outputs/test5.txt"], \
                        ]:
        single_test(in_f, out_f,"fitting", any_alignment.AlignmentStrategyFitting(1, False, scoring_match_value = 1, scoring_mismatch_value = -1))

def overlap_tests():
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
        single_test(in_f, out_f,"overlap", any_alignment.AlignmentStrategyOverlap(2, False, scoring_match_value = 1, scoring_mismatch_value = -2))

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