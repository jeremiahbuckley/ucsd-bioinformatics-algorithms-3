#! /usr/bin/python3

import sys
import time
#from week2.any_alignment import AlignmentStrategyGlobal

import affine_gap as any_alignment
import space_efficient_global_match

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

def test_verify_space_efficient(results, verifyfile):
    with open(verifyfile) as f:
        edge = f.readline()

    print("calcualted result:")
    print(results)

    print("expected result:")
    print(edge)

    assert results == edge


def single_test(input_file, expected_results_file, alignment_strategy):
    print(input_file)
    print(expected_results_file)


    with open(input_file) as f:
        nucleotide_h = f.readline().rstrip()
        nucleotide_w = f.readline().rstrip()

    results = any_alignment.lcsbacktrack(nucleotide_h, nucleotide_w, alignment_strategy)
    test_verify(results, expected_results_file)

def single_test_space_efficient(input_file, expected_results_file, alignment_strategy):
    print(input_file)
    print(expected_results_file)


    with open(input_file) as f:
        nucleotide_h = f.readline().rstrip()
        nucleotide_w = f.readline().rstrip()

    results = space_efficient_global_match.find_middle_edge(nucleotide_h, nucleotide_w, alignment_strategy)
    test_verify_space_efficient(results, expected_results_file)

def global_tests():
    for in_f, out_f in [ \
                        ["test/globalalignment.txt", "test/globalalignment_expected_result.txt"], \
                        ["test/global_alignment_test_dataset_2.txt","test/global_alignment_test_dataset_2_expected_result.txt"], \
                        ["test/dataset_247_3.txt","test/dataset_247_3_expected_result.txt"] \
                        ]:
        single_test(in_f, out_f, any_alignment.AlignmentStrategyGlobal(5, True, scoring_filename = "./BLOSUM62.txt"))

def local_tests():
    for in_f, out_f in [ \
                        ["test/localalignment.txt", "test/localalignment_expected_result.txt"], \
                        ["test/local_alignment_test_2.txt","test/local_alignment_test_2_expected_result.txt"], \
                        ["test/local_alignment_test_3.txt","test/local_alignment_test_3_expected_result.txt"], \
                        ["test/local_alignment_test_4.txt","test/local_alignment_test_4_expected_result.txt"], \
                        ["test/local_alignment_test_5.txt","test/local_alignment_test_5_expected_result.txt"], \
                        ]:
        single_test(in_f, out_f, any_alignment.AlignmentStrategyLocal(5, True, scoring_filename = "./PAM250_scoring.txt"))

def fitting_tests():
    for in_f, out_f in [ \
                        ["test/08_FittingAlignment/inputs/test1.txt", "test/08_FittingAlignment/outputs/test1.txt"], \
                        ["test/08_FittingAlignment/inputs/test2.txt", "test/08_FittingAlignment/outputs/test2.txt"], \
                        ["test/08_FittingAlignment/inputs/test3.txt", "test/08_FittingAlignment/outputs/test3.txt"], \
                        ["test/08_FittingAlignment/inputs/test4.txt", "test/08_FittingAlignment/outputs/test4.txt"], \
#                        ["test/08_FittingAlignment/inputs/test5.txt", "test/08_FittingAlignment/outputs/test5.txt"], \
                        ]:
        single_test(in_f, out_f, any_alignment.AlignmentStrategyFitting(1, False, scoring_match_value = 1, scoring_mismatch_value = -1))

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
        single_test(in_f, out_f, any_alignment.AlignmentStrategyOverlap(2, False, scoring_match_value = 1, scoring_mismatch_value = -2))

def affine_gap_tests():
    for in_f, out_f in [ \
                        ["test/dataset_249_8.txt", "test/dataset_249_8_expected_result.txt"], \
                        ["test/dataset_249_8_2.txt", "test/dataset_249_8_2_expected_result.txt"], \
                        ["test/affine_gap.txt", "test/affine_gap_expected_result.txt"], \
                        ["test/affine_gap_2.txt", "test/affine_gap_2_expected_result.txt"], \
                        ["test/10_AffineGapPenalties/inputs/test1.txt","test/10_AffineGapPenalties/outputs/test1.txt"], \
                        ["test/10_AffineGapPenalties/inputs/test2.txt","test/10_AffineGapPenalties/outputs/test2.txt"], \
                        ["test/10_AffineGapPenalties/inputs/test3.txt","test/10_AffineGapPenalties/outputs/test3.txt"], \
                        ["test/10_AffineGapPenalties/inputs/test4.txt","test/10_AffineGapPenalties/outputs/test4.txt"], \
                        ["test/10_AffineGapPenalties/inputs/test5.txt","test/10_AffineGapPenalties/outputs/test5.txt"], \
                        ["test/10_AffineGapPenalties/inputs/test6.txt","test/10_AffineGapPenalties/outputs/test6.txt"], \
                        ["test/10_AffineGapPenalties/inputs/test7.txt","test/10_AffineGapPenalties/outputs/test7.txt"], \
                        ["test/10_AffineGapPenalties/inputs/test8.txt","test/10_AffineGapPenalties/outputs/test8.txt"], \
                        ["test/10_AffineGapPenalties/inputs/sample.txt","test/10_AffineGapPenalties/outputs/sample.txt"], \
                        ["test/10_AffineGapPenalties/inputs/sample_2.txt","test/10_AffineGapPenalties/outputs/sample_2.txt"], \
                        ["test/10_AffineGapPenalties/inputs/sample_3.txt","test/10_AffineGapPenalties/outputs/sample_3.txt"], \
                        ]:
        single_test(in_f, out_f, any_alignment.AlignmentStrategyAffineGap(1, 11, True, scoring_filename = "./BLOSUM62.txt"))
        # single_test(in_f, out_f, any_alignment.AlignmentStrategyAffineGap(-4, 0, True, scoring_filename = "./BLOSUM62.txt"))

def affine_gap_stop_and_think():
    for in_f, out_f in [ \
                        ["test/affine_gap_stop_and_think_sample.txt", "test/affine_gap_stop_and_think_4_4_expected_result.txt"], \
                        ]:
        single_test(in_f, out_f, any_alignment.AlignmentStrategyAffineGap(4, 4, True, scoring_filename = "./BLOSUM62.txt"))
    for in_f, out_f in [ \
                        ["test/affine_gap_stop_and_think_sample.txt", "test/affine_gap_stop_and_think_1_11_expected_result.txt"], \
                        ]:
        single_test(in_f, out_f, any_alignment.AlignmentStrategyAffineGap(1, 11, True, scoring_filename = "./BLOSUM62.txt"))
    for in_f, out_f in [ \
                        ["test/affine_gap_stop_and_think_sample.txt", "test/affine_gap_stop_and_think_2_10_expected_result.txt"], \
                        ]:
        single_test(in_f, out_f, any_alignment.AlignmentStrategyAffineGap(2, 10, True, scoring_filename = "./BLOSUM62.txt"))
    for in_f, out_f in [ \
                        ["test/affine_gap_stop_and_think_sample.txt", "test/affine_gap_stop_and_think_1_15_expected_result.txt"], \
                        ]:
        single_test(in_f, out_f, any_alignment.AlignmentStrategyAffineGap(1, 15, True, scoring_filename = "./BLOSUM62.txt"))
    for in_f, out_f in [ \
                        ["test/affine_gap_stop_and_think_sample.txt", "test/affine_gap_stop_and_think_2_4_expected_result.txt"], \
                        ]:
        single_test(in_f, out_f, any_alignment.AlignmentStrategyAffineGap(2, 4, True, scoring_filename = "./BLOSUM62.txt"))

def affine_gap_local_tests():
    for in_f, out_f in [ \
                        ["test/affine_gap_local.txt","test/affine_gap_local_expected_result.txt"], \
                        ]:
        single_test(in_f, out_f, any_alignment.AlignmentStrategyAffineGapLocal(1, 11, True, scoring_filename = "./BLOSUM62.txt"))

def space_efficient_global_match_tests():
    for in_f, out_f in [ \
#                        ["test/11_MiddleEdge/inputs/sample.txt","test/11_MiddleEdge/outputs/sample.txt"], \
                        ["test/11_MiddleEdge/inputs/sample_2.txt","test/11_MiddleEdge/outputs/sample_2.txt"], \
                        ["test/11_MiddleEdge/inputs/sample_3.txt","test/11_MiddleEdge/outputs/sample_3.txt"], \
                        ["test/11_MiddleEdge/inputs/sample_4.txt","test/11_MiddleEdge/outputs/sample_4.txt"], \
                        ["test/11_MiddleEdge/inputs/sample_5.txt","test/11_MiddleEdge/outputs/sample_5.txt"], \
                        ["test/11_MiddleEdge/inputs/sample_6.txt","test/11_MiddleEdge/outputs/sample_6.txt"], \
                        ["test/11_MiddleEdge/inputs/sample_7.txt","test/11_MiddleEdge/outputs/sample_7.txt"], \
                        ]:
        single_test_space_efficient(in_f, out_f, space_efficient_global_match.AlignmentStrategyGlobal(5, True, scoring_filename = "./BLOSUM62.txt"))



def run_all_tests():

    # global_tests()
    
    # local_tests()

    # fitting_tests()

    # overlap_tests()

    # affine_gap_tests()

    # affine_gap_local_tests()

    # affine_gap_stop_and_think()

    space_efficient_global_match_tests()

    return

if __name__ == '__main__':
    start = time.process_time()
    run_all_tests()
    end = time.process_time()
    print("Time: {0}".format(end-start))