#! /bin/python3

import pytest

from middle_edge import find_middle_edge_data, find_middle_edge, Scoring

@pytest.fixture
def sample_gaga_gat():
    nucleotide_vert = "GAGA"
    nucleotide_horz = "GAT"
    scoring = Scoring(1, -1, -2)

    top_left = [0,0]
    bottom_right = [len(nucleotide_vert), len(nucleotide_horz)]

    return nucleotide_vert, nucleotide_horz, top_left, bottom_right, scoring

@pytest.fixture
def sample_gat_gaga():
    nucleotide_vert = "GAT"
    nucleotide_horz = "GAGA"
    scoring = Scoring(1, -1, -2)

    top_left = [0,0]
    bottom_right = [len(nucleotide_vert), len(nucleotide_horz)]

    return nucleotide_vert, nucleotide_horz, top_left, bottom_right, scoring

def init_test_from_file(path):
    with open(path) as f:
        int_params = f.readline().rstrip()
        try:
            x = int(int_params.split()[0])
            gotint = True
        except ValueError:
            gotint = False
        if gotint:
            match_reward = int(int_params.split()[0])
            mismatch_penalty = -1 * int(int_params.split()[1])
            indel_penalty = -1 * int(int_params.split()[2])
            scoring = Scoring(match_reward, mismatch_penalty, indel_penalty)
        else:
            scoring_file = int_params.split()[0]            
            indel_penalty = -1 * int(int_params.split()[1])
            scoring = Scoring(0, 0, indel_penalty, scoring_file)

        nucleotide_horizontal = f.readline().rstrip()
        nucleotide_vertical = f.readline().rstrip()
        # nucleotide_vertical = f.readline().rstrip()
        # nucleotide_horizontal = f.readline().rstrip()

    results = find_middle_edge(nucleotide_vertical, nucleotide_horizontal, scoring)
    return results

def file_based_test(in_file, assert_file):
    results = init_test_from_file(in_file)
    with open(assert_file) as f:
        assert_val = f.readline().rstrip()

    assert results == assert_val

def test_sample_0(sample_gaga_gat):
    init_vals = sample_gaga_gat
    me_data = find_middle_edge_data(init_vals[0], init_vals[1], init_vals[2], init_vals[3], init_vals[4])

    middle_column = me_data[0]
    middle_edge_vals_from_source = me_data[1]
    middle_edge_vals_backwards_from_sink = me_data[2]
    middle_edge_direction_towards_sink = me_data[3]

    assert middle_column == 1

    assert_vals = [-2, 1, -1, -3, -5]
    assert len(middle_edge_vals_from_source) == len(assert_vals)
    for i in range(len(middle_edge_vals_from_source)):
        assert middle_edge_vals_from_source[i] == assert_vals[i]

    assert_vals = [-4, -2, -2, -1, -4]
    assert len(middle_edge_vals_backwards_from_sink) == len(assert_vals)
    for i in range(len(middle_edge_vals_backwards_from_sink)):
        assert middle_edge_vals_backwards_from_sink[i] == assert_vals[i]

    assert_vals = [['V'],['D'],['D'],['D'],['H']]
    assert len(middle_edge_direction_towards_sink) == len(assert_vals)
    for i in range(len(middle_edge_direction_towards_sink)):
        assert len(middle_edge_direction_towards_sink[i]) == len(assert_vals[i])
        for j in range(len(middle_edge_direction_towards_sink[i])):
            assert middle_edge_direction_towards_sink[i][j] == assert_vals[i][j]




def test_sample_0_swap(sample_gat_gaga):
    init_vals = sample_gat_gaga
    me_data = find_middle_edge_data(init_vals[0], init_vals[1], init_vals[2], init_vals[3], init_vals[4])

    middle_column = me_data[0]
    middle_edge_vals_from_source = me_data[1]
    middle_edge_vals_backwards_from_sink = me_data[2]
    middle_edge_direction_towards_sink = me_data[3]

    assert middle_column == 2

    assert_vals = [-4, -1, 2, 0]
    assert len(middle_edge_vals_from_source) == len(assert_vals)
    for i in range(len(middle_edge_vals_from_source)):
        assert middle_edge_vals_from_source[i] == assert_vals[i]

    assert_vals = [0, -2, -3, -4]
    assert len(middle_edge_vals_backwards_from_sink) == len(assert_vals)
    for i in range(len(middle_edge_vals_backwards_from_sink)):
        assert middle_edge_vals_backwards_from_sink[i] == assert_vals[i]

    assert_vals = [['D'],['D'],['D','H'],['H']]
    assert len(middle_edge_direction_towards_sink) == len(assert_vals)
    for i in range(len(middle_edge_direction_towards_sink)):
        assert len(middle_edge_direction_towards_sink[i]) == len(assert_vals[i])
        for j in range(len(middle_edge_direction_towards_sink[i])):
            assert middle_edge_direction_towards_sink[i][j] == assert_vals[i][j]

def test_txt1():
    root = "test/11_MiddleEdge/"
    r_in = "inputs/"
    r_out = "outputs/"

    file_based_test(root + r_in + "test1.txt.old", root + r_out + "test1.txt")

def test_txt2():
    root = "test/11_MiddleEdge/"
    r_in = "inputs/"
    r_out = "outputs/"

    file_based_test(root + r_in + "test2.txt.old", root + r_out + "test2.txt")

# def test_txt3():
#     root = "test/11_MiddleEdge/"
#     r_in = "inputs/"
#     r_out = "outputs/"

#     file_based_test(root + r_in + "test3.txt.old", root + r_out + "test3.txt")

def test_txt4():
    root = "test/11_MiddleEdge/"
    r_in = "inputs/"
    r_out = "outputs/"

    file_based_test(root + r_in + "test4.txt.old", root + r_out + "test4.txt")

def test_txt5():
    root = "test/11_MiddleEdge/"
    r_in = "inputs/"
    r_out = "outputs/"

    file_based_test(root + r_in + "test5.txt.old", root + r_out + "test5.txt")

# def test_txt6():
#     root = "test/11_MiddleEdge/"
#     r_in = "inputs/"
#     r_out = "outputs/"

#     file_based_test(root + r_in + "test6.txt.old", root + r_out + "test6.txt")

def test_sample5():
    root = "test/11_MiddleEdge/"
    r_in = "inputs/"
    r_out = "outputs/"

    file_based_test(root + r_in + "sample_5.txt", root + r_out + "sample_5.txt")

def test_sample5():
    root = "test/11_MiddleEdge/"
    r_in = "inputs/"
    r_out = "outputs/"

    file_based_test(root + r_in + "dataset_250_12_success.txt", root + r_out + "dataset_250_12_success.txt")

