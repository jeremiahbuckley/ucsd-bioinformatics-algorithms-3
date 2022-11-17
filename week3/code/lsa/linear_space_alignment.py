#! /usr/bin/python3

import sys
import time
import math
import pdb

import my_utils

from middle_edge import find_middle_edge


def find_alignment(nucleotide_vertical, nucleotide_horizontal, top_left, bottom_right, scoring):
    path_ret = linear_space_alignment(nucleotide_vertical, nucleotide_horizontal, top_left, bottom_right, scoring, 0)
    path = path_ret[0]
    path_val = path_ret[1]

    if len(path) == 0:
        ValueError("unexpected path length = 0")
    print(path)

    previous_node = path[0][0]
    top_path = ""
    bottom_path = ""
    for edge in path:
        if (previous_node[0] != edge[0][0]) or (previous_node[1] != edge[0][1]):
            ValueError("unexpected break in path. previous ({0}, {1}), current({2}, {3}).".format(previous_node[0], previous_node[1], edge[0][0], edge[0][1]))
        next_node = edge[1]
        if next_node[0] == previous_node[0]:
            top_path += "-"
        else:
            top_path += nucleotide_vertical[next_node[0] - 1]

        if next_node[1] == previous_node[1]:
            bottom_path += "-"
        else:
            bottom_path += nucleotide_horizontal[next_node[1] - 1]
        previous_node = next_node
    
    return path_val, top_path, bottom_path

def linear_space_alignment(nucleotide_vertical, nucleotide_horizontal, top_left, bottom_right, scoring, otab):
    if (top_left[1] == bottom_right[1]):
        path = []
        prev_node = top_left
        for i in range(top_left[0] + 1, bottom_right[0] + 1):
            next_node = [i, top_left[1]]
            path.append([prev_node, next_node])
            prev_node = next_node
        if my_utils._debug_:
            print("  " * otab + str(path))
            print("  " * otab + str(scoring.indel_penalty * (bottom_right[0] - top_left[0])))
        return path, scoring.indel_penalty * (bottom_right[0] - top_left[0])
    if (top_left[0] == bottom_right[0]):
        path = []
        prev_node = top_left
        for i in range(top_left[1] + 1, bottom_right[1] + 1):
            next_node = [top_left[0], i]
            path.append([prev_node, next_node])
            prev_node = next_node
        if my_utils._debug_:
            print("  " * otab + str(path))
            print("  " * otab + str(scoring.indel_penalty * (bottom_right[1] - top_left[1])))
        return path, scoring.indel_penalty * (bottom_right[1] - top_left[1])


    middle_edge_ret = find_middle_edge(nucleotide_vertical, nucleotide_horizontal, top_left, bottom_right, scoring)
    middle_edge = middle_edge_ret[0]
    middle_edge_val = middle_edge_ret[1]

    if my_utils._debug_:
        print("  " * otab + str(middle_edge) + " (me)")
        print("  " * otab + str(middle_edge_val) + " (me)")
    left_path_ret = linear_space_alignment(nucleotide_vertical, nucleotide_horizontal, top_left, middle_edge[0], scoring, otab+1)
    left_path = left_path_ret[0]
    left_path_val = left_path_ret[1]
    path =[]
    for edge in left_path:
        path.append(edge)
    path.append(middle_edge)
    right_path_ret = linear_space_alignment(nucleotide_vertical, nucleotide_horizontal, middle_edge[1], bottom_right, scoring, otab+1)
    right_path = right_path_ret[0]
    right_path_val = right_path_ret[1]
    for edge in right_path:
        path.append(edge)

    if my_utils._debug_:
        print("  " * otab + str(path))
        print("  " * otab + str(left_path_val + middle_edge_val + right_path_val))

    return path, left_path_val + middle_edge_val + right_path_val




if __name__ == '__main__':
    start = time.process_time()

    if len(sys.argv) < 1:
        print("Expected input:\n[str: filename path]\n\nfile contents:\n[int: match reward] [int: mismatch penalty] [int: indel penalty],\n[string: nucleotide]\n[string: nucleotide]\n-v = verbose, -vv = debug")

    with open(sys.argv[1]) as f:
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
            scoring = my_utils.Scoring(match_reward, mismatch_penalty, indel_penalty)
        else:
            scoring_file = int_params.split()[0]            
            indel_penalty = -1 * int(int_params.split()[1])
            scoring = my_utils.Scoring(0, 0, indel_penalty, scoring_file)

        # nucleotide_horizontal = f.readline().rstrip()
        # nucleotide_vertical = f.readline().rstrip()
        nucleotide_vertical = f.readline().rstrip()
        nucleotide_horizontal = f.readline().rstrip()

    for a_idx in range(2, 3, 1):
        if len(sys.argv) > a_idx:
            if sys.argv[a_idx] == "-v":
                my_utils._verbose_ = True
            elif sys.argv[a_idx] == "-vv":
                my_utils._verbose_ = True
                my_utils._debug_ = True

    top_left = [0, 0]
    bottom_right = [len(nucleotide_vertical), len(nucleotide_horizontal)]

    results = find_alignment(nucleotide_vertical, nucleotide_horizontal, top_left, bottom_right, scoring)
 
    print(results[0])
    print(results[1])
    print(results[2])

    end = time.process_time()
    print("Time: {0}".format(end-start))
