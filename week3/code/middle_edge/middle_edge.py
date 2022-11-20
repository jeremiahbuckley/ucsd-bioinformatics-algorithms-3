#! /usr/bin/python3

import sys
import time
import math
import pdb

import my_utils

_next_edge_right_ = "H"
_next_edge_down_ = "V"
_next_edge_diagonal_ = "D"

def find_middle_edge_vals_from_source(nucleotide_vertical, nucleotide_horizontal, middle_column, top_left, bottom_right, scoring):
    top = top_left[0]
    left = top_left[1]
    bottom = bottom_right[0]
    height = bottom - top
    v1_vals = []
    for i in range(height+1):
        v1_vals.append(i * scoring.indel_penalty)    

    if my_utils._debug_:
        pstr = ""
        for i in v1_vals:
            if len(pstr) == 0:
                pstr += str(i)
            else:
                pstr += ", " + str(i)
        print(pstr)

    v2_vals = v1_vals
    v2_edge_direction = []
    for i in v1_vals:
        v2_edge_direction.append([_next_edge_down_])
    for horizontal_idx in range(middle_column):
        v2_vals = []
        v2_edge_direction = []
        for i in range(height + 1):
            if i == 0:
                v2_vals.append(v1_vals[i] + scoring.indel_penalty)
                v2_edge_direction.append([_next_edge_right_])
            else:
                diag = scoring.get_match_val(nucleotide_vertical[top + i-1], nucleotide_horizontal[left + horizontal_idx])

                val_down = v2_vals[i-1] + scoring.indel_penalty
                val_diag = v1_vals[i-1] + diag
                val_right = v1_vals[i] + scoring.indel_penalty

                if val_down > val_diag and val_down > val_right:
                    v2_edge_direction.append([_next_edge_down_])
                elif val_diag > val_down and val_diag > val_right:
                    v2_edge_direction.append([_next_edge_diagonal_])
                elif val_right > val_down and val_right > val_diag:
                    v2_edge_direction.append([_next_edge_right_])
                elif val_down > val_right and val_down == val_diag:
                    v2_edge_direction.append([_next_edge_diagonal_, _next_edge_down_])
                elif val_down > val_diag and val_down == val_right:
                    v2_edge_direction.append([_next_edge_right_, _next_edge_down_])
                elif val_diag == val_right and val_diag > val_down:
                    v2_edge_direction.append([_next_edge_diagonal_, _next_edge_right_])
                elif val_diag == val_down and val_diag == val_right:
                    v2_edge_direction.append([_next_edge_diagonal_, _next_edge_right_, _next_edge_down_])
                else:
                    ValueError("shouldn't be here, all val_down / val_right / val_diag inequalities should be handled")

                v2_vals.append(max(val_down, val_diag, val_right))

        if my_utils._debug_:
            pstr = ""
            for i in v2_vals:
                if len(pstr) == 0:
                    pstr += str(i)
                else:
                    pstr += ", " + str(i)
            print(pstr)

        v1_vals = v2_vals

    return v2_vals, v2_edge_direction, v2_edge_direction


def find_middle_edge_data(nucleotide_vertical, nucleotide_horizontal, top_left, bottom_right, scoring):

    middle_column = math.floor((top_left[1] + bottom_right[1]) / 2) + top_left[1]

    if my_utils._verbose_:
        print("horz: {0}, length: {1}, mid: {2}".format(nucleotide_horizontal, top_left[1] + bottom_right[1], middle_column))
        print(middle_column)
        print()

    if my_utils._verbose_:
        print("forward")

    middle_edge_vals_from_source = find_middle_edge_vals_from_source(nucleotide_vertical, nucleotide_horizontal, middle_column, top_left, bottom_right, scoring)[0]
    if my_utils._verbose_:
        print("reverse")
    middle_edge_vals_backwards_from_sink_ret_vals = find_middle_edge_vals_from_source(nucleotide_vertical[::-1], nucleotide_horizontal[::-1], \
        len(nucleotide_horizontal) - middle_column, [len(nucleotide_vertical) - bottom_right[0], len(nucleotide_horizontal) - bottom_right[1]], [len(nucleotide_vertical) - 0, len(nucleotide_horizontal) - middle_column], scoring)
    middle_edge_vals_backwards_from_sink = middle_edge_vals_backwards_from_sink_ret_vals[0][::-1]
    middle_edge_direction_towards_sink = middle_edge_vals_backwards_from_sink_ret_vals[1][::-1]

    return middle_column, middle_edge_vals_from_source, middle_edge_vals_backwards_from_sink, middle_edge_direction_towards_sink

def find_middle_edge(nucleotide_vertical, nucleotide_horizontal, top_left, bottom_right, scoring):
    if my_utils._verbose_:
        scoring.print_scoring()
        print("virt: {0}".format(nucleotide_vertical))
        print("horz: {0}".format(nucleotide_horizontal))

    if len(nucleotide_horizontal) < 2:
        ValueError("not worth doing this for a 1-string horizontal")

    me_data = find_middle_edge_data(nucleotide_vertical, nucleotide_horizontal, top_left, bottom_right, scoring)
    
    middle_column = me_data[0]
    middle_edge_vals_from_source = me_data[1]
    middle_edge_vals_backwards_from_sink = me_data[2]
    middle_edge_direction_towards_sink = me_data[3]

    if len(middle_edge_vals_from_source) != len(middle_edge_vals_backwards_from_sink):
        ValueError("middle edge arrays must be same length")

    max_path = -math.inf
    max_node_idx = []
    max_node_idx_edge_directions = []
    for i in range(len(middle_edge_vals_from_source)):
        path_val = middle_edge_vals_from_source[i] + middle_edge_vals_backwards_from_sink[i]
        if path_val > max_path:
            max_path = path_val
            max_node_idx = []
            max_node_idx_edge_directions = []
        if path_val == max_path:
            max_node_idx.append([i, middle_column])
            max_node_idx_edge_directions = middle_edge_direction_towards_sink[i]
        
        # print(path_val)

    if my_utils._debug_:
        pstr = ""
        for i in range(len(middle_edge_vals_from_source)):
            if len(pstr) == 0:
                pstr += str(middle_edge_vals_from_source[i] + middle_edge_vals_backwards_from_sink[i])
            else:
                pstr += ", " + str(middle_edge_vals_from_source[i] + middle_edge_vals_backwards_from_sink[i])
        print("totals: " + pstr)


    if my_utils._verbose_:
        print("middle nodes")

    edges = []
    for node in max_node_idx:
        for j in max_node_idx_edge_directions:
            edge = []
            edge.append(node)

            if j == "H":
                edge.append([node[0], node[1] + 1])
            elif j == "V":
                edge.append([node[0] + 1, node[1]])
            elif j == "D":
                edge.append([node[0] + 1, node[1] + 1])
            else:
                ValueError("unexpected value: {0}".format(j))
            print()

            edges.append(edge)

    for edge in edges:
        print("{0} {1}".format(edge[0][0], edge[0][1]))
        print("{0} {1}".format(edge[1][0], edge[1][1]))
        print()

    return edges[0]

def find_middle_edge_for_test(nucleotide_vertical, nucleotide_horizontal, scoring):
    top_left = [0, 0]
    bottom_right = [len(nucleotide_vertical), len(nucleotide_horizontal)]

    edge = find_middle_edge(nucleotide_vertical, nucleotide_horizontal, top_left, bottom_right, scoring)

    testing_str = ""
    if len(edge) > 0:
        testing_str = "({0}, {1}) ({2}, {3})".format(edge[0][0], edge[0][1], edge[1][0], edge[1][1])
    
    return testing_str



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

        nucleotide_horizontal = f.readline().rstrip()
        nucleotide_vertical = f.readline().rstrip()
        # nucleotide_vertical = f.readline().rstrip()
        # nucleotide_horizontal = f.readline().rstrip()

    for a_idx in range(2, 3, 1):
        if len(sys.argv) > a_idx:
            if sys.argv[a_idx] == "-v":
                my_utils._verbose_ = True
            elif sys.argv[a_idx] == "-vv":
                my_utils._verbose_ = True
                my_utils._debug_ = True

    top_left = [0, 0]
    bottom_right = [len(nucleotide_vertical), len(nucleotide_horizontal)]

    results = find_middle_edge(nucleotide_vertical, nucleotide_horizontal, top_left, bottom_right, scoring)
 
    #print(results)

    end = time.process_time()
    print("Time: {0}".format(end-start))
