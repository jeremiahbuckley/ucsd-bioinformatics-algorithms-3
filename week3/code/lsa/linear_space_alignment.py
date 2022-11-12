#! /usr/bin/python3

import sys
import time
import math
import pdb
from matrix_print import print_matrix, print_matrices

_prev_node_is_up_ = '↑'
_prev_node_is_left_ = '←'
_prev_node_is_diagonal_ = '↖︎'
_prev_node_is_upper_level_up_ = 'U'
_prev_node_is_lower_level_left_ = 'L'
_prev_node_is_middle_level_ = 'M'
_prev_node_is_zero_ = '0'

_segment_gap_ = [-1, -1]

class AlignmentStrategyGlobal:
    def __init__(self, indel_penalty, \
                 scoring_from_file, scoring_filename = "", \
                 scoring_match_value = 1, scoring_mismatch_value = -1):
        self.indel_penalty = indel_penalty
        if scoring_from_file:
            self.scoring = self.load_scoring(scoring_filename)
        else:
            self.scoring = self.simple_scoring(scoring_match_value, scoring_mismatch_value)
        self.three_level_dag = False
        print("indel: {0}".format(str(self.indel_penalty)))
        print()

    def create_matrices(self):
        backtrack = []
        max_val = []
        return [[max_val, backtrack]]

    def init_matrices(self, nucleotide_h, nucleotide_w):
        mv_bt_pairs = self.create_matrices()
        backtracks = []
        max_vals_list = []

        level = -1
        for mv_bt in mv_bt_pairs:
            max_vals = mv_bt[0]
            backtrack = mv_bt[1]
            level += 1
            for i in range(len(nucleotide_h)+1):
                vals = []
                bt_vals = []

                for j in range(len(nucleotide_w)+1):
                    vals.append(0)
                    bt_vals.append('x')
                max_vals.append(vals)
                backtrack.append(bt_vals)

            for i in range(len(nucleotide_h)+1):
                max_vals[i][0] = self.init_vertical_start_row_value(level, i)
                backtrack[i][0] = _prev_node_is_zero_
            for i in range(len(nucleotide_w)+1):
                max_vals[0][i] = self.init_horizontal_start_row_value(level, i)
                backtrack[0][i] = _prev_node_is_zero_

            max_vals_list.append(max_vals)
            backtracks.append(backtrack)

        #print()
        print_matrices(max_vals_list, backtracks, " " + nucleotide_h, " " + nucleotide_w)

        return backtracks, max_vals_list
        
    def get_last_node(self, v_len, h_len):
        return [v_len, h_len, 0] # default to only 1 level of matrix

    def init_horizontal_start_row_value(self, matrix_level, idx):
        return -1 * self.indel_penalty * idx

    def init_vertical_start_row_value(self, matrix_level, idx):
        return -1 * self.indel_penalty * idx

    def get_max_val_for_loc(self, max_val_matrices, \
                            v_idx, h_idx, match):
        middle_matrix_vals = max_val_matrices[0]
        return [-math.inf, \
                max(middle_matrix_vals[v_idx-1][h_idx] - self.indel_penalty, \
                    middle_matrix_vals[v_idx-1][h_idx-1] + match, \
                    middle_matrix_vals[v_idx][h_idx-1] - self.indel_penalty), \
                -math.inf]

    def assign_matrix_vals(self, max_vals_matrices, backtrack_matrices, max_results, i, j, match):
        middle_max_vals = max_vals_matrices[0]
        middle_backtrack = backtrack_matrices[0]

        middle_max_vals[i][j] = max_results[1]

        if middle_max_vals[i][j] == middle_max_vals[i-1][j] - self.indel_penalty:
            middle_backtrack[i][j] = _prev_node_is_up_
        elif middle_max_vals[i][j] == middle_max_vals[i][j-1] - self.indel_penalty:
            middle_backtrack[i][j] = _prev_node_is_left_
        elif middle_max_vals[i][j] == middle_max_vals[i-1][j-1] + match:
            middle_backtrack[i][j] = _prev_node_is_diagonal_
        elif middle_max_vals[i][j] == 0:
            middle_backtrack[i][j] = _prev_node_is_zero_
        else:
            ValueError("Unexpected value, middle_max_vals[{0}][{1}] == {2}". \
                       format(str(i), str(j), str(middle_max_vals[i][j])))

    def set_end(self, i, j, max_val_matrices, end_val, end_val_loc, all_ends, last_node_loc):
        if i > 0 and j > 0:
            new_end_val = self.found_new_end_value( \
                               max_val_matrices[0][i][j], end_val, [i, j, 0], last_node_loc)
            if new_end_val:
                if max_val_matrices[0][i][j] > end_val:
                    all_ends = []

                end_val_loc = [i, j, 0]
                end_val = max_val_matrices[0][i][j]
                all_ends.append(end_val_loc)
        return all_ends, end_val, end_val_loc

    def found_new_end_value(self, current_val, end_val, end_val_loc, last_node_loc):
        return current_val >= end_val and end_val_loc[1] == last_node_loc[1]

    def load_scoring(self, scoring_file_loc):
        print(scoring_file_loc)
        scoring_strs = []
        with open(scoring_file_loc) as f:
            scoring_strs = f.readlines()

        scoring = {}
        scoring_fast_index = []
        horz_keys = "".join(scoring_strs[0].split())
        virt_keys = ""

        for i in range(1, len(scoring_strs)):
            row = scoring_strs[i].split()
            key = row[0]
            virt_keys = virt_keys + key
            row_dict = {}
            row_fast_index = []
            for j in range(1, len(row)):
                horz_key = horz_keys[j-1]
                row_dict[horz_key] = int(row[j])
                row_fast_index.append(int(row[j]))
            scoring[key] = row_dict
            scoring_fast_index.append(row_fast_index)

        print_matrix(scoring_fast_index, virt_keys, horz_keys, 5)

        return scoring

    def simple_scoring(self, match, mismatch):

        scoring = {}
        virt_keys = 'ACDEFGHIKLMNPQRSTVWY'
        for v_key in virt_keys:
            row_dict = {}
            for h_key in virt_keys:
                row_dict[h_key] = match if v_key == h_key else mismatch
            scoring[v_key] = row_dict

        self.print_scoring(scoring, virt_keys)

        return scoring

    def outputlcs(self, backtrack_matrices, start_loc, end_loc, nucleotide_h, nucleotide_w):

        align_h = ""
        align_w = ""
        path = [[end_loc[0], end_loc[1]]]
        #print(max_location)
        #print(nucleotide_h)
        #print(nucleotide_w)
        backtrack = backtrack_matrices[0]
        height = end_loc[0]
        width = end_loc[1]
        while height != start_loc[0] or width != end_loc[1]:
            if backtrack[height][width] == _prev_node_is_up_:
                path.insert(0, [height-1, width])
                align_h = nucleotide_h[height-1] + align_h
                align_w = "-" + align_w
                height -= 1
            elif backtrack[height][width] == _prev_node_is_left_:
                path.insert(0, [height, width-1])
                align_h = "-" + align_h
                align_w = nucleotide_w[width-1] + align_w
                width -= 1
            elif backtrack[height][width] == _prev_node_is_diagonal_:
                path.insert(0, [height-1, width-1])
                align_h = nucleotide_h[height-1] + align_h
                align_w = nucleotide_w[width-1] + align_w
                height -= 1
                width -= 1
            else:
                ValueError("Unexpected value. Level {0}, Value: {1}". \
                           format(str(1), backtrack[height][width]))
        print("outputlcs path")
        print(path)

        return path, align_h, align_w

def find_end_edge(nucleotide_h, nucleotide_w, alignment_strategy, starting_loc, ending_loc):
    backtrack_matrices, max_vals_matrices = alignment_strategy.init_matrices(nucleotide_h, nucleotide_w)
    last_node_loc = ending_loc

    end_val = -math.inf
    end_val_loc = [0,0,0]
    all_possible_ends = []
    for i in range(starting_loc[0], ending_loc[0]+1):
        for j in range(starting_loc[1], ending_loc[1]+1):
            #print("i:{0} j:{1}".format(str(i), str(j)))
            match = alignment_strategy.scoring[nucleotide_w[j-1]][nucleotide_h[i-1]]

            max_results = alignment_strategy.get_max_val_for_loc(max_vals_matrices,
                                                                 i, j, match )
            alignment_strategy.assign_matrix_vals(max_vals_matrices, backtrack_matrices, \
                                                  max_results, i, j, match)

            all_possible_ends, end_val, end_val_loc = \
                alignment_strategy.set_end(i, j, max_vals_matrices, \
                                           end_val, end_val_loc, all_possible_ends, last_node_loc)

        #print()
        #print_matrices(max_vals_matrices, backtrack_matrices, nucleotide_h, nucleotide_w)

    #print()
    print_matrices(max_vals_matrices, backtrack_matrices, " " + nucleotide_h, " " + nucleotide_w)
    print("end: ")
    print(end_val_loc)

    prior_end_loc = []
    if backtrack_matrices[0][end_val_loc[0]][end_val_loc[1]] == _prev_node_is_left_:
        prior_end_loc = [end_val_loc[0], end_val_loc[1]-1]
    elif backtrack_matrices[0][end_val_loc[0]][end_val_loc[1]] == _prev_node_is_up_:
        prior_end_loc = [end_val_loc[0]-1, end_val_loc[1]]
    elif backtrack_matrices[0][end_val_loc[0]][end_val_loc[1]] == _prev_node_is_diagonal_:
        prior_end_loc = [end_val_loc[0]-1, end_val_loc[1]-1]
    return prior_end_loc, end_val_loc

def find_path(nucleotide_h, nucleotide_w, starting_loc, ending_loc, alignment_strategy):
    backtrack_matrices, max_vals_matrices = alignment_strategy.init_matrices(nucleotide_h, nucleotide_w)
    last_node_loc = ending_loc

    print(starting_loc)
    print(ending_loc)
    end_val = -math.inf
    end_val_loc = [0,0,0]
    all_possible_ends = []
    for i in range(starting_loc[0], ending_loc[0]+1):
        for j in range(starting_loc[1], ending_loc[1]+1):
            print("i:{0} j:{1}".format(str(i), str(j)))
            match = alignment_strategy.scoring[nucleotide_w[j-1]][nucleotide_h[i-1]]

            max_results = alignment_strategy.get_max_val_for_loc(max_vals_matrices,
                                                                 i, j, match )
            alignment_strategy.assign_matrix_vals(max_vals_matrices, backtrack_matrices, \
                                                  max_results, i, j, match)

            all_possible_ends, end_val, end_val_loc = \
                alignment_strategy.set_end(i, j, max_vals_matrices, \
                                           end_val, end_val_loc, all_possible_ends, last_node_loc)

        #print()
        print_matrices(max_vals_matrices, backtrack_matrices, " " + nucleotide_h, " " + nucleotide_w)

    #print()
    print("start: ")
    print (starting_loc)
    print("end: ")
    print(end_val_loc)

    path, align_h, align_w = alignment_strategy.outputlcs(backtrack_matrices, starting_loc, end_val_loc, nucleotide_h, nucleotide_w)
    return path, align_h, align_w


def find_middle(nucleotide_h, nucleotide_w, start_loc, end_loc, alignment_strategy):
    bisected_corner_w = math.ceil((start_loc[1] + end_loc[1]) / 2)

    bisected_corner = [end_loc[0], bisected_corner_w]

    # if (start_loc[1] +1) >= bisected_corner[1]:
    #     return [start_loc, bisected_corner]
    # else:
    end_edge_start, end_edge_end = find_end_edge(nucleotide_h, nucleotide_w, alignment_strategy, start_loc, bisected_corner)
    print(end_edge_start)
    print(end_edge_end)
    return [end_edge_start, end_edge_end]

def find_alignment_segment(nucleotide_h, nucleotide_w, segment, alignment_strategy):
    start_loc = segment[0]
    end_loc = segment[1]

    if (start_loc[0] == end_loc[0] and start_loc[1] == end_loc[1]):
        print("same")
        print(start_loc)
        return [start_loc]
    elif (start_loc[1] + 1) >= end_loc[1]:
        print("neighbors")
        path_segment = find_path(nucleotide_h, nucleotide_w, start_loc, end_loc, alignment_strategy)
        print(path_segment)
        return path_segment
    else:
        print("distant")
        #TODO - if middle + 1 == start, find_middle could output all the info for path1 to be created using outputlcs
        middle_nodes = find_middle(nucleotide_h, nucleotide_w, start_loc, end_loc, alignment_strategy)
        print(middle_nodes)
        path_segment1 = find_alignment_segment(nucleotide_h, nucleotide_w, [start_loc, middle_nodes[0]], alignment_strategy)
        path_segment2 = find_alignment_segment(nucleotide_h, nucleotide_w, [middle_nodes[1], end_loc], alignment_strategy)
        print(path_segment1)
        print(path_segment2)
        print("oin")
        p1 = path_segment1[0:len(path_segment1)]
        p2 = path_segment2[0:len(path_segment2)]
        print(p1)
        print(p2)
        full_path = path_segment1[0:len(path_segment1)] + path_segment2[0:len(path_segment2)]
        print(full_path)
        return full_path




def find_alignment(nucleotide_h, nucleotide_w, alignment_strategy):
    start_loc = [1, 1]
    end_loc = [len(nucleotide_h), len(nucleotide_w)]
    return find_alignment_segment(nucleotide_h, nucleotide_w, [start_loc, end_loc], alignment_strategy)




if __name__ == '__main__':
    start = time.process_time()

    if len(sys.argv) < 3:
        print("Expected input:\n 1 type of alignment <global>,\n 2 param,file. Format '<str_nucleotide_1>\n<str_nucleotide_2>")

    alignment_type = sys.argv[1]

    scoring = {}
    if alignment_type == "global":
        alignment_strategy = AlignmentStrategyGlobal(5, \
                                 True, scoring_filename = "./BLOSUM62.txt")
    else:
        raise ValueError("Unexpected alignment type: {0}".format(alignment_type))

    with open(sys.argv[2]) as f:
        nucleotide_h = f.readline().rstrip()
        nucleotide_w = f.readline().rstrip()

    results = find_alignment(nucleotide_h, nucleotide_w, alignment_strategy)
 
    #print(results)

    end = time.process_time()
    print("Time: {0}".format(end-start))
