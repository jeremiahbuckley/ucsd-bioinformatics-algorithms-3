#! /usr/bin/python3

import sys
import time
import math
import pdb

_prev_node_is_up_ = '↑'
_prev_node_is_left_ = '←'
_prev_node_is_diagonal_ = '↖︎'
_prev_node_is_upper_level_up_ = 'U'
_prev_node_is_lower_level_left_ = 'L'
_prev_node_is_middle_level_ = 'M'
_prev_node_is_zero_ = '0'


class AlignmentStrategy:
    def __init__(self, indel_penalty, scoring_from_file, scoring_filename = "", scoring_match_value = 1, scoring_mismatch_value = -1):
        self.indel_penalty = indel_penalty
        if scoring_from_file:
            self.scoring = self.load_scoring(scoring_filename)
        else:
            self.scoring = self.simple_scoring(scoring_match_value, scoring_mismatch_value)
        self.three_level_dag = False
        print("indel: {0}".format(str(self.indel_penalty)))
    
    def init_matrixes(self, nucleotide_h, nucleotide_w):
        lower_max_vals = []
        middle_max_vals = []
        upper_max_vals = []
        lower_backtrack = []
        middle_backtrack = []
        upper_backtrack = []

        for i in range(len(nucleotide_h)+1):
            lower_vals = []
            middle_vals = []
            upper_vals = []
            lower_bv = []
            middle_bv = []
            upper_bv = []
            for j in range(len(nucleotide_w)+1):
                lower_vals.append(0)
                middle_vals.append(0)
                upper_vals.append(0)
                lower_bv.append('x')
                middle_bv.append('x')
                upper_bv.append('x')
            lower_max_vals.append(lower_vals)
            middle_max_vals.append(middle_vals)
            upper_max_vals.append(upper_vals)
            lower_backtrack.append(lower_bv)
            middle_backtrack.append(middle_bv)
            upper_backtrack.append(upper_bv)

        for i in range(len(nucleotide_h)+1):
            lower_max_vals[i][0] = self.init_vertical_start_row_value(0, i)
            middle_max_vals[i][0] = self.init_vertical_start_row_value(1, i)
            upper_max_vals[i][0] = self.init_vertical_start_row_value(2, i)
            lower_backtrack[i][0] = _prev_node_is_zero_
            middle_backtrack[i][0] = _prev_node_is_zero_
            upper_backtrack[i][0] = _prev_node_is_zero_
        for i in range(len(nucleotide_w)+1):
            lower_max_vals[0][i] = self.init_horizontal_start_row_value(0, i)
            middle_max_vals[0][i] = self.init_horizontal_start_row_value(1, i)
            upper_max_vals[0][i] = self.init_horizontal_start_row_value(2, i)
            lower_backtrack[0][i] = _prev_node_is_zero_
            middle_backtrack[0][i] = _prev_node_is_zero_
            upper_backtrack[0][i] = _prev_node_is_zero_
            
        #print()
        #print_matrix(middle_max_vals)
        #print_matrix(backtrack)

        return lower_backtrack, middle_backtrack, upper_backtrack, lower_max_vals, middle_max_vals, upper_max_vals

    def init_horizontal_start_row_value(self, matrix_level, idx):
        return 0
    
    def init_vertical_start_row_value(self, matrix_level, idx):
        return 0

    def get_max_val_for_loc(self, lower_matrix_vals, middle_matrix_vals, upper_matrix_vals, v_idx, h_idx, match):
        return [-math.inf, max(middle_matrix_vals[0] - self.indel_penalty, middle_matrix_vals[1] - self.indel_penalty, middle_matrix_vals[2] + match), -math.inf]

    def found_new_max_value(self, current_value, max_value, loc, horizontal_leng, vertical_leng):
        return current_value >= max_value

    def pad_alignment_strings(self, align_h, align_w, height, width, nucleotide_h, nucleotide_w):
        return align_h, align_w

    def print_scoring(self, scoring, keys):
        print("      " + "    ".join([k for k in keys]))
        for i in range(len(keys)):
            val_str = ""
            for j in range(len(keys)):
                val_str += " " + str(scoring[keys[i]][keys[j]]).rjust(4)
            print("{0} {1}".format(keys[i], val_str))

    def load_scoring(self, scoring_file_loc):
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

        self.print_scoring(scoring, virt_keys)

        return scoring

    def simple_scoring(self, match, mismatch):

        scoring = {}
        virt_keys = 'ACDEFGHIKLMNPQRSTVWY'
        for i in range(len(virt_keys)):
            key = virt_keys[i]
            row_dict = {}
            for j in range(len(virt_keys)):
                virt_key = virt_keys[j]
                row_dict[virt_key] = match if i == j else mismatch
            scoring[key] = row_dict

        self.print_scoring(scoring, virt_keys)

        return scoring



class AlignmentStrategyGlobal(AlignmentStrategy):
    def init_horizontal_start_row_value(self, matrix_level, idx):
        if matrix_level == 1:
            return -1 * self.indel_penalty * idx
        return -math.inf
    
    def init_vertical_start_row_value(self, matrix_level, idx):
        if matrix_level == 1:
            return -1 * self.indel_penalty * idx
        return -math.inf

    def pad_alignment_strings(self, align_h, align_w, height, width, nucleotide_h, nucleotide_w):
        if height == 0 and width != 0:
            align_h = "-"*(width) + align_h
            align_w = nucleotide_w[0:width] + align_w
        if width == 0 and height != 0:
            align_h = nucleotide_h[0:height] + align_h
            align_w = "-"*(height) + align_w
        return align_h, align_w

class AlignmentStrategyLocal(AlignmentStrategy):
    def get_max_val_for_loc(self, lower_matrix_vals, middle_matrix_vals, upper_matrix_vals, v_idx, h_idx, match):
        m = max(middle_matrix_vals[0] - self.indel_penalty, middle_matrix_vals[1] - self.indel_penalty, middle_matrix_vals[2] + match)
        return [-math.inf, max(m, 0), -math.inf]

class AlignmentStrategyFitting(AlignmentStrategy):
    def init_horizontal_start_row_value(self, matrix_level, idx):
        if matrix_level == 1:
            return -1 * self.indel_penalty * idx
        return -math.inf

    def get_max_val_for_loc(self, lower_matrix_vals, middle_matrix_vals, upper_matrix_vals, v_idx, h_idx, match):
        mv = max(middle_matrix_vals[0] - self.indel_penalty, middle_matrix_vals[1] - self.indel_penalty, middle_matrix_vals[2] + match)
        if h_idx == 1:
            mv = max(0, mv)
        return [-math.inf, mv, -math.inf]

    def found_new_max_value(self, current_value, max_value, loc, horizontal_len, vertical_len):
        return loc[1] == horizontal_len and current_value >= max_value

class AlignmentStrategyOverlap(AlignmentStrategy):
    def init_horizontal_start_row_value(self, matrix_level, idx):
        if matrix_level == 1:
            return -1 * self.indel_penalty * idx
        return -math.inf

    def found_new_max_value(self, current_value, max_value, loc, horizontal_len, veritical_len):
        return loc[0] == veritical_len and current_value >= max_value

class AlignmentStrategyAffineGap(AlignmentStrategy):
    def __init__(self, indel_penalty, gap_initiation_penalty, scoring_from_file, scoring_filename = "", scoring_match_value = 1, scoring_mismatch_value = -1):
        super().__init__(indel_penalty, scoring_from_file, scoring_filename, scoring_match_value, scoring_mismatch_value)
        self.three_level_dag = True
        self.gap_initation_penatlty = gap_initiation_penalty
        print("gap init: {0}".format(str(self.gap_initation_penatlty)))

    def get_max_val_for_loc(self, lower_matrix_vals, middle_matrix_vals, upper_matrix_vals, v_idx, h_idx, match):
            # max_results = alignment_strategy.get_max_val_for_loc([lower_max_vals[i-1][j], lower_max_vals[i][j]], \
            #                                                      [middle_max_vals[i-1][j], \
            #                                                       middle_max_vals[i][j-1], \
            #                                                       middle_max_vals[i-1][j-1]], \
            #                                                      [upper_max_vals[i][j-1], upper_max_vals[i][j]], \
            #                                                      i, j, match )
        max_for_lower = max(lower_matrix_vals[0] - self.indel_penalty, middle_matrix_vals[0] - self.gap_initation_penatlty)
        max_for_upper = max(upper_matrix_vals[0] - self.indel_penalty, middle_matrix_vals[1] - self.gap_initation_penatlty)
        max_for_middle = max(max_for_lower, max_for_upper, middle_matrix_vals[2] + match)
        return [max_for_lower, max_for_middle, max_for_upper]

    def init_horizontal_start_row_value(self, matrix_level, idx):
        if matrix_level > 0:
            val = 0
            if idx > 0:
                val -= self.gap_initation_penatlty
            if idx > 1:
                val -= (self.indel_penalty * (idx-1))
            return val
        else:
            return -math.inf
    
    def init_vertical_start_row_value(self, matrix_level, idx):
        if matrix_level < 2:
            val = 0
            if idx > 0:
                val -= self.gap_initation_penatlty
            if idx > 1:
                val -= (self.indel_penalty * (idx-1))
            return val
        else:
            return -math.inf


def print_matrix(list_of_lists, str_h, str_v):
    height = len(list_of_lists)
    width = len(list_of_lists[0])

    if height > 80 or width > 80:
        print("Matrix will be too large to analyze. Not printing.")
        return

    #print()
    #print(len(str_v))
    #print(len(list_of_lists))
    print("           {0}".format("  ".join([ch.rjust(5) for ch in str_h])))
    for i in range(height):
        #print(list_of_lists[i])
        #print(str_v[i])
        #print()
        out_str = "  ".join([str(n).rjust(5) for n in list_of_lists[i]])
#        out_str = out_str.replace("↖︎", " ↖︎") #on the terminal this character is incorrectly right-justified by one-too-few spaces, this is not true n other formats
        if i == 0:
            print("    {0}".format(out_str))
        else:
            print("{0}   {1}".format(str_v[i-1], out_str))





def outputlcs(backtrack, height, width, nucleotide_h, nucleotide_w, alignment_strategy):

    align_h = ""
    align_w = ""
    #print(max_overall_loc)
    #print(nucleotide_h)
    #print(nucleotide_w)
    while height != 0 and width != 0:
        if backtrack[height][width] == _prev_node_is_zero_:
            height = 0
            width = 0
        elif backtrack[height][width] == _prev_node_is_up_:
            align_h = nucleotide_h[height-1] + align_h
            align_w = "-" + align_w
            height -= 1
        elif backtrack[height][width] == _prev_node_is_left_:
            align_h = "-" + align_h
            align_w = nucleotide_w[width-1] + align_w
            width -= 1
        else: #backtrack[height][width] == _prev_node_is_diagonal_:
            align_h = nucleotide_h[height-1] + align_h
            align_w = nucleotide_w[width-1] + align_w
            height -= 1
            width -= 1

    align_h, align_w = alignment_strategy.pad_alignment_strings(align_h, align_w, height, width, nucleotide_h, nucleotide_w)

    #print()
    return align_h, align_w

def lcsbacktrack(nucleotide_h, nucleotide_w, alignment_strategy):
    lower_backtrack, middle_backtrack, upper_backtrack, lower_max_vals, middle_max_vals, upper_max_vals = alignment_strategy.init_matrixes(nucleotide_h, nucleotide_w)

    max_overall = 0
    max_overall_loc = (0, 0)
    all_maxes = []
    for i in range(1, len(nucleotide_h)+1):
        for j in range(1, len(nucleotide_w)+1):
            #print("i:{0} j:{1}".format(str(i), str(j)))
            match = alignment_strategy.scoring[nucleotide_w[j-1]][nucleotide_h[i-1]]

            max_results = alignment_strategy.get_max_val_for_loc([lower_max_vals[i-1][j], lower_max_vals[i][j]], \
                                                                 [middle_max_vals[i-1][j], \
                                                                  middle_max_vals[i][j-1], \
                                                                  middle_max_vals[i-1][j-1]], \
                                                                 [upper_max_vals[i][j-1], upper_max_vals[i][j]], \
                                                                 i, j, match )
            middle_max_vals[i][j] = max_results[1]

            if alignment_strategy.three_level_dag:
                lower_max_vals[i][j] = max_results[0]
                upper_max_vals[i][j] = max_results[2]

            if middle_max_vals[i][j] == middle_max_vals[i-1][j] - alignment_strategy.indel_penalty:
                middle_backtrack[i][j] = _prev_node_is_up_
            elif middle_max_vals[i][j] == middle_max_vals[i][j-1] - alignment_strategy.indel_penalty:
                middle_backtrack[i][j] = _prev_node_is_left_
            elif middle_max_vals[i][j] == middle_max_vals[i-1][j-1] + match:
                middle_backtrack[i][j] = _prev_node_is_diagonal_
            else: # max_vals[i][j] == 0:
                middle_backtrack[i][j] = _prev_node_is_zero_

            new_max_val = alignment_strategy.found_new_max_value(middle_max_vals[i][j], max_overall, (i, j), len(nucleotide_w), len(nucleotide_h))
            if new_max_val:
                if middle_max_vals[i][j] > max_overall:
                    all_maxes = []

                #print(max_vals[i][j])
                max_overall = middle_max_vals[i][j]
                max_overall_loc = (i, j)
                all_maxes.append(max_overall_loc)
                # still have to figure out backtrack

        #print()
        #print_matrix(max_vals)
        #print_matrix(backtrack)
    

    #print()
    #print(all_maxes)
    f_score = middle_max_vals[max_overall_loc[0]][max_overall_loc[1]]
    out_h, out_w = outputlcs(middle_backtrack, max_overall_loc[0], max_overall_loc[1], nucleotide_h, nucleotide_w, alignment_strategy)

    if alignment_strategy.three_level_dag:
        print_matrix(lower_max_vals, nucleotide_w, nucleotide_h)
    print_matrix(middle_max_vals, nucleotide_w, nucleotide_h)
    if alignment_strategy.three_level_dag:
        print_matrix(upper_max_vals, nucleotide_w, nucleotide_h)
    print()
    if alignment_strategy.three_level_dag:
        print_matrix(lower_backtrack, nucleotide_w, nucleotide_h)
    print_matrix(middle_backtrack, nucleotide_w, nucleotide_h)
    if alignment_strategy.three_level_dag:
        print_matrix(middle_backtrack, nucleotide_w, nucleotide_h)
    print(f_score)
    print(out_h)
    print(out_w)

    return f_score, out_h, out_w




if __name__ == '__main__':
    start = time.process_time()

    if len(sys.argv) < 3:
        print("Expected input:\n 1 type of alignment <global,local,fitting,overlap>,\n 2 param,file. Format '<str_nucleotide_1>\n<str_nucleotide_2>")

    alignment_type = sys.argv[1]

    scoring = {}
    if alignment_type == "global":
        alignment_strategy = AlignmentStrategyGlobal(5, True, scoring_filename = "./BLOSUM62.txt")
    elif alignment_type == "local":
        alignment_strategy = AlignmentStrategyLocal(5, True, scoring_filename = "./PAM250_scoring.txt")
    elif alignment_type == "fitting":
        alignment_strategy = AlignmentStrategyFitting(1, False, scoring_match_value = 1, scoring_mismatch_value = -1)
    elif alignment_type == "overlap":
        alignment_strategy = AlignmentStrategyOverlap(2, False, scoring_match_value = 1, scoring_mismatch_value = -2)
    elif alignment_type == "affine_gap":
        alignment_strategy = AlignmentStrategyAffineGap(1, 11, True, scoring_filename = "./BLOSUM62.txt")
    else:
        raise ValueError("Unexpected alignment type: {0}".format(alignment_type))

    with open(sys.argv[2]) as f:
        nucleotide_h = f.readline().rstrip()
        nucleotide_w = f.readline().rstrip()

    results = lcsbacktrack(nucleotide_h, nucleotide_w, alignment_strategy)

    for i in range(3):
        print(results[i])

    end = time.process_time()
    print("Time: {0}".format(end-start))
