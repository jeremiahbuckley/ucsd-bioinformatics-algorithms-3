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

    def create_matrices(self):
        backtrack = []
        max_val = []
        return [[max_val, backtrack]]

    def init_matrixes(self, nucleotide_h, nucleotide_w):
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
        #print_matrix(middle_max_vals)
        #print_matrix(backtrack)

        return backtracks, max_vals_list

    def init_horizontal_start_row_value(self, matrix_level, idx):
        return 0

    def init_vertical_start_row_value(self, matrix_level, idx):
        return 0

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
            ValueError("Unexpected value, middle_max_vals[{0}][{1}] == {2}".format(str(i), str(j), str(middle_max_vals[i][j])))


    def found_new_max_value(self, current_value, max_value, loc, horizontal_len, vertical_len):
        return current_value >= max_value

    def pad_alignment_strings(self, align_h, align_w, height, width, nucleotide_h, nucleotide_w):
        return align_h, align_w

    def print_scoring(self, scoring, keys):
        print("      " + "    ".join(list(keys)))
        for v_key in keys:
            val_str = ""
            for h_key in keys:
                val_str += " " + str(scoring[v_key][h_key]).rjust(4)
            print("{0} {1}".format(v_key, val_str))

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
        for v_key in virt_keys:
            row_dict = {}
            for h_key in virt_keys:
                row_dict[h_key] = match if v_key == h_key else mismatch
            scoring[v_key] = row_dict

        self.print_scoring(scoring, virt_keys)

        return scoring

    def outputlcs(self, backtrack_matrices, max_location, nucleotide_h, nucleotide_w):

        align_h = ""
        align_w = ""
        #print(max_location)
        #print(nucleotide_h)
        #print(nucleotide_w)
        backtrack = backtrack_matrices[0]
        height = max_location[0]
        width = max_location[1]
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
            elif backtrack[height][width] == _prev_node_is_diagonal_:
                align_h = nucleotide_h[height-1] + align_h
                align_w = nucleotide_w[width-1] + align_w
                height -= 1
                width -= 1
            else:
                ValueError("Unexpected value. Level {0}, Value: {1}".format(str(1), middle_backtrack[height][width]))

        align_h, align_w = self.pad_alignment_strings(align_h, align_w, height, width, nucleotide_h, nucleotide_w)

        return align_h, align_w

class AlignmentStrategyGlobal(AlignmentStrategy):
    def init_horizontal_start_row_value(self, matrix_level, idx):
        return -1 * self.indel_penalty * idx

    def init_vertical_start_row_value(self, matrix_level, idx):
        return -1 * self.indel_penalty * idx

    def pad_alignment_strings(self, align_h, align_w, height, width, nucleotide_h, nucleotide_w):
        if height == 0 and width != 0:
            align_h = "-"*(width) + align_h
            align_w = nucleotide_w[0:width] + align_w
        if width == 0 and height != 0:
            align_h = nucleotide_h[0:height] + align_h
            align_w = "-"*(height) + align_w
        return align_h, align_w

class AlignmentStrategyLocal(AlignmentStrategy):
    def get_max_val_for_loc(self, max_val_matrices, v_idx, h_idx, match):
        middle_matrix_vals = max_val_matrices[0]
        m = max(middle_matrix_vals[v_idx-1][h_idx] - self.indel_penalty, \
                middle_matrix_vals[v_idx-1][h_idx-1] + match, \
                middle_matrix_vals[v_idx][h_idx-1] - self.indel_penalty)
        return [-math.inf, max(m, 0), -math.inf]

class AlignmentStrategyFitting(AlignmentStrategy):
    def init_horizontal_start_row_value(self, matrix_level, idx):
        return -1 * self.indel_penalty * idx

    def get_max_val_for_loc(self, max_val_matrices, v_idx, h_idx, match):
        middle_matrix_vals = max_val_matrices[0]
        m = max(middle_matrix_vals[v_idx-1][h_idx] - self.indel_penalty, \
                middle_matrix_vals[v_idx-1][h_idx-1] + match, \
                middle_matrix_vals[v_idx][h_idx-1] - self.indel_penalty)
        if h_idx == 1:
            m = max(0, m)
        return [-math.inf, m, -math.inf]

    def found_new_max_value(self, current_value, max_value, loc, horizontal_len, vertical_len):
        return loc[1] == horizontal_len and current_value >= max_value

class AlignmentStrategyOverlap(AlignmentStrategy):
    def init_horizontal_start_row_value(self, matrix_level, idx):
        return -1 * self.indel_penalty * idx

    def found_new_max_value(self, current_value, max_value, loc, horizontal_len, vertical_len):
        return loc[0] == vertical_len and current_value >= max_value

class AlignmentStrategyAffineGap(AlignmentStrategy):
    def __init__(self, indel_penalty, gap_initiation_penalty, scoring_from_file, scoring_filename = "", scoring_match_value = 1, scoring_mismatch_value = -1):
        super().__init__(indel_penalty, scoring_from_file, scoring_filename, scoring_match_value, scoring_mismatch_value)
        self.three_level_dag = True
        self.gap_initiation_penalty = gap_initiation_penalty
        print("gap init: {0}".format(str(self.gap_initiation_penalty)))

    def create_matrices(self):
        mv_bt_pairs = []
        for i in range(3):
            max_vals = []
            backtrack = []
            mv_bt_pairs.append([max_vals, backtrack])
        return mv_bt_pairs
        
    def get_max_val_for_loc(self, max_val_matrices, v_idx, h_idx, match):
        lower_matrix_vals = max_val_matrices[0]
        middle_matrix_vals = max_val_matrices[1]
        upper_matrix_vals = max_val_matrices[2]
        max_for_lower = max(lower_matrix_vals[v_idx][h_idx-1] - self.indel_penalty, \
                            middle_matrix_vals[v_idx][h_idx-1] - self.gap_initiation_penalty)
        max_for_upper = max(upper_matrix_vals[v_idx-1][h_idx] - self.indel_penalty, \
                            middle_matrix_vals[v_idx-1][h_idx] - self.gap_initiation_penalty)
        max_for_middle = max(max_for_lower, max_for_upper, \
                             middle_matrix_vals[v_idx-1][h_idx-1] + match)
        return [max_for_lower, max_for_middle, max_for_upper]

    def init_horizontal_start_row_value(self, matrix_level, idx):
        if matrix_level > 0:
            val = 0
            if idx > 0:
                val -= (self.gap_initiation_penalty + (self.indel_penalty * (idx - 1)))
            return val

        return -math.inf

    def init_vertical_start_row_value(self, matrix_level, idx):
        if matrix_level < 2:
            val = 0
            if idx > 0:
                val -= (self.gap_initiation_penalty + (self.indel_penalty * (idx - 1)))
            return val

        return -math.inf

    def pad_alignment_strings(self, align_h, align_w, height, width, nucleotide_h, nucleotide_w):
        print("padding " + str(height) + " " + str(width))
        if height == 0 and width != 0:
            align_h = "-"*(width) + align_h
            align_w = nucleotide_w[0:width] + align_w
        if width == 0 and height != 0:
            align_h = nucleotide_h[0:height] + align_h
            align_w = "-"*(height) + align_w
        return align_h, align_w

    def found_new_max_value(self, current_value, max_value, loc, horizontal_len, vertical_len):
        possible = super().found_new_max_value(current_value, max_value, loc, horizontal_len, vertical_len)
        if loc[1] == horizontal_len and loc[0] == vertical_len:
            return True
        return possible

    def outputlcs(self, backtrack_matrices, max_location, nucleotide_h, nucleotide_w):

        align_h = ""
        align_w = ""
        #print(max_location)
        #print(nucleotide_h)
        #print(nucleotide_w)
        lower_backtrack = backtrack_matrices[0]
        middle_backtrack = backtrack_matrices[1]
        upper_backtrack = backtrack_matrices[2]
        height = max_location[0]
        width = max_location[1]
        level = max_location[2]
        while height != 0 and width != 0:
            # pdb.set_trace()
            if level == 1:
                if middle_backtrack[height][width] == _prev_node_is_zero_:
                    height = 0
                    width = 0
                elif middle_backtrack[height][width] == _prev_node_is_upper_level_up_:
                    level = 2
                elif middle_backtrack[height][width] == _prev_node_is_lower_level_left_:
                    level = 0
                elif middle_backtrack[height][width] == _prev_node_is_diagonal_:
                    align_h = nucleotide_h[height-1] + align_h
                    align_w = nucleotide_w[width-1] + align_w
                    height -= 1
                    width -= 1
                else:
                    ValueError("Unexpected value. Level {0}, Value: {1}" \
                               .format(str(level), middle_backtrack[height][width]))
            elif level == 0:
                align_h = "-" + align_h
                align_w = nucleotide_w[width-1] + align_w
                if lower_backtrack[height][width] == _prev_node_is_middle_level_:
                    width -= 1
                    level = 1
                elif lower_backtrack[height][width] == _prev_node_is_left_:
                    width -= 1
                else:
                    ValueError("Unexpected value. Level {0}, Value: {1}" \
                               .format(str(level), middle_backtrack[height][width]))
            elif level == 2:
                align_h = nucleotide_h[height-1] + align_h
                align_w = "-" + align_w
                if upper_backtrack[height][width] == _prev_node_is_middle_level_:
                    height -= 1
                    level = 1
                elif upper_backtrack[height][width] == _prev_node_is_up_:
                    height -= 1
                else:
                    ValueError("Unexpected value. Level {0}, Value: {1}" \
                               .format(str(level), middle_backtrack[height][width]))
            else:
                ValueError("Unexpected value. Level {0}, Value: {1}" \
                           .format(str(level), middle_backtrack[height][width]))

        align_h, align_w = self.pad_alignment_strings(align_h, align_w, height, width, \
                                                      nucleotide_h, nucleotide_w)

        return align_h, align_w

    def assign_matrix_vals(self, max_vals_matrices, backtrack_matrices, max_results, i, j, match):
        lower_max_vals = max_vals_matrices[0]
        middle_max_vals = max_vals_matrices[1]
        upper_max_vals = max_vals_matrices[2]
        lower_backtrack = backtrack_matrices[0]
        middle_backtrack = backtrack_matrices[1]
        upper_backtrack = backtrack_matrices[2]

        middle_max_vals[i][j] = max_results[1]
        lower_max_vals[i][j] = max_results[0]
        upper_max_vals[i][j] = max_results[2]

        if middle_max_vals[i][j] == middle_max_vals[i-1][j-1] + match:
            middle_backtrack[i][j] = _prev_node_is_diagonal_
        elif middle_max_vals[i][j] == upper_max_vals[i][j]:
            middle_backtrack[i][j] = _prev_node_is_upper_level_up_
        elif middle_max_vals[i][j] == lower_max_vals[i][j]:
            middle_backtrack[i][j] = _prev_node_is_lower_level_left_
        elif middle_max_vals[i][j] == 0:
            middle_backtrack[i][j] = _prev_node_is_zero_
        else:
            ValueError("Unexpected value, middle_max_vals[{0}][{1}] == {2}".format(str(i), str(j), str(middle_max_vals[i][j])))

        if upper_max_vals[i][j] == upper_max_vals[i-1][j] - self.indel_penalty:
            upper_backtrack[i][j] = _prev_node_is_up_
        elif upper_max_vals[i][j] == middle_max_vals[i-1][j] - self.gap_initiation_penalty:
            upper_backtrack[i][j] = _prev_node_is_middle_level_
        else:
            ValueError("Unexpected value, uppper_max_vals[{0}][{1}] == {2}".format(str(i), str(j), str(upper_max_vals[i][j])))

        if lower_max_vals[i][j] == lower_max_vals[i][j-1] - self.indel_penalty:
            lower_backtrack[i][j] = _prev_node_is_left_
        elif lower_max_vals[i][j] == middle_max_vals[i][j-1] - self.gap_initiation_penalty:
            lower_backtrack[i][j] = _prev_node_is_middle_level_
        else:
            ValueError("Unexpected value, lower_max_vals[{0}][{1}] == {2}".format(str(i), str(j), str(lower_max_vals[i][j])))


def print_matrices(max_vals_matrices, backtrack_matrices, nucleotide_h, nucleotide_w):
    for m in max_vals_matrices:
        print()
        print_matrix(m, nucleotide_h, nucleotide_w)
    for m in backtrack_matrices:
        print()
        print_matrix(m, nucleotide_h, nucleotide_w)


def print_matrix(list_of_lists, str_v, str_h):
    height = len(list_of_lists)
    width = len(list_of_lists[0])

    if height > 20 or width > 20:
        print("Matrix will be too large to analyze.")
        print("Truncating to two subgraphs of first 10 and last 10.")

        print("first 10")
        out_str = "           {0}".format("  ".join([ch.rjust(5) for ch in str_h]))
        print(out_str[:79])
        for i in range(10):
            #print(list_of_lists[i])
            #print(str_v[i])
            #print()
            out_str = "  ".join([str(n).rjust(5) for n in list_of_lists[i]])
            out_str = out_str[:79]
            out_str = out_str.replace("↖︎", " ↖︎") #on the terminal this character is incorrectly right-justified by one-too-few spaces, this is not true n other formats
            if i == 0:
                print("    {0}".format(out_str))
            else:
                print("{0}   {1}".format(str_v[i-1], out_str))

        print("last 10")
        out_str = "           {0}".format("  ".join([ch.rjust(5) for ch in str_h]))
        print("    " + out_str[-75:])
        for i in range(height-10, height):
            #print(list_of_lists[i])
            #print(str_v[i])
            #print()
            out_str = "  ".join([str(n).rjust(5) for n in list_of_lists[i]])
            out_str = out_str[-75:]
            out_str = out_str.replace("↖︎", " ↖︎") #on the terminal this character is incorrectly right-justified by one-too-few spaces, this is not true n other formats
            if i == 0:
                print("    {0}".format(out_str))
            else:
                print("{0}   {1}".format(str_v[i-1], out_str))

    #print()
    #print(len(str_v))
    #print(len(list_of_lists))
    else:
        print("           {0}".format("  ".join([ch.rjust(5) for ch in str_h])))
        for i in range(height):
            #print(list_of_lists[i])
            #print(str_v[i])
            #print()
            out_str = "  ".join([str(n).rjust(5) for n in list_of_lists[i]])
            out_str = out_str.replace("↖︎", " ↖︎") #on the terminal this character is incorrectly right-justified by one-too-few spaces, this is not true n other formats
            if i == 0:
                print("    {0}".format(out_str))
            else:
                print("{0}   {1}".format(str_v[i-1], out_str))


def lcsbacktrack(nucleotide_h, nucleotide_w, alignment_strategy):
    backtrack_matrices, max_vals_matrices = alignment_strategy.init_matrixes(nucleotide_h, nucleotide_w)
    lower_max_vals=middle_max_vals=upper_max_vals = []
    if alignment_strategy.three_level_dag:
        lower_max_vals = max_vals_matrices[0]
        middle_max_vals = max_vals_matrices[1]
        upper_max_vals = max_vals_matrices[2]
    else:
        middle_max_vals = max_vals_matrices[0]

    max_overall = -math.inf
    max_overall_loc = []
    all_maxes = []
    for i in range(1, len(nucleotide_h)+1):
        for j in range(1, len(nucleotide_w)+1):
            #print("i:{0} j:{1}".format(str(i), str(j)))
            match = alignment_strategy.scoring[nucleotide_w[j-1]][nucleotide_h[i-1]]

            max_results = alignment_strategy.get_max_val_for_loc(max_vals_matrices,
                                                                 i, j, match )

            alignment_strategy.assign_matrix_vals(max_vals_matrices, backtrack_matrices, max_results, i, j, match)


#            all_maxes, max_overall, max_overall_loc = alignment_strategy.set_max(i, j, max_val_matrices, max_overall, len(nucleotide_h), len(nucleotide_w))

            if i > 0 and j > 0:
                new_max_val = alignment_strategy.found_new_max_value(middle_max_vals[i][j], max_overall, (i, j), len(nucleotide_w), len(nucleotide_h))
                if new_max_val:
                    if middle_max_vals[i][j] > max_overall:
                        all_maxes = []

                    max_overall_loc = [i, j]
                    #print(max_vals[i][j])
                    if alignment_strategy.three_level_dag:
                        max_overall = max(lower_max_vals[i][j], middle_max_vals[i][j], upper_max_vals[i][j])
                        if max_overall == upper_max_vals[i][j]:
                            max_overall_loc.append(2)
                        elif max_overall == middle_max_vals[i][j]:
                            max_overall_loc.append(1)
                        elif max_overall == lower_max_vals[i][j]:
                            max_overall_loc.append(0)
                        else:
                            ValueError("Unexpected max value {0}. Possible values were: {1}, {2}, {3}".format(str(max_overall), str(lower_max_vals[i][j]), str(middle_max_vals[i][j]), str(upper_max_vals[i][j])))
                    else:
                        max_overall = middle_max_vals[i][j]
                    all_maxes.append(max_overall_loc)

        #print()
        #print_matrix(max_vals)
        #print_matrix(backtrack)


    #print()
    #print(all_maxes)
    f_score = middle_max_vals[max_overall_loc[0]][max_overall_loc[1]]

    out_h, out_w = alignment_strategy.outputlcs(backtrack_matrices, max_overall_loc, nucleotide_h, nucleotide_w)


    print_matrices(max_vals_matrices, backtrack_matrices, nucleotide_h, nucleotide_w)

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
        alignment_strategy = AlignmentStrategyGlobal(5, \
                                 True, scoring_filename = "./BLOSUM62.txt")
    elif alignment_type == "local":
        alignment_strategy = AlignmentStrategyLocal(5, \
                                 True, scoring_filename = "./PAM250_scoring.txt")
    elif alignment_type == "fitting":
        alignment_strategy = AlignmentStrategyFitting(1, \
                                 False, scoring_match_value = 1, scoring_mismatch_value = -1)
    elif alignment_type == "overlap":
        alignment_strategy = AlignmentStrategyOverlap(2, \
                                 False, scoring_match_value = 1, scoring_mismatch_value = -2)
    elif alignment_type == "affine_gap":
        alignment_strategy = AlignmentStrategyAffineGap(1, 11, \
                                 True, scoring_filename = "./BLOSUM62.txt")
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
