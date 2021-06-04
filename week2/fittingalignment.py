#! /bin/python3

import sys
import time

_prev_node_is_up_ = 1
_prev_node_is_left_ = 2
_prev_node_is_diagonal_ = 3
_prev_node_is_zero_ = 4
_prev_node_is_end_local_ = 5

def build_min_coins_list(money, coins):
    min_coins_list = [0]

    for m in range(1, money+1):
        min_coins_list.append(sys.maxsize)
        for c in coins[::-1]:
            if m >= c:
                if min_coins_list[m-c] + 1 < min_coins_list[m]:
                    min_coins_list[m] = min_coins_list[m -c] + 1

    return min_coins_list

def dpchange(money, coins):
    min_coins_list = build_min_coins_list(money, coins)
    return min_coins_list[money]


def print_matrix(list_of_lists):
    height = len(list_of_lists)
    width = len(list_of_lists[0])

    print()
    for i in range(height):
        print("  ".join([str(n).rjust(3) for n in list_of_lists[i]]))


def print_matrix_enhanced(list_of_lists, nuc_h, nuc_w):
    height = len(list_of_lists)
    width = len(list_of_lists[0])

    print()
    print("        {0}".format("  ".join(c.rjust(3) for c in nuc_w)))
    for i in range(height):
        prefix = " "
        if i > 0:
            prefix = nuc_h[i-1]
        print("{0}  {1}".format(prefix, "  ".join([str(n).rjust(3) for n in list_of_lists[i]])))

def manhattan_tourist(height, width, down_weights, right_weights):
    print(height)
    print(width)
    print(down_weights)
    print(right_weights)

    max_vals = []
    for i in range(height+1):
        row = []
        for j in range(width+1):
            row.append(0)
        max_vals.append(row)
    print_matrix(max_vals)

    for i in range(height):
        max_vals[i+1][0] = max_vals[i][0] + down_weights[i][0]    
    print_matrix(max_vals)

    for i in range(width):
        max_vals[0][i+1] = max_vals[0][i] + right_weights[0][i]
    print_matrix(max_vals)

    for k in range(1, width+1):
        print(k)
        for i in range(1, height+1):
            for j in range(k, k+1):
                from_left = max_vals[i][j-1] + right_weights[i][j-1]
                from_above = max_vals[i-1][j] + down_weights[i-1][j]

                if from_left > from_above:
                    max_vals[i][j] = from_left
                else:
                    max_vals[i][j] = from_above

            print_matrix(max_vals)
    return max_vals[height][width]

def outputlcs(backtrack, v, height, width, nucleotide_h, nucleotide_w, max_overall_loc):

    output = ""
    align_h = ""
    align_w = ""
    #print(max_overall_loc)
    #print(nucleotide_h)
    #print(nucleotide_w)
    while height != 0 and width != 0:
        if backtrack[height][width] == _prev_node_is_zero_:
            height = 0
            width = 0
        elif backtrack[height][width] == _prev_node_is_end_local_:
            align_h = nucleotide_h[max_overall_loc[0]-1] + align_h
            align_w = nucleotide_w[max_overall_loc[1]-1] + align_w
            height = max_overall_loc[0]
            width = max_overall_loc[1]
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
            output = v[height] + output

    #if height == 0 and width != 0:
    #    align_h = "-"*(width) + align_h
    #    align_w = nucleotide_w[0:width] + align_w
    #if width == 0 and height != 0:
    #    align_h = nucleotide_h[0:height] + align_h
    #    align_w = "-"*(height) + align_w

    #print()
    print(align_h)
    print(align_w)
    return output

def outputlcs_recursive(backtrack, v, height, width):
    if height == 0 or width == 0:
        return ""
    if backtrack[height][width] == _prev_node_is_up_:
        return outputlcs(backtrack, v, height-1, width)
    elif backtrack[height][width] == _prev_node_is_left_:
        return outputlcs(backtrack, v, height, width-1)
    elif backtrack[height][width] == _prev_node_is_diagonal_:
        return outputlcs(backtrack, v, height-1, width-1) + v[height]
    else: #_prev_node_is_zero_
        return v[height]

def lcsbacktrack(nucleotide_h, nucleotide_w, scoring):
    indel_penalty = 1
    max_vals = []
    backtrack = []
    alignment_h = ""
    alignment_w = ""
    score = 0
    for i in range(len(nucleotide_h)+1):
        vals = []
        bv = []
        for j in range(len(nucleotide_w)+1):
            vals.append(100)
            bv.append(100)
        max_vals.append(vals)
        backtrack.append(bv)

    print()
    print_matrix(max_vals)
    print_matrix(backtrack)

    for i in range(len(nucleotide_h)+1):
        max_vals[i][0] = 0# - (indel_penalty*(i))
        backtrack[i][0] = _prev_node_is_zero_
    for j in range(len(nucleotide_w)+1):
        max_vals[0][j] = 0 - (indel_penalty*(j))
        backtrack[0][j] = _prev_node_is_left_

    print()
    print_matrix(max_vals)
    print_matrix(backtrack)

    max_overall = 0
    max_overall_loc = (0, 0)
    for i in range(1, len(nucleotide_h)+1):
        for j in range(1, len(nucleotide_w)+1):
            #print("i:{0} j:{1}".format(str(i), str(j)))
            match = -1
            if nucleotide_h[i-1] == nucleotide_w[j-1]:
                match = 1
            
            next_val = max(max_vals[i-1][j] - indel_penalty, \
                           max_vals[i][j-1] - indel_penalty, \
                           max_vals[i-1][j-1]+match)
            if j == 1:
                next_val = max(0, next_val)
            max_vals[i][j] = next_val

            if max_vals[i][j] == 0: # max_vals[i][j] == 0
                backtrack[i][j] = _prev_node_is_zero_
                score = 0
                alignment_h = ""
                alignment_w = ""
            elif max_vals[i][j] == max_vals[i-1][j] - indel_penalty:
                backtrack[i][j] = _prev_node_is_up_
                alignment_h += nucleotide_h[i-1]
                alignment_w += "-"
                score -= indel_penalty
            elif max_vals[i][j] == max_vals[i][j-1] - indel_penalty:
                backtrack[i][j] = _prev_node_is_left_
                alignment_h += "-"
                alignment_w += nucleotide_w[j-1]
                score -= indel_penalty
            else: # max_vals[i][j] == max_vals[i-1][j-1] + match:
                backtrack[i][j] = _prev_node_is_diagonal_
                score += match
                alignment_h += nucleotide_h[i-1]
                alignment_w += nucleotide_w[j-1]

            if max_vals[i][j] >= max_overall and j == len(nucleotide_w):
                print(max_vals[i][j])
                max_overall = max_vals[i][j]
                max_overall_loc = (i, j)
                # still have to figure out backtrack

    if max_overall_loc[0] != len(nucleotide_h) or max_overall_loc[1] != len(nucleotide_w):
        backtrack[len(nucleotide_h)][len(nucleotide_w)] = _prev_node_is_end_local_
        #print()
        #print(score)
        #print(alignment_h)
        #print(alignment_w)
        #print_matrix(max_vals)
        #print_matrix(backtrack)
    

    #print(score)
    #print(alignment_h)
    #print(alignment_w)
    print()
    print_matrix_enhanced(max_vals, nucleotide_h, nucleotide_w)
    print_matrix_enhanced(backtrack, nucleotide_h, nucleotide_w)
    #print(max_vals[len(nucleotide_h)][len(nucleotide_w)])
    print(nucleotide_h)
    print(nucleotide_w)
    print(max_vals[max_overall_loc[0]][max_overall_loc[1]])

    return outputlcs(backtrack, nucleotide_h, max_overall_loc[0], max_overall_loc[1], nucleotide_h, nucleotide_w, max_overall_loc)

if __name__ == '__main__':
    start = time.process_time()

    if len(sys.argv) < 2:
        print("Expected input: 1 param,file. Format '<str_nucleotide_1>\n<str_nucleotide_2>")

    with open(sys.argv[1]) as f:
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

    match = lcsbacktrack(nucleotide_h, nucleotide_w, scoring)
    print(match)

    end = time.process_time()
    print("Time: {0}".format(end-start))
