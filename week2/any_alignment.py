#! /usr/bin/python3

import sys
import time

_prev_node_is_up_ = 1
_prev_node_is_left_ = 2
_prev_node_is_diagonal_ = 3
_prev_node_is_zero_ = 4
_prev_node_is_end_local_ = 5

def print_matrix(list_of_lists, str_h, str_v):
    height = len(list_of_lists)
    width = len(list_of_lists[0])

    #print()
    #print(len(str_v))
    #print(len(list_of_lists))
    print("           {0}".format("  ".join([ch.rjust(5) for ch in str_h])))
    for i in range(height):
        #print(list_of_lists[i])
        #print(str_v[i])
        #print()
        if i == 0:
            print("    {0}".format("  ".join([str(n).rjust(5) for n in list_of_lists[i]])))
        else:
            print("{0}   {1}".format(str_v[i-1], "  ".join([str(n).rjust(5) for n in list_of_lists[i]])))

def outputlcs(backtrack, v, height, width, nucleotide_h, nucleotide_w, max_overall_loc, alignment_type):

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

    if alignment_type == "global":
        if height == 0 and width != 0:
            align_h = "-"*(width) + align_h
            align_w = nucleotide_w[0:width] + align_w
        if width == 0 and height != 0:
            align_h = nucleotide_h[0:height] + align_h
            align_w = "-"*(height) + align_w

    #print()
    print(align_h)
    print(align_w)
    return align_h, align_w

def lcsbacktrack(nucleotide_h, nucleotide_w, scoring, alignment_type):
    max_vals = []
    backtrack = []

    indel_penalty = 0
    if alignment_type == "global":
        indel_penalty = 5
    elif alignment_type == "local":
        indel_penalty = 5
    elif alignment_type == "fitting":
        indel_penalty = 1
    elif alignment_type == "overlap":
        indel_penalty = 2
    else:
        raise ValueError("Unknown alingment_type {0}".format(alignment_type))

    for i in range(len(nucleotide_h)+1):
        vals = []
        bv = []
        for j in range(len(nucleotide_w)+1):
            vals.append(100)
            bv.append(100)
        max_vals.append(vals)
        backtrack.append(bv)

    #print()
    #print_matrix(max_vals)
    #print_matrix(backtrack)

    for i in range(len(nucleotide_h)+1):
        init_val = 0
        # TODO: should this if include "overlap"?
        if alignment_type == "global":
            init_val -= indel_penalty*i
        max_vals[i][0] = init_val
        backtrack[i][0] = _prev_node_is_zero_
    for j in range(len(nucleotide_w)+1):
        init_val = 0
        if alignment_type == "global" or alignment_type == "fitting" or alignment_type == "overlap":
            init_val -= indel_penalty*j
        max_vals[0][j] = init_val
        backtrack[0][j] = _prev_node_is_zero_

    #print()
    #print_matrix(max_vals)
    #print_matrix(backtrack)

    max_overall = 0
    max_overall_loc = (0, 0)
    all_maxes = []
    for i in range(1, len(nucleotide_h)+1):
        for j in range(1, len(nucleotide_w)+1):
            #print("i:{0} j:{1}".format(str(i), str(j)))
            # TODO: change fitting, overlap to a scoring matrix for readability
            match = -100
            if alignment_type == "global" or alignment_type == "local":
                match = scoring[nucleotide_w[j-1]][nucleotide_h[i-1]]
            elif alignment_type == "fitting":
                match = -1
                if nucleotide_h[i-1] == nucleotide_w[j-1]:
                    match = 1
            elif alignment_type == "overlap":
                match = -2
                if nucleotide_h[i-1] == nucleotide_w[j-1]:
                    match = 1
            else:
                raise ValueError("Unknown alingment_type {0}".format(alignment_type))

            max_vals[i][j] = max(max_vals[i-1][j] - indel_penalty, \
                                 max_vals[i][j-1] - indel_penalty, \
                                 max_vals[i-1][j-1]+match)
            
            if alignment_type == "local":
                max_vals[i][j] = max(0, max_vals[i][j])
            
            #if (alignment_type == "fitting" or alignment_type == "overlap") and j == 1 and nucleotide_h[i-1] != nucleotide_w[j-1]:
            if alignment_type == "fitting" and j == 1 and nucleotide_h[i-1] != nucleotide_w[j-1]:
                max_vals[i][j] = 0

            if max_vals[i][j] == max_vals[i-1][j] - indel_penalty:
                backtrack[i][j] = _prev_node_is_up_
            elif max_vals[i][j] == max_vals[i][j-1] - indel_penalty:
                backtrack[i][j] = _prev_node_is_left_
            elif max_vals[i][j] == max_vals[i-1][j-1] + match:
                backtrack[i][j] = _prev_node_is_diagonal_
            else: # max_vals[i][j] == 0:
                backtrack[i][j] = _prev_node_is_zero_

            if max_vals[i][j] > max_overall:
                all_maxes = []

            new_max_val = max_vals[i][j] >= max_overall
            if alignment_type == "fitting":
                new_max_val = new_max_val and j == len(nucleotide_w)
            if alignment_type == "overlap":
                new_max_val = new_max_val and i == len(nucleotide_h)
            if new_max_val:
                #print(max_vals[i][j])
                max_overall = max_vals[i][j]
                max_overall_loc = (i, j)
                all_maxes.append(max_overall_loc)
                # still have to figure out backtrack

    if max_overall_loc[0] != len(nucleotide_h) or max_overall_loc[1] != len(nucleotide_w):
        backtrack[len(nucleotide_h)][len(nucleotide_w)] = _prev_node_is_end_local_
        #print()
        #print_matrix(max_vals)
        #print_matrix(backtrack)
    

    #print()
    print_matrix(max_vals, nucleotide_w, nucleotide_h)
    print_matrix(backtrack, nucleotide_w, nucleotide_h)
    #print(all_maxes)
    print(max_vals[len(nucleotide_h)][len(nucleotide_w)])
    f_score = max_vals[max_overall_loc[0]][max_overall_loc[1]]
    out_h, out_w = outputlcs(backtrack, nucleotide_h, max_overall_loc[0], max_overall_loc[1], nucleotide_h, nucleotide_w, max_overall_loc, alignment_type)

    return f_score, out_h, out_w

if __name__ == '__main__':
    start = time.process_time()

    if len(sys.argv) < 3:
        print("Expected input:\n 1 type of alignment <global,local,fitting,overlap>,\n 2 param,file. Format '<str_nucleotide_1>\n<str_nucleotide_2>")

    alignment_type = sys.argv[1]

    scoring_strs = []
    with open("./PAM250_scoring.txt") as f:
        scoring_strs = f.readlines()

    if alignment_type == "global":
        with open("./BLOSUM62.txt") as f:
            scoring_strs = f.readlines()
    elif alignment_type == "local":
        with open("./PAM250_scoring.txt") as f:
            scoring_strs = f.readlines()
    elif alignment_type == "fitting":
        pass
    elif alignment_type == "overlap":
        pass
    else:
        raise ValueError("Unexpected alignment type: {0}".format(alignment_type))

    scoring = {}
    virt_keys = 'ACDEFGHIKLMNPQRSTVWY'

    if alignment_type == "global" or alignment_type == "local":
        for i in range(1, len(scoring_strs)):
            row = scoring_strs[i].split()
            key = row[0]
            row_dict = {}
            for j in range(1, len(row)):
                virt_key = virt_keys[j-1]
                row_dict[virt_key] = int(row[j])
            scoring[key] = row_dict
    elif alignment_type == "fitting" or alignment_type == "overlap":
        mismatch = -1
        if alignment_type == "overlap":
            mismatch = -2
        match = 1

        for i in range(len(virt_keys)):
            key = virt_keys[i]
            row_dict = {}
            for j in range(len(virt_keys)):
                virt_key = virt_keys[j]
                row_dict[virt_key] = match if i == j else mismatch
            scoring[key] = row_dict
    else:
        raise ValueError("Unexpected alignment type: {0}".format(alignment_type))
    
    print(scoring)


    with open(sys.argv[2]) as f:
        nucleotide_h = f.readline().rstrip()
        nucleotide_w = f.readline().rstrip()

    results = lcsbacktrack(nucleotide_h, nucleotide_w, scoring, alignment_type)
    for i in range(3):
        print(results[i])

    end = time.process_time()
    print("Time: {0}".format(end-start))
