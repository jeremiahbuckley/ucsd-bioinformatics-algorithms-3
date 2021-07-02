#! /usr/bin/python3

import sys
import time


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

def manhattan_tourist(height, width, down_weights, right_weights):
    #print(height)
    #print(width)
    #print(down_weights)
    #print(right_weights)

    max_vals = []
    for i in range(height+1):
        row = []
        for j in range(width+1):
            row.append(0)
        max_vals.append(row)
    #print_matrix(max_vals)

    for i in range(height):
        max_vals[i+1][0] = max_vals[i][0] + down_weights[i][0]    
    #print_matrix(max_vals)

    for i in range(width):
        max_vals[0][i+1] = max_vals[0][i] + right_weights[0][i]
    #print_matrix(max_vals)

    for k in range(1, width+1):
        #print(k)
        for i in range(1, height+1):
            for j in range(k, k+1):
                from_left = max_vals[i][j-1] + right_weights[i][j-1]
                from_above = max_vals[i-1][j] + down_weights[i-1][j]

                if from_left > from_above:
                    max_vals[i][j] = from_left
                else:
                    max_vals[i][j] = from_above

            #print_matrix(max_vals)
    return max_vals[height][width]

def outputlcs(backtrack, v, height, width, nucleotide_h, nucleotide_w):

    output = ""
    align_h = ""
    align_w = ""
    while height != 0 and width != 0:
        if backtrack[height][width] == -1:
            align_h = nucleotide_h[height-1] + align_h
            align_w = "-" + align_w
            height -= 1
        elif backtrack[height][width] == 1:
            align_h = "-" + align_h
            align_w = nucleotide_w[width-1] + align_w
            width -= 1
        else:
            align_h = nucleotide_h[height-1] + align_h
            align_w = nucleotide_w[width-1] + align_w
            height -= 1
            width -= 1
            output = v[height] + output
    if height == 0 and width != 0:
        align_h = "-"*(width) + align_h
        align_w = nucleotide_w[0:width] + align_w
    if width == 0 and height != 0:
        align_h = nucleotide_h[0:height] + align_h
        align_w = "-"*(height) + align_w

    #print()
    #print(align_h)
    #print(align_w)
    return align_h, align_w

def outputlcs_recursive(backtrack, v, height, width):
    if height == 0 or width == 0:
        return ""
    if backtrack[height][width] == -1:
        return outputlcs(backtrack, v, height-1, width)
    elif backtrack[height][width] == 1:
        return outputlcs(backtrack, v, height, width-1)
    else:
        return outputlcs(backtrack, v, height-1, width-1) + v[height]

def lcsbacktrack(nucleotide_h, nucleotide_w, scoring):
    indel_penalty = 5
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

    #print()
    #print_matrix(max_vals)
    #print_matrix(backtrack)

    for i in range(len(nucleotide_h)+1):
        max_vals[i][0] = 0 - (indel_penalty*(i))
    for j in range(len(nucleotide_w)+1):
        max_vals[0][j] = 0 - (indel_penalty*(j))

    #print()
    #print_matrix(max_vals)
    #print_matrix(backtrack)

    for i in range(1, len(nucleotide_h)+1):
        for j in range(1, len(nucleotide_w)+1):
            #print("i:{0} j:{1}".format(str(i), str(j)))
            match = scoring[nucleotide_w[j-1]][nucleotide_h[i-1]]
            max_vals[i][j] = max(max_vals[i-1][j] - indel_penalty, \
                                 max_vals[i][j-1] - indel_penalty, \
                                 max_vals[i-1][j-1]+match)

            if max_vals[i][j] == max_vals[i-1][j] - indel_penalty:
                backtrack[i][j] = -1
                alignment_h += nucleotide_h[i-1]
                alignment_w += "-"
                score -= indel_penalty
            elif max_vals[i][j] == max_vals[i][j-1] - indel_penalty:
                backtrack[i][j] = 1
                alignment_h += "-"
                alignment_w += nucleotide_w[j-1]
                score -= indel_penalty
            else:
                backtrack[i][j] = 0
                score += match
                alignment_h += nucleotide_h[i-1]
                alignment_w += nucleotide_w[j-1]
        #print()
        #print(score)
        #print(alignment_h)
        #print(alignment_w)
        #print_matrix(max_vals)
        #print_matrix(backtrack)
    print_matrix(max_vals)
    print_matrix(backtrack)

    #print(score)
    #print(alignment_h)
    #print(alignment_w)
    print(max_vals[len(nucleotide_h)][len(nucleotide_w)])
    f_score = max_vals[len(nucleotide_h)][len(nucleotide_w)]

    out_h, out_w = outputlcs(backtrack, nucleotide_h, len(nucleotide_h), len(nucleotide_w), nucleotide_h, nucleotide_w)
    return f_score, out_h, out_w

if __name__ == '__main__':
    start = time.process_time()

    if len(sys.argv) < 2:
        print("Expected input: 1 param,file. Format '<str_nucleotide_1>\n<str_nucleotide_2>")

    with open(sys.argv[1]) as f:
        nucleotide_h = f.readline().rstrip()
        nucleotide_w = f.readline().rstrip()

    scoring_strs = []
    with open("./BLOSUM62.txt") as f:
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

    f_score, out_h, out_w = lcsbacktrack(nucleotide_h, nucleotide_w, scoring)
    print(f_score)
    print(out_h)
    print(out_w)

    end = time.process_time()
    print("Time: {0}".format(end-start))
