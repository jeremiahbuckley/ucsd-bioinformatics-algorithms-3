#! /bin/python3

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

def outputlcs(backtrack, v, height, width):

    output = ""
    while height != 0 and width != 0:
        if backtrack[height][width] == -1:
            height -= 1
        elif backtrack[height][width] == 1:
            width -= 1
        else:
            height -= 1
            width -= 1
            output = v[height] + output
    return output

def outputlcs_recursive(backtrack, v, height, width):
    if height == 0 or width == 0:
        return ""
    if backtrack[height][width] == -1:
        return outputlcs(backtrack, v, height-1, width)
    elif backtrack[height][width] == 1:
        return outputlcs(backtrack, v, height, width-1)
    else:
        return outputlcs(backtrack, v, height-1, width-1) + v[height]

def lcsbacktrack(nucleotide_h, nucleotide_w):
    max_vals = []
    backtrack = []
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
        max_vals[i][0] = 0
    for j in range(len(nucleotide_w)+1):
        max_vals[0][j] = 0

    #print()
    #print_matrix(max_vals)
    #print_matrix(backtrack)

    for i in range(1, len(nucleotide_h)+1):
        for j in range(1, len(nucleotide_w)+1):
            print("i:{0} j:{1}".format(str(i), str(j)))
            match = 0
            if nucleotide_h[i-1] == nucleotide_w[j-1]:
                match = 1
            max_vals[i][j] = max(max_vals[i-1][j], max_vals[i][j-1], max_vals[i-1][j-1]+match)

            if max_vals[i][j] == max_vals[i-1][j]:
                backtrack[i][j] = -1
            elif max_vals[i][j] == max_vals[i][j-1]:
                backtrack[i][j] = 1
            else:
                backtrack[i][j] = 0
        #print()
        #print_matrix(max_vals)
        #print_matrix(backtrack)


    return outputlcs(backtrack, nucleotide_h, len(nucleotide_h), len(nucleotide_w))

def read_edges(edges):
    es = []
    for e in edges:
        print(e)
        str1 = e.split("->")
        str2 = str1[1].split(":")
        start = str1[0]
        end = str2[0]
        weight = int(str2[1])
        es.append((start, (end, weight)))
    print(es)
    return es

def dag_backtrack(backtrack, start_node, end_node):
    if end_node == start_node:
        return str(start_node)
    else:
        return dag_backtrack(backtrack, start_node, backtrack[end_node]) + "->" + str(end_node) 

def dag_longest_path(start_node, end_node, edges):
    max_weight = {start_node:0}
    backtrack = {start_node:-1}
    es = read_edges(edges)
    nodes_to_visit = []
    for e in es:
        if e[0] == start_node:
            max_weight[e[1][0]] = e[1][1]
            backtrack[e[1][0]] = start_node
            nodes_to_visit.append(e[1][0])
    print(max_weight)
    print(backtrack)
    print(nodes_to_visit)

    cutoff = 1000
    while len(nodes_to_visit) > 0 and cutoff > 0:
        next_node = nodes_to_visit.pop(0)
        for e in es:
            if e[0] == next_node:
                weight_target = max_weight.get(e[1][0], 0)
                if max_weight[next_node] + e[1][1] > weight_target:
                    max_weight[e[1][0]] = max_weight[next_node] + e[1][1]
                    backtrack[e[1][0]] = next_node
                if e[1][0] not in nodes_to_visit:
                    nodes_to_visit.append(e[1][0])
        nodes_to_visit.sort()
        print()
        print(max_weight)
        print(backtrack)
        print(nodes_to_visit)
        cutoff -= 1
    
    return max_weight[end_node], dag_backtrack(backtrack, start_node, end_node)

if __name__ == '__main__':
    start = time.process_time()

    if len(sys.argv) < 2:
        print("Expected input: 1 param,file. Format '<int_start_node>\n<int_end_node><node_list>\nnode format for node list: {start_node}->{end_node}:{weight}")

    with open(sys.argv[1]) as f:
        start_node = f.readline().rstrip()
        end_node = f.readline().rstrip()
        edges = f.readlines()

    weight, path = dag_longest_path(start_node, end_node, edges)
    print(weight)
    print(path)

    end = time.process_time()
    print("Time: {0}".format(end-start))
