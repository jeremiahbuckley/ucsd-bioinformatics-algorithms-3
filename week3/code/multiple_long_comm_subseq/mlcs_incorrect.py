#! /usr/bin/python3

import sys
import time
import math
import pdb

_verbose_ = False
_debug_ = False

_not_calculated_ = -8

class Scoring:
    def __init__ (self, match_reward, mismatch_penalty, indel_penalty):
        self.match_reward = match_reward
        self.mismatch_penalty = mismatch_penalty
        self.indel_penalty = indel_penalty

    def get_match_val(self, protein_list):
        # match = True
        # if len(protein_list) > 1:
        #     prev_p = protein_list[0]
        #     for i in range(1, len(protein_list)):
        #         next_p = protein_list[i]
        #         match = match and (prev_p == next_p)
        
        # if match:
        #     return self.match_reward
        # return self.indel_penalty
        pmatches = {}
        for i in protein_list:
            if i in pmatches:
                pmatches[i] += 1
            else:
                pmatches[i] = 1
        
        max = 0
        for k,v in pmatches.items():
            if v > max:
                max = v

        if max > 1:
            return max * self.match_reward
        return self.indel_penalty


def walk_backwards(score_keeper, incoming_direction_keeper, nucs, existing_paths, node, path, score, i, j, k, tabs):
    if _debug_:
        print("  " * tabs + node)
    # print(existing_paths)
    if node in existing_paths:
        return existing_paths[node]

    paths1 = paths2 = paths3 = paths4 = paths5 = paths6 = paths7 = paths8 = []
    # print(" " * tabs + node + " " + incoming_direction_keeper[i][j][k])
    for dir in incoming_direction_keeper[i][j][k]:
        if dir == "1":
            paths1 = walk_backwards(score_keeper, incoming_direction_keeper, nucs, existing_paths, str(i-1) + str(j-1) + str(k-1), str(i-1) + str(j-1) + str(k-1) + path, score + score_keeper[i-1][j-1][k-1], i-1, j-1, k-1, tabs+1)
        elif dir == "2":
            paths2 = walk_backwards(score_keeper, incoming_direction_keeper, nucs, existing_paths, str(i-1) + str(j-1) + str(k), str(i-1) + str(j-1) + str(k) + path, score + score_keeper[i-1][j-1][k], i-1, j-1, k, tabs+1)
        elif dir == "3":
            paths3 = walk_backwards(score_keeper, incoming_direction_keeper, nucs, existing_paths, str(i) + str(j-1) + str(k-1), str(i) + str(j-1) + str(k-1) + path, score + score_keeper[i][j-1][k-1], i, j-1, k-1, tabs+1)
        elif dir == "4":
            paths4 = walk_backwards(score_keeper, incoming_direction_keeper, nucs, existing_paths, str(i) + str(j-1) + str(k), str(i) + str(j-1) + str(k) + path, score + score_keeper[i][j-1][k], i, j-1, k, tabs+1)
        elif dir == "5":
            paths5 = walk_backwards(score_keeper, incoming_direction_keeper, nucs, existing_paths, str(i-1) + str(j) + str(k-1), str(i-1) + str(j) + str(k-1) + path, score + score_keeper[i-1][j][k-1], i-1, j, k-1, tabs+1)
        elif dir == "6":
            paths6 = walk_backwards(score_keeper, incoming_direction_keeper, nucs, existing_paths, str(i-1) + str(j) + str(k), str(i-1) + str(j) + str(k) + path, score + score_keeper[i-1][j][k], i-1, j, k, tabs+1)
        elif dir == "7":
            paths7 = walk_backwards(score_keeper, incoming_direction_keeper, nucs, existing_paths, str(i) + str(j) + str(k-1), str(i) + str(j) + str(k-1) + path, score + score_keeper[i][j][k-1], i, j, k-1, tabs+1)
        else:
            # existing_paths[node] = [node]
            return [[(i,j,k)]]

    # print("  " * tabs + ",".join(paths1) + " " ".".join(paths2) + " " + "/".join(paths3) + " " + "|".join(paths4) + " " + "\\".join(paths5) + " " + "`".join(paths6) + " " "'".join(paths7))
    node_paths = []
    p1 = []
    for p in paths1:
        p1.append(p)
    for p in paths2:
        p1.append(p)
    for p in paths3:
        p1.append(p)
    for p in paths4:
        p1.append(p)
    for p in paths5:
        p1.append(p)
    for p in paths6:
        p1.append(p)
    for p in paths7:
        p1.append(p)
    # print(" " * tabs + ": " + ",".join(p1))
    print(p1)
    for paths_collection in p1:
    # for paths_collection in paths1,paths2,paths3,paths4,paths5,paths6,paths7:
        print(paths_collection)
        print('hi')
        path = paths_collection
        print(path)
        path.append((i,j,k))
        print(path)
        if path not in node_paths:
            node_paths.append(path)
        # for pc in paths_collection:
        #     print(" " * tabs + 'pc')
        #     print(" " * tabs + pc + " " + node)
        #     for path in pc:
        #         pathplus = path + node
        #         print("  " * tabs + pathplus)
        #         if pathplus not in node_paths:
        #             node_paths.append(pathplus)

    existing_paths[(i,j,k)] = node_paths
    return node_paths


def print_space(score_keeper, incoming_direction_keeper, nucs):
    str1 = ""
    instr1 = ""
    for i in range(len(nucs[0])+1):
        str2 = ""
        instr2 = ""
        for j in range(len(nucs[1])+1):
            str3 = ", ".join([str(x).rjust(3) for x in score_keeper[i][j]])
            instr3 = ", ".join([x.rjust(5) for x in incoming_direction_keeper[i][j]])
            if len(str2) > 0:
                str2 += "\n"
                instr2 +="\n"
            str2 += str3
            instr2 += instr3
        if len(str1) > 0:
            str1 +="\n\n"
            instr1 +="\n\n"
        str1 += str2
        instr1 += instr2

    print(str1)
    print(instr1)

def find_common_subseq(nucs, scoring):
    score_keeper = []
    incoming_direction_keeper = []
    for i in range(len(nucs[0]) + 1):
        plane = []
        in_dir_plane = []
        for j in range(len(nucs[1])+1):
            edge = []
            in_dir_edge = []
            for k in range(len(nucs[2])+1):
                edge.append(_not_calculated_)
                in_dir_edge.append("*")
            plane.append(edge)
            in_dir_plane.append(in_dir_edge)
        score_keeper.append(plane)
        incoming_direction_keeper.append(in_dir_plane)
    
    # directions example
    # three strings
    # ATCGT
    # ATAGT
    # ATTGT
    # if you ware walking backwards through the matrix, reconstructing the alignments. Let's say you have gotten to this point:
    # "...GT"
    # "...GT"
    # "...GT"
    # and you are trying to figure out the 3rd-from-the-right alignment (first-from-right and second-from-right both lined up perfectly! awesome!)
    # 1 means:
    # "..CGT"
    # "..AGT"
    # "..TGT"
    # 2 means:
    # "..CGT"
    # "..-GT"
    # "..TGT"
    # 3 means:
    # "..-GT"
    # "..AGT"
    # "..TGT"
    # 4 means:
    # "..-GT"
    # "..-GT"
    # "..TGT"
    # 5 means:
    # "..CGT"
    # "..AGT"
    # "..-GT"
    # 6 means:
    # "..CGT"
    # "..-GT"
    # "..-GT"
    # 7 means:
    # "..-GT"
    # "..AGT"
    # "..-GT"

    
    # make this for-real n-dimensional
    for i in range(len(nucs[0]) + 1):
        for j in range(len(nucs[1]) + 1):
            score_keeper[i][j][0] = scoring.indel_penalty
            if i == 0 and j == 0:
                incoming_direction_keeper[i][j][0] = "-"
            elif i > 0 and j == 0:
                incoming_direction_keeper[i][j][0] = "6"
            elif i == 0 and j > 0:
                incoming_direction_keeper[i][j][0] = "4"
            else:
                incoming_direction_keeper[i][j][0] = "246"
    for i in range(len(nucs[1]) + 1):
        for j in range(len(nucs[2]) + 1):
            score_keeper[0][i][j] = scoring.indel_penalty
            if i == 0 and j == 0:
                incoming_direction_keeper[0][i][j] = "-"
            elif i > 0 and j == 0:
                incoming_direction_keeper[0][i][j] = "4"
            elif i == 0 and j > 0:
                incoming_direction_keeper[0][i][j] = "7"
            else:
                incoming_direction_keeper[0][i][j] = "347"
    for i in range(len(nucs[0]) + 1):
        for j in range(len(nucs[2]) + 1):
            score_keeper[i][0][j] = scoring.indel_penalty
            if i == 0 and j == 0:
                incoming_direction_keeper[i][0][j] = "-"
            elif i > 0 and j == 0:
                incoming_direction_keeper[i][0][j] = "6"
            elif i == 0 and j > 0:
                incoming_direction_keeper[i][0][j] = "7"
            else:
                incoming_direction_keeper[i][0][j] = "567"
                
    if _verbose_:
        print("prep")
        print_space(score_keeper, incoming_direction_keeper, nucs)
    

    for i in range(1, len(nucs[0]) + 1):
        plane = []
        for j in range(1, len(nucs[1]) + 1):
            edge = []
            for k in range(1, len(nucs[2]) + 1):
                score = max(score_keeper[i-1][j-1][k - 1] + scoring.get_match_val([nucs[0][i-1], nucs[1][j-1], nucs[2][k-1]]), \
                                score_keeper[i-1][j-1][k] + scoring.indel_penalty, score_keeper[i-1][j][k-1] + scoring.indel_penalty, \
                                score_keeper[i][j-1][k-1] + scoring.indel_penalty, \
                                score_keeper[i-1][j][k] + scoring.indel_penalty, score_keeper[i][j-1][k] + scoring.indel_penalty, score_keeper[i][j][k-1] + scoring.indel_penalty)
                score_keeper[i][j][k] = score
                dir = ""
                if score == score_keeper[i-1][j-1][k-1] + scoring.get_match_val([nucs[0][i-1], nucs[1][j-1], nucs[2][k-1]]):
                    dir += "1"
                if score == score_keeper[i-1][j-1][k] + scoring.indel_penalty:
                    dir += "2"
                if score == score_keeper[i][j-1][k-1] + scoring.indel_penalty:
                    dir += "3"
                if score == score_keeper[i][j-1][k] + scoring.indel_penalty:
                    dir += "4"
                if score == score_keeper[i-1][j][k-1] + scoring.indel_penalty:
                    dir += "5"
                if score == score_keeper[i-1][j][k] + scoring.indel_penalty:
                    dir += "6"
                if score == score_keeper[i][j][k-1] + scoring.indel_penalty:
                    dir += "7"
                incoming_direction_keeper[i][j][k] = dir
    

    if _verbose_:
        print("complete")
        print_space(score_keeper, incoming_direction_keeper, nucs)
    
    existing_paths = {}
    paths = walk_backwards(score_keeper, incoming_direction_keeper, nucs, existing_paths, str(len(nucs[0])) + str(len(nucs[1])) + str(len(nucs[2])), "", score_keeper[len(nucs[0])][len(nucs[1])][len(nucs[2])], len(nucs[0]), len(nucs[1]), len(nucs[2]), 0)
    for i in range(5):
        print()
    
    print("score: " + str(score_keeper[len(nucs[0])][len(nucs[1])][len(nucs[2])]))

    for i in range(5):
        print()
    
    candidates = []
    all_candidates_longest_subseq = 0
    for x in paths:
        if _debug_:
            print(x)
        i = 0
        j = 0
        k = 0
        str1 = str2 = str3 = ""
        for y in range(1, len(x)):
            nuc_idx = x[y][0]
            if nuc_idx == i + 1:
                str1 += nucs[0][i]
            else:
                str1 += "-"
            i = nuc_idx

            nuc_idx = x[y][1]
            if nuc_idx == j + 1:
                str2 += nucs[1][j]
            else:
                str2 += "-"
            j = nuc_idx

            nuc_idx = x[y][2]
            if nuc_idx == k + 1:
                str3 += nucs[2][k]
            else:
                str3 += "-"
            k = nuc_idx
            if x == "00011112123223334444454454565676787898910910101010":
                print("y: {0} i: {1} j: {2} k: {3} str1: {4} str2: {5} str3: {6}".format(str(y), str(i), str(j), str(k), str1, str2, str3))

        largest_common_subseq = 0
        current_subseqence = 0
        for idx in range(len(str1)):
            if str1[idx] == "-" or str2[idx] == "-" or str3[idx] == "-":
                current_subseqence = 0
            else:
                current_subseqence += 1
            if current_subseqence > largest_common_subseq:
                largest_common_subseq = current_subseqence

        if largest_common_subseq > all_candidates_longest_subseq:
            all_candidates_longest_subseq = largest_common_subseq

        candidates.append([largest_common_subseq, str1, str2, str3])

        if _debug_:
            print(str1)
            print(str2)
            print(str3)
            print()

    winners = []
    for candidate in candidates:
        winners.append([candidate[1],candidate[2],candidate[3]])

    # return all_candidates_longest_subseq, winners
    return score_keeper[len(nucs[0])][len(nucs[1])][len(nucs[2])], winners


if __name__ == '__main__':
    start = time.process_time()

    if len(sys.argv) < 1:
        print("Expected input:\n[str: filename path]\n\nfile contents:\n[string: nucleotide]\n[string: nucleotide]\n[string: nucleotide]\n-v = verbose, -vv = debug")

    with open(sys.argv[1]) as f:
        nuc_1 = f.readline().rstrip()
        nuc_2 = f.readline().rstrip()
        nuc_3 = f.readline().rstrip()

    match_reward = 1
    mismatch_penalty = 0
    indel_penalty = 0
    scoring = Scoring(match_reward, mismatch_penalty, indel_penalty)


    for a_idx in range(2, 3, 1):
        if len(sys.argv) > a_idx:
            if sys.argv[a_idx] == "-v":
                _verbose_ = True
            elif sys.argv[a_idx] == "-vv":
                _verbose_ = True
                _debug_ = True

    if _verbose_:
        print(nuc_1)
        print(nuc_2)
        print(nuc_3)
    results = find_common_subseq([nuc_1, nuc_2, nuc_3], scoring)
 

    end = time.process_time()
    print("Time: {0}".format(end-start))
