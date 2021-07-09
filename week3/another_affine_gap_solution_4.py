#! /usr/bin/python3
"""
1.1 Penalizing Insertions and Deletions in Sequence Alignment
Code Challenge: Solve the Alignment with Affine Gap Penalties Problem.

Input: Two amino acid strings v and w (each of length at most 100).
Output: The maximum alignment score between v and w, followed by an alignment of v and w achieving this maximum score. Use the BLOSUM62 scoring matrix, a gap opening penalty of 11, and a gap extension penalty of 1.
"""


import numpy as np
import copy
import math


def score_from_file(file):
    """
    build a map of pair scores from file
    """
    with open("blosum62.txt") as file:
        lines = file.read().splitlines()

        chars = lines[0].split()

        score_key = {}

        for l in lines[1:]:
            linevals = l.split()
            for i, s in enumerate(linevals[1:]):
                pair = sorted([linevals[0], chars[i]])
                score_key[":".join(pair)] = int(s)

        return score_key


def path_from_pair(v, w, gap_open_s, gap_ext_e, score_matrix):
    """
    given strings v and w, construct a path matrix
    """
    # use y x coords for score (nodes), use i j coords for backtrack (paths)

    low = 0
    mid = 1
    upp = 2

    score = np.full((3, len(v) + 1, len(w) + 1), -math.inf)
    backtrack = np.full((3, len(v), len(w)), " ")

    # middle edges
    score[mid, 0, 0] = 0
    score[mid, 1, 0] = gap_open_s
    score[mid, 0, 1] = gap_open_s
    for y in range(2, len(v) + 1):
        score[mid, y, 0] = score[mid, y - 1, 0] + gap_ext_e
    for x in range(2, len(w) + 1):
        score[mid, 0, x] = score[mid, 0, x - 1] + gap_ext_e

    score[low, :, 0] = copy.copy(score[mid, :, 0])
    score[upp, 0, :] = copy.copy(score[mid, 0, :])

    for i in range(len(v)):
        for j in range(len(w)):
            y, x = i + 1, j + 1
            pair = [v[i], w[j]]
            pair = ":".join(sorted(pair))

            # lower
            south = score[low, y - 1, x] + gap_ext_e
            down = score[mid, y - 1, x] + gap_open_s
            choices = [south, down]
            paths = ["↓", "v"]
            score[low, y, x] = max(choices)
            backtrack[low, i, j] = paths[choices.index(score[low, y, x])]

            # upper
            east = score[upp, y, x - 1] + gap_ext_e
            up = score[mid, y, x - 1] + gap_open_s
            choices = [east, up]
            paths = ["→", "^"]
            score[upp, y, x] = max(choices)
            backtrack[upp, i, j] = paths[choices.index(score[upp, y, x])]

            # middle
            up = score[low, y, x]
            diag = score[mid, y - 1, x - 1] + score_matrix[pair]
            down = score[upp, y, x]
            choices = [up, diag, down]
            paths = ["u", "↘", "d"]
            score[mid, y, x] = max(choices)
            backtrack[mid, i, j] = paths[choices.index(score[mid, y, x])]

    print("score\n", score)
    print("backtrack\n", backtrack)
    return score, backtrack


def back_path(backtrack, v, w, i, j):
    """
    output path created by backtrack
    """
    vv = list(v)
    ww = list(w)

    istring = []
    jstring = []

    low = 0
    mid = 1
    upp = 2

    lvl = mid

    path = np.full((3, len(v), len(w)), " ")

    while i > -1 and j > -1:
        path[lvl, i, j] = "x"
        b = backtrack[lvl, i, j]
        print(["low", "mid", "upp"][lvl], i, j, b)
        if b == "↓":
            istring = [vv[i]] + istring
            jstring = ["-"] + jstring
            i -= 1
            continue
        if b == "→":
            istring = ["-"] + istring
            jstring = [ww[j]] + jstring
            j -= 1
            continue
        if b == "↘":
            istring = [vv[i]] + istring
            jstring = [ww[j]] + jstring
            i -= 1
            j -= 1
            continue
        if b == "u":
            lvl = low
            continue
        if b == "d":
            lvl = upp
            continue
        if b == "^":
            istring = ["-"] + istring
            jstring = [ww[j]] + jstring
            j -= 1
            lvl = mid
            continue
        if b == "v":
            istring = [vv[i]] + istring
            jstring = ["-"] + jstring
            i -= 1
            lvl = mid
            continue

        raise "error A"

    print("path\n", path)
    print("v", istring)
    print("w", jstring)
    return ["".join(istring), "".join(jstring)]


if __name__ == "__main__":

    score_matrix = score_from_file("blosum62.txt")
    print("score_matrix\n", score_matrix)

    v = "PRTEINS"
    w = "PRTWPSEIN"

    v = "INGQTFAVMPRAQIAIHCYWTFAFNYQVCQYATYFLMDTYDCAEERMECPEFLQAMCYYSNTITNWPMKNEG"
    w = "INGQIFAVMPWTFAFNHEQLPMIYTKYFFAACIRCVTPEKNQAWEEVWYPRWVAKQAMCYNTIHNWPMKNEG"

    # v = "YHFDVPDCWAHRYWVENPQAIAQMEQICFNWFPSMMMKQPHVFKVDHHMSCRWLPIRGKKCSSCCTRMRVRTVWE"
    # w = "YHEDVAHEDAIAQMVNTFGFVWQICLNQFPSMMMKIYWIAVLSAHVADRKTWSKHMSCRWLPIISATCARMRVRTVWE"

    indel = -5
    node_scores, path_matrix = path_from_pair(v, w, gap_open_s=-11, gap_ext_e=-1, score_matrix=score_matrix)

    out = back_path(path_matrix, v, w, len(v) - 1, len(w) - 1)

    mid = 1
    print(int(node_scores[mid, -1, -1]))
    for o in out:
        print("".join(o))