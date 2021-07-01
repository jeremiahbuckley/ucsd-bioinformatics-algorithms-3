#! /usr/bin/python3

import sys
import math

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



def FittingAlignmentBackTrack(seq1, seq2, \
                              *, reward=1, penalty_subst=-1, penalty_indel=-1):
    
    # create backtrack with blank initial values
    backtrack = [['_']*len(seq2) for _ in range(len(seq1))]
    
    # Create score matrix.
    # In the first row I have accumulated penalties for insertions,
    # but in the first column zeros since consider deletions from the ends as "free rides".
    s = [[j*penalty_indel for j in range(len(seq2)+1)]] + \
        [[0] + ['_']*len(seq2) for i in range(len(seq1))]
    
    # max value in the last column and index of row where it is achieved
    max_value, max_i = -math.inf, 0
    
    for i, ch1 in enumerate(seq1, start=1):
        for j, ch2 in enumerate(seq2, start=1):
            match = reward if ch1 == ch2 else penalty_subst
            s[i][j] = max(s[i-1][j] + penalty_indel,
                          s[i][j-1] + penalty_indel,
                          s[i-1][j-1] + match)
            if s[i][j] == s[i-1][j] + penalty_indel:    # if move from top
                backtrack[i-1][j-1] = '↓'         
            elif s[i][j] == s[i][j-1] + penalty_indel:  # if move from left
                backtrack[i-1][j-1] = '→'
            else:                                       # if diagonal movement
                backtrack[i-1][j-1] = '↘'
                
        if s[i][j] > max_value:   # if we found value bigger then maximum before
            max_value = s[i][j]
            max_i = i
                
    print_matrix(s, seq2, seq1)        
    return backtrack, s, max_i-1


def FittingAlignmentOutput(backtrack, max_i, seq1, seq2,
                           *, reward=1, penalty_subst=-1, penalty_indel=-1):
    i = max_i
    j = len(seq2) - 1
    seq1_mod = ''      # modified version of string seq1
    seq2_mod = ''      # modified version of string seq2
    lcs = ''           # calculated longest common sequence
    sc = 0             # score of calculated longest common sequence
    
    while j != -1:
        ch = backtrack[i][j]
        if ch == '↓':
            seq1_mod += seq1[i]
            seq2_mod += '-'
            sc += penalty_indel
            i -= 1
        elif ch == '→':
            seq1_mod += '-'
            seq2_mod += seq2[j]
            sc += penalty_indel
            j -= 1
        elif ch == '↘':
            seq1_mod += seq1[i]
            seq2_mod += seq2[j]
            if seq1[i] == seq2[j]:
                lcs += seq1[i]
                sc += reward
            else:
                sc += penalty_subst
            i -= 1
            j -= 1
        else:
            raise ValueError('Unrecognize symbol in matrix: ' + ch)

    return seq1_mod[::-1], seq2_mod[::-1], lcs[::-1], sc


seq1 = 'GTAGGCTTAAGGTTA'
seq2 = 'TAGATA'
seq1 = 'ACGACAGAG'
seq2 = 'CGAGAGGTT'

backtrack, _, max_i = FittingAlignmentBackTrack(seq1, seq2)
seq1_mod, seq2_mod, lcs, sc  = FittingAlignmentOutput(backtrack, max_i, seq1, seq2)

print(sc)        #> 2
print(seq1_mod)  #> TAGGCTTA
print(seq2_mod)  #> TAGA-T-A
print(lcs)       #> TAGTA (redundant for given task)