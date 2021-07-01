#! /usr/bin/python3

import sys
import time

PAM250 = {'A': {'A': 2, 'C': -2, 'E': 0, 'D': 0, 'G': 1, 'F': -3, 'I': -1, 'H': -1, 'K': -1, 'M': -1, 'L': -2, 'N': 0, 'Q': 0, 'P': 1, 'S': 1, 'R': -2, 'T': 1, 'W': -6, 'V': 0, 'Y': -3}, \
'C': {'A': -2, 'C': 12, 'E': -5, 'D': -5, 'G': -3, 'F': -4, 'I': -2, 'H': -3, 'K': -5, 'M': -5, 'L': -6, 'N': -4, 'Q': -5, 'P': -3, 'S': 0, 'R': -4, 'T': -2, 'W': -8, 'V': -2, 'Y': 0}, \
'E': {'A': 0, 'C': -5, 'E': 4, 'D': 3, 'G': 0, 'F': -5, 'I': -2, 'H': 1, 'K': 0, 'M': -2, 'L': -3, 'N': 1, 'Q': 2, 'P': -1, 'S': 0, 'R': -1, 'T': 0, 'W': -7, 'V': -2, 'Y': -4}, \
'D': {'A': 0, 'C': -5, 'E': 3, 'D': 4, 'G': 1, 'F': -6, 'I': -2, 'H': 1, 'K': 0, 'M': -3, 'L': -4, 'N': 2, 'Q': 2, 'P': -1, 'S': 0, 'R': -1, 'T': 0, 'W': -7, 'V': -2, 'Y': -4}, \
'G': {'A': 1, 'C': -3, 'E': 0, 'D': 1, 'G': 5, 'F': -5, 'I': -3, 'H': -2, 'K': -2, 'M': -3, 'L': -4, 'N': 0, 'Q': -1, 'P': 0, 'S': 1, 'R': -3, 'T': 0, 'W': -7, 'V': -1, 'Y': -5}, \
'F': {'A': -3, 'C': -4, 'E': -5, 'D': -6, 'G': -5, 'F': 9, 'I': 1, 'H': -2, 'K': -5, 'M': 0, 'L': 2, 'N': -3, 'Q': -5, 'P': -5, 'S': -3, 'R': -4, 'T': -3, 'W': 0, 'V': -1, 'Y': 7}, \
'I': {'A': -1, 'C': -2, 'E': -2, 'D': -2, 'G': -3, 'F': 1, 'I': 5, 'H': -2, 'K': -2, 'M': 2, 'L': 2, 'N': -2, 'Q': -2, 'P': -2, 'S': -1, 'R': -2, 'T': 0, 'W': -5, 'V': 4, 'Y': -1}, \
'H': {'A': -1, 'C': -3, 'E': 1, 'D': 1, 'G': -2, 'F': -2, 'I': -2, 'H': 6, 'K': 0, 'M': -2, 'L': -2, 'N': 2, 'Q': 3, 'P': 0, 'S': -1, 'R': 2, 'T': -1, 'W': -3, 'V': -2, 'Y': 0}, \
'K': {'A': -1, 'C': -5, 'E': 0, 'D': 0, 'G': -2, 'F': -5, 'I': -2, 'H': 0, 'K': 5, 'M': 0, 'L': -3, 'N': 1, 'Q': 1, 'P': -1, 'S': 0, 'R': 3, 'T': 0, 'W': -3, 'V': -2, 'Y': -4}, \
'M': {'A': -1, 'C': -5, 'E': -2, 'D': -3, 'G': -3, 'F': 0, 'I': 2, 'H': -2, 'K': 0, 'M': 6, 'L': 4, 'N': -2, 'Q': -1, 'P': -2, 'S': -2, 'R': 0, 'T': -1, 'W': -4, 'V': 2, 'Y': -2}, \
'L': {'A': -2, 'C': -6, 'E': -3, 'D': -4, 'G': -4, 'F': 2, 'I': 2, 'H': -2, 'K': -3, 'M': 4, 'L': 6, 'N': -3, 'Q': -2, 'P': -3, 'S': -3, 'R': -3, 'T': -2, 'W': -2, 'V': 2, 'Y': -1}, \
'N': {'A': 0, 'C': -4, 'E': 1, 'D': 2, 'G': 0, 'F': -3, 'I': -2, 'H': 2, 'K': 1, 'M': -2, 'L': -3, 'N': 2, 'Q': 1, 'P': 0, 'S': 1, 'R': 0, 'T': 0, 'W': -4, 'V': -2, 'Y': -2}, \
'Q': {'A': 0, 'C': -5, 'E': 2, 'D': 2, 'G': -1, 'F': -5, 'I': -2, 'H': 3, 'K': 1, 'M': -1, 'L': -2, 'N': 1, 'Q': 4, 'P': 0, 'S': -1, 'R': 1, 'T': -1, 'W': -5, 'V': -2, 'Y': -4}, \
'P': {'A': 1, 'C': -3, 'E': -1, 'D': -1, 'G': 0, 'F': -5, 'I': -2, 'H': 0, 'K': -1, 'M': -2, 'L': -3, 'N': 0, 'Q': 0, 'P': 6, 'S': 1, 'R': 0, 'T': 0, 'W': -6, 'V': -1, 'Y': -5}, \
'S': {'A': 1, 'C': 0, 'E': 0, 'D': 0, 'G': 1, 'F': -3, 'I': -1, 'H': -1, 'K': 0, 'M': -2, 'L': -3, 'N': 1, 'Q': -1, 'P': 1, 'S': 2, 'R': 0, 'T': 1, 'W': -2, 'V': -1, 'Y': -3}, \
'R': {'A': -2, 'C': -4, 'E': -1, 'D': -1, 'G': -3, 'F': -4, 'I': -2, 'H': 2, 'K': 3, 'M': 0, 'L': -3, 'N': 0, 'Q': 1, 'P': 0, 'S': 0, 'R': 6, 'T': -1, 'W': 2, 'V': -2, 'Y': -4}, \
'T': {'A': 1, 'C': -2, 'E': 0, 'D': 0, 'G': 0, 'F': -3, 'I': 0, 'H': -1, 'K': 0, 'M': -1, 'L': -2, 'N': 0, 'Q': -1, 'P': 0, 'S': 1, 'R': -1, 'T': 3, 'W': -5, 'V': 0, 'Y': -3}, \
'W': {'A': -6, 'C': -8, 'E': -7, 'D': -7, 'G': -7, 'F': 0, 'I': -5, 'H': -3, 'K': -3, 'M': -4, 'L': -2, 'N': -4, 'Q': -5, 'P': -6, 'S': -2, 'R': 2, 'T': -5, 'W': 17, 'V': -6, 'Y': 0}, \
'V': {'A': 0, 'C': -2, 'E': -2, 'D': -2, 'G': -1, 'F': -1, 'I': 4, 'H': -2, 'K': -2, 'M': 2, 'L': 2, 'N': -2, 'Q': -2, 'P': -1, 'S': -1, 'R': -2, 'T': 0, 'W': -6, 'V': 4, 'Y': -2}, \
'Y': {'A': -3, 'C': 0, 'E': -4, 'D': -4, 'G': -5, 'F': 7, 'I': -1, 'H': 0, 'K': -4, 'M': -2, 'L': -1, 'N': -2, 'Q': -4, 'P': -5, 'S': -3, 'R': -4, 'T': -3, 'W': 0, 'V': -2, 'Y': 10}}


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

def LocalAlignmentBackTrack(seq1, seq2, aa_dict, *, penalty=-5):
    
    # create backtrack with blank initial values
    backtrack = [['_']*len(seq2) for _ in range(len(seq1))]
    
    # score matrix with considering of free rides (no negative values)
    fr = [[0 for j in range(len(seq2)+1)]] + \
         [[0] + ['_']*len(seq2) for i in range(len(seq1))]
    
    # coordinates of max value in 'free ride' matrix
    max_i, max_j, max_value = 0, 0, 0
        
    # fill in matrices cell by cell
    for i, ch1 in enumerate(seq1, start=1):
        for j, ch2 in enumerate(seq2, start=1):
            fr_i_j = max(fr[i-1][j] + penalty, fr[i][j-1] + penalty, fr[i-1][j-1] + aa_dict[ch1][ch2])
            if fr_i_j == fr[i-1][j] + penalty:    # if move from top
                backtrack[i-1][j-1] = '↓'         
            elif fr_i_j == fr[i][j-1] + penalty:  # if move from left
                backtrack[i-1][j-1] = '→'
            else:                                 # if diagonal movement
                backtrack[i-1][j-1] = '↘'
            if fr_i_j < 0:
                fr_i_j = 0
            if fr_i_j > max_value:
                max_value = fr_i_j
                max_i = i
                max_j = j
            fr[i][j] = fr_i_j

    print_matrix(fr, seq2, seq1)
    return backtrack, fr, max_i-1, max_j-1

def LocalAlignmentOutput(backtrack, score, max_i, max_j, seq1, seq2, aa_dict, *, penalty=-5):
    i = max_i
    j = max_j
    seq1_mod = ''   # modified version of string seq1
    seq2_mod = ''   # modified version of string seq2
    lcs = ''        # calculated longestg common sequence
    sc = 0          # score of calculated longest local sequence
    
    while score[i+1][j+1] != 0:
        ch = backtrack[i][j]
        if ch == '↓':
            seq1_mod += seq1[i]
            seq2_mod += '-'
            sc += penalty
            i -= 1
        elif ch == '→':
            seq1_mod += '-'
            seq2_mod += seq2[j]
            sc += penalty
            j -= 1
        elif ch == '↘':
            seq1_mod += seq1[i]
            seq2_mod += seq2[j]
            if seq1[i] == seq2[j]:
                lcs += seq1[i]
            sc += aa_dict[seq1[i]][seq2[j]]
            i -= 1
            j -= 1
        else:
            raise ValueError('Unrecognize symbol in matrix: ' + ch)
            
    return seq1_mod[::-1], seq2_mod[::-1], lcs[::-1], sc


seq1 = 'CYDFHYTKPNNKRDAYPSIKQQFDQLIAILPHWMYSPAVRHNDDWKVHMDEIWFQHAGSSHMGCRAFDDSPPQGVRSSHQCNVMFDNPTLALTWIAEKMRPHDTVHLCHHVMWMIFSDDFWYYQLCWDSEMKRFLRTNSGGMAGHVDKMKHENMTDPCWCTNFTREIAKKTRLYNAKQTCHSPSWADTSQCEMKCHNDLHKNSGQAMRMCGESMLEMHNDCKKWTINAAGYWNLKTRMDISFPEIMSPMQNLSNNGRTPGRLQTDVVDFPLGVIEDQNNSWFDAHGTMNSNLYGHNQCEIHNYGHHMVPCCQDNVVENVPACYWKSGTQEDNIEKNNNGWFIYRLVKSQFDSDIGNVISCATCQYAPYYRAATEEKQVSDFLYKWARSWCTMMNPWAPLEFCINGIKHNYNELHSGIGGHNETNTAYRYMCWGYQWHFRDHYKFSWSPTITMSVADWHNRYPELCRDRNWHAEDFLATDYFTLMFNWGHGMAWESRAAVHFRDIGPFDFRPMPGDEAEECDCYFAWSWASPCFSERFLWWWYVHGDQERNGYKEIRYVNRHGFCNCMPYGNHTRTMPWGMRFFTYRNRQMNHCSVRNIVIFPSRHFWKSDPVDASPHYGNKCWFHIDWMYWLYGQLEEDSPAVDCGGKRQWMNISSMPLALKGMCLTPLRDRYSYNDDCYVINNKRCNFSENKAHWWWCRYGMFWNDGFVTMDYNVYKHVFDMRIVQQFRHQRHWGVDPFNHLIKTHIVVEQTGHEIGIAVHEIKRQTEAVSKEDSLQNATQCQSMWRHDPEWCDYHVFNGNIMWHIRHPYLGFQYITQWNEILEKKGSNHLDADSPCDNYILVKSQNNPLQHNQNFPVMEHCGRMYECRWAAASMKTLVRSYNLLCAKHERVGERVKYE'
seq2 = 'KDGFHLLAWYEFTTKRSAVGQEQYDVHAYMPYGGHAFMWHHRIRDDFTRWFRLLSQACHNWLEEAFHTWRVFSFMAGWGKRWSMISQLGNSRSGECGEAWCLDCRLSVWMKMSIAMTNFYCFFSGDYRLSMSHICNMMVCMATERPQGNFVHWKWCCMWDFATSWKSWMCIIICFRIDSMAMAQFTRFKQSHFETRPSTCWEGLMDMRCFCKCRCHGPWVIYREFTGRYIGHHHWFNYIHLKGTEWFVTLQIYKDMAWQERHPKFITEDQTMNSNLYGHNQCEINCGHHMVVENDNRRMVWPACRQTPRYVVKSGTQEQGGIKNKICNQDHTKNNVAAIWPPVKEPDQFDSDGGNNISCATCYNAPYYRAVSDFLYCAMMNPWAPLGFCRNGIKHNYNELHYETSARLHGIGGHNETNTAYYQWHFLRDEHWADHYQFSWSPTITMSVADWHNRYCVELNRRNWHAEDFLATDYFTLMAYMINRNWGYGMWWESRYIYYIPFGDEAEECTNRLKRSRQCYFTSWTFSDGVWSDDFGASFSERFLWWWYVHTDFWQAEGNPQVQPLLNGYKAIRFMHFGLQWAVFHCHFRHGFCNCMPYGRQTLVMPPWGMRFFTYRKDRVMNHIASVRNIVIFPSRAWKSDPITLKGMYDPASPHYGNKTWFHCDQMYWPYGWCNPNRCCMNGFWTDAWYKDLTYPGAFETNKGMEGSGKMKHRDCILGRFWFTEQNTLYAGCPIQKHYVSTCHEPIFTAKTLPVICKTDDMPEQESHRMHDLPVSENAEEESTEAYDHSPWFLIKCECIRCKMFRYQQMALEWWGHMMEEATQPKFWNPVHTTSTEAMRPCSDSGGQRRLYQNLNAYVYHSSYLSKSRGRCHDVSVEIIIPWQGNTDTPTWHFMWINAALGHEHTCCYKMWCRNRCDVFDFTCHARQPLRHDFC'

seq1 = 'ATGVWYYYYC'
seq2 = 'GYYCFFF'

seq1 = 'AAAAAASSSSSSSVVVVVVVVTTTTTTTLLLLLL'
seq2 = 'SSSSLLLLLTTTTTTTWWWWWWWWWWWPPPPPPPPPP'

seq1 = 'QTVHQIWMKRLASFFFSMMMRRLLQQSSSTTTFQQDWLLLLLSSSCCSPQPQPQTYTYTYTFFRRRSSMMMLNNNDNDKVWSKLLWPPPPQRMS'
seq2 = 'TQDFFSMLRRAAQFFMMRLCCDPPPATATATFPFPMMMDSNSSSDNSQRRQRQRQTFTDTSTCCQVVMNFRMNFRVDFPLLK'
aa_dict = PAM250

backtrack, score, max_i, max_j = LocalAlignmentBackTrack(seq1, seq2, aa_dict)
# print('Max value of local subpath:', max_i, max_j, score[max_i+1][max_j+1])

seq1_mod, seq2_mod, lcs, sc  = LocalAlignmentOutput(backtrack, score, max_i, max_j, seq1, seq2, aa_dict)

print(sc)        #> 15
print(seq1_mod)  #> EANL-Y
print(seq2_mod)  #> ENALTY
print(lcs)       #> ELY (redundant for given task)