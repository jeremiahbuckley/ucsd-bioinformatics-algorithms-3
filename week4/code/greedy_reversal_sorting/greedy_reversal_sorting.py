#! /usr/bin/python3

import sys
import time
import math

_verbose_ = False
_timed_output_ = False
_debug_ = False

def print_perumutation(working_set):
    pstr = ""
    for w in working_set:
        if w > 0:
            pstr += " +" + str(w)
        else:
            pstr += " " + str(w)
    pstr = pstr.lstrip()
    if _verbose_:
        print(pstr)
        print()
    return pstr

def greedy_sort(permutation):
    reversal_distance = 0
    working_set = permutation.copy()
    reversal_steps = []
    if _verbose_:
        print_perumutation(working_set)
    for k in range(len(permutation)):
        if k+1 != abs(working_set[k]): # checking value (k+1)  against intended position abs(working_set[k]) also checking that k is positive w/abs()
            try:
                # if this is found, then +k is in the set, do a local reverse, and change the sign... which means you have to do a 2nd revers on k alone
                # if it is not found, then look for -k in the except block
                # print("here4")
                z = working_set.index(k+1)
                if k == 0:
                    working_set = working_set[0:k] + [-y for y in working_set[z::-1]] + working_set[z+1:]
                else:
                    working_set = working_set[0:k] + [-y for y in working_set[z:k-1:-1]] + working_set[z+1:]
                reversal_steps.append(print_perumutation(working_set))
                working_set[k] = working_set[k] * -1
                reversal_steps.append(print_perumutation(working_set))
                reversal_distance += 2
            except ValueError:
                try:
                    # if this is found then -k is in the set, do a local reverse, and change the sign ... which means that the resulting k is positive, as desired, no more change needed
                    z = working_set.index(-(k+1))
                    # print("here3")
                    if k == 0:
                        working_set = working_set[0:k] + [-y for y in working_set[z::-1]] + working_set[z+1:]
                    else:
                        working_set = working_set[0:k] + [-y for y in working_set[z:k-1:-1]] + working_set[z+1:]
                    reversal_steps.append(print_perumutation(working_set))
                    reversal_distance += 1
                except ValueError:
                    # if k is not found, then the permutation is wonky, since 1...n (+/-) has to be in the set of length n
                    raise ValueError("{0} and -{0} not found in permutation".format(str(k)))
        if k+1 != working_set[k]:
            # print("here2")
            working_set[k] = working_set[k] * -1
            reversal_steps.append(print_perumutation(working_set))
            reversal_distance += 1

    if len(working_set) != working_set[len(working_set)-1]:
        # print("here {0} {1}".format(str(len(working_set)-1), str(working_set[len(working_set)-1])))
        working_set[len(working_set)-1] *= -1
        reversal_distance += 1
        reversal_steps.append(print_perumutation(working_set))

    if _verbose_:
        print(reversal_distance)
    return reversal_steps



if __name__ == '__main__':
    start = time.process_time()
    if len(sys.argv) < 1:
        print("Expected input:\n[str: filename path]\n\nfile contents:\n[string: nucleotide permutation]\n-v = verbose, -vv = regular output -vvv = debug")

    with open(sys.argv[1]) as f:
        permstr = f.readline().rstrip()

    permutation = []
    for x in permstr.split():
        permutation.append(int(x))

    for a_idx in range(2, 3, 1):
        if len(sys.argv) > a_idx:
            if sys.argv[a_idx] == "-v":
                _verbose_ = True
            elif sys.argv[a_idx] == "-vv":
                _verbose_ = True
                _timed_output_ = True
            elif sys.argv[a_idx] == "-vvv":
                _verbose_ = True
                _timed_output_ = True
                _debug_ = True

    results = greedy_sort(permutation)
    for r in results:
        print(r)

    end = time.process_time()
    print("Time: {0}".format(end-start))

