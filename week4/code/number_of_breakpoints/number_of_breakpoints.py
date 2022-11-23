#! /usr/bin/python3

import sys
import time
import math

_verbose_ = False
_timed_output_ = False
_debug_ = False


def num_of_breakpoints(permutation):
    prev_val = permutation[0]
    bps = 0
    for next_val in permutation[1:]:
        if next_val - prev_val != 1:
            print("{0} {1}".format(prev_val, next_val))
            bps += 1
        prev_val = next_val
    
    return bps


if __name__ == '__main__':
    start = time.process_time()
    if len(sys.argv) < 1:
        print("Expected input:\n[str: filename path]\n\nfile contents:\n[string: nucleotide permutation]\n-v = verbose, -vv = regular output -vvv = debug")

    with open(sys.argv[1]) as f:
        permstr = f.readline().rstrip()

    permutation = [0]
    for x in permstr.split():
        permutation.append(int(x))
    permutation.append(len(permutation))

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

    results = num_of_breakpoints(permutation)
    print(results)

    end = time.process_time()
    print("Time: {0}".format(end-start))
