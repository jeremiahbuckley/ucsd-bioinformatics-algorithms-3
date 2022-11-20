#! /usr/bin/python3

import sys
import time

if __name__ == '__main__':
    start = time.process_time()

    if len(sys.argv) < 3:
        print("Expected input:\n[str: test files path]\n[str: filename]\n\nfile contents:\n[string: nucleotide]\n[string: nucleotide]\n[string: nucleotide]\n-v = verbose, -vv = debug")

    filename = sys.argv[2]
    path = sys.argv[1] + "/inputs/" + sys.argv[2]
    path = path.replace("//", "/")
    with open(path) as f:
        input_n_1 = f.readline().rstrip()
        input_n_2 = f.readline().rstrip()
        input_n_3 = f.readline().rstrip()

    
    path = path.replace("input", "output")
    with open(path) as f:
        output_n_0 = f.readline().rstrip() #score,for this purpose, discard
        output_n_1 = f.readline().rstrip()
        output_n_2 = f.readline().rstrip()
        output_n_3 = f.readline().rstrip()

    out_n_1 = out_n_2 = out_n_3 = ""
    for i in range(len(output_n_1)):
        if output_n_1[i] != "-":
            out_n_1 += output_n_1[i]

    for i in range(len(output_n_2)):
        if output_n_2[i] != "-":
            out_n_2 += output_n_2[i]

    for i in range(len(output_n_3)):
        if output_n_3[i] != "-":
            out_n_3 += output_n_3[i]

    if input_n_1 != out_n_1:
        print("error in {0}. no match for string {1}\n{2}\n{3}".format(sys.argv[1], "1", input_n_1, out_n_1))
    if input_n_2 != out_n_2:
        print("error in {0}. no match for string {1}\n{2}\n{3}".format(sys.argv[1], "2", input_n_2, out_n_2))
    if input_n_3 != out_n_3:
        print("error in {0}. no match for string {1}\n{2}\n{3}".format(sys.argv[1], "3", input_n_3, out_n_3))
  
    end = time.process_time()
    print("Time: {0}".format(end-start))
