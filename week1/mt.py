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
        print("  ".join([str(n) for n in list_of_lists[i]]))

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

if __name__ == '__main__':
    start = time.process_time()

    if len(sys.argv) < 2:
        print("Expected input: 1 param, format '<int_width> <int_height>\n<int n x (m+1) matrix>\n'-'\n<int (n+1) x m matrix>'")

    height = 0
    width = 0
    down_weights = []
    right_weights = []
    with open(sys.argv[1]) as f:
        hw_str = f.readline()
        height = [int(i) for i in hw_str.split(' ')][0]
        width = [int(i) for i in hw_str.split(' ')][1]
        down_weights = []
        for i in range(height):
            w_str = f.readline()
            down_weights_i = [int(i) for i in w_str.split(' ')]
            down_weights.append(down_weights_i)

        f.readline()

        right_weights = []
        for i in range(height+1):
            w_str = f.readline()
            right_weights_i = [int(i) for i in w_str.split(' ')]
            right_weights.append(right_weights_i)


    print("{0} {1}".format(str(height), str(width)))
    for r in down_weights:
        print(" ".join([str(i) for i in r]))
    print("---")
    for r in right_weights:
        print(" ".join([str(i) for i in r]))

    longest_path = manhattan_tourist(height, width, down_weights, right_weights)
    print(longest_path)

    end = time.process_time()
    print("Time: {0}".format(end-start))
