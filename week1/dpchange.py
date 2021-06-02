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


if __name__ == '__main__':
    start = time.process_time()

    if len(sys.argv) < 2:
        print("Expected input: 1 param, format '<int_money>\n<list_int_coins>'")

    with open(sys.argv[1]) as f:
        money = int(f.readline())
        coins_str = f.readline()
        coins = [int(i) for i in coins_str.split(',')]

    min_coins = dpchange(money, coins)
    print(min_coins)

    end = time.process_time()
    print("Time: {0}".format(end-start))
