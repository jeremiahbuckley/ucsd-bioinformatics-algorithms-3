def print_matrices(max_vals_matrices, backtrack_matrices, nucleotide_h, nucleotide_w):
    print("values matrices")
    for m in max_vals_matrices:
        print()
        print_matrix(m, nucleotide_h, nucleotide_w)
    print("backtrack matrices")
    for m in backtrack_matrices:
        print()
        print_matrix(m, nucleotide_h, nucleotide_w)


def print_matrix(list_of_lists, str_v, str_h, truncate_height = 20):
    height = len(list_of_lists)
    width = len(list_of_lists[0])
    column_gap = 4

    truncate_width = truncate_height

    if height > truncate_height or width > truncate_width:
        print("Matrix too large for useful print. Truncating to two subgraphs of first " + str(truncate_height) + " and last " + str(truncate_height) + ".")
        print(str_h)
        print(str_v)

        print("first " + str(truncate_height))
        out_str = " {0}".format("".join([ch.rjust(column_gap) for ch in str_h[0:truncate_width]]))
        print(out_str[:79]) # this is to keep formating at < 80 chars
        for i in range(truncate_height):
            #print(list_of_lists[i])
            #print(str_v[i])
            #print()
            out_str = "".join([str(n).rjust(column_gap) for n in list_of_lists[i][:truncate_width]])
            out_str = out_str[:79] # this is to keep formating at < 80 chars
            out_str = out_str.replace("↖︎", " ↖︎") #on the terminal this character is incorrectly right-justified by one-too-few spaces, this is not true in other formats
            print("{0}{1}".format(str_v[i], out_str))
            # if i == 0:
            #     print("    {0}".format(out_str))
            # else:
            #     print("{0}   {1}".format(str_v[i-1], out_str))

        print("last " + str(truncate_height))
        out_str = " {0}".format("".join([ch.rjust(column_gap) for ch in str_h[-1 * truncate_width:]]))
        print(out_str[:79]) # this is to keep formating at < 80 chars
        for i in range(height-truncate_height, height):
            #print(list_of_lists[i])
            #print(str_v[i])
            #print()
            out_str = "".join([str(n).rjust(column_gap) for n in list_of_lists[i][-1 * truncate_width:]])
            out_str = out_str[:79] # this is to keep formating at < 80 chars
            out_str = out_str.replace("↖︎", " ↖︎") #on the terminal this character is incorrectly right-justified by one-too-few spaces, this is not true in other formats
            print("{0}{1}".format(str_v[i], out_str))
            # if i == 0:
            #     print("    {0}".format(out_str))
            # else:
            #     print("{0}   {1}".format(str_v[i-1], out_str))

    else:
        print(" {0}".format("".join([ch.rjust(column_gap) for ch in str_h])))
        for i in range(height):
            # print(list_of_lists[i])
            # print(str_v[i])
            # print()
            out_str = "".join([str(n).rjust(column_gap) for n in list_of_lists[i]])
            out_str = out_str.replace("↖︎", " ↖︎") #on the terminal this character is incorrectly right-justified by one-too-few spaces, this is not true in other formats
            print("{0}{1}".format(str_v[i], out_str))
            # if i == 0:
            #     print("    {0}".format(out_str))
            # else:
            #     print("{0}   {1}".format(str_v[i-1], out_str))

