
_verbose_ = False
_debug_ = False

class Scoring:
    def __init__ (self, match_reward, mismatch_penalty, indel_penalty, scoring_matrix_file=""):
        self.match_reward = match_reward
        self.mismatch_penalty = mismatch_penalty
        self.indel_penalty = indel_penalty
        self.using_scoring_matrix = (len(scoring_matrix_file.rstrip()) != 0)

        if self.using_scoring_matrix:
            self.scoring_matrix = self.load_scoring(scoring_matrix_file)
    
    def print_scoring(self):
        print("{0} {1} {2}".format(str(self.match_reward), str(self.mismatch_penalty), str(self.indel_penalty)))

    def get_match_val(self, protien1, protien2):
        if self.using_scoring_matrix:
            return self.scoring_matrix[protien1][protien2]
        elif protien1 == protien2:
            return self.match_reward
        else:
            return self.mismatch_penalty

    def load_scoring(self, scoring_file_loc):
        if _debug_:
            print(scoring_file_loc)

        scoring_strs = []
        with open(scoring_file_loc) as f:
            scoring_strs = f.readlines()

        scoring = {}
        scoring_fast_index = []
        horz_keys = "".join(scoring_strs[0].split())
        virt_keys = ""

        for i in range(1, len(scoring_strs)):
            row = scoring_strs[i].split()
            key = row[0]
            virt_keys = virt_keys + key
            row_dict = {}
            row_fast_index = []
            for j in range(1, len(row)):
                horz_key = horz_keys[j-1]
                row_dict[horz_key] = int(row[j])
                row_fast_index.append(int(row[j]))
            scoring[key] = row_dict
            scoring_fast_index.append(row_fast_index)

        if _verbose_:
            print_matrix(scoring_fast_index, virt_keys, horz_keys, 5)

        return scoring


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
        print("Matrix will be too large to analyze.")
        print("Truncating to two subgraphs of first " + str(truncate_height) + " and last " + str(truncate_height) + ".")
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

