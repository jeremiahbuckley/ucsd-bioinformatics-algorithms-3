#! /usr/bin/python3

import sys
import time
import math

_verbose_ = False
_timed_output_ = False
_debug_ = False


def check_graph(graph):
    if len(graph) < 2:
        return

    start_cycle = graph[0][0]
    prev = graph[0][1]
    for i in range(1, len(graph)):

        next_node_offset = 1 if prev % 2 == 1 else -1

        if prev + next_node_offset != graph[i][0]:
            if prev + next_node_offset != start_cycle:
                print(i)
                print(prev)
                print(graph)
                raise ValueError("idx: {0}. {1} needs to correspond with {2}".format(str(i), prev, graph[i][0]))
            elif i < len(graph)-1:
                start_cycle = graph[i][0]
                prev = graph[i][1]
                if _debug_:
                    print("cycle closed at {0}. beginning new cycle.".format(str(i)))

        prev = graph[i][1]
    
    next_node_offset = 1 if prev % 2 == 1 else -1
    if prev + next_node_offset != start_cycle:
        raise ValueError("idx: {0}. {1} needs to correspond with {2}".format("end-start", prev, start_cycle))

def create_graph(genome):
    gen_graph = []
    for gen_cycle in genome:
        if len(gen_cycle) > 0:
            start = gen_cycle[0]
            end = gen_cycle[len(gen_cycle)-1]

            prev = None
            for synteny in gen_cycle:
                if synteny == start:
                    start_second_endpoint_offset = 0 if start > 0 else -1
                    prev = abs(start) * 2 + start_second_endpoint_offset
                    continue

                synteny_first_endpoint_offset = -1 if synteny > 0 else 0
                synteny_second_endpoint_offset = 0 if synteny > 0 else -1
                node = (prev, abs(synteny) * 2 + synteny_first_endpoint_offset)
                prev = abs(synteny) * 2 + synteny_second_endpoint_offset

                gen_graph.append(node)

            # link end node second-endpoint to start node first-endpoint
            end_second_endpoint_offset = 0 if end > 0 else -1
            start_first_endpoint_offset = -1 if start > 0 else 0
            node = (abs(end) * 2 + end_second_endpoint_offset, abs(start) * 2 + start_first_endpoint_offset)

            gen_graph.append(node)
    if _debug_:
        print("created graph")
        print(genome)
        print(gen_graph)
    check_graph(gen_graph)
    return gen_graph

def create_graphs(genomes):
    graphs = []
    for gen in genomes:
        graphs.append(create_graph(gen))
    if _debug_:
        print(graphs)
    return graphs

def analyze_p_q_cycles(genome_graphs, synteny_count):
    genome_visited_nodes = [{}, {}]

    unvisited_nodes = len(genome_graphs[0]) > 0

    if _debug_:
        print("analyze_p_q_cycles")
        print(genome_graphs[0])
        print(genome_graphs[1])
        print(synteny_count)
        print()
    p_q_cycles = []
    while unvisited_nodes:
        current_node = None
        p_q_cycle = []
        for node in genome_graphs[0]:
            if node not in genome_visited_nodes[0]:
                current_node = node
                break

        if current_node is None:
            raise ValueError("Should have found node. Some problem with loop constants.")
        
        start_node = current_node
        genome_graph_idx = 1
        current_node_other_endpoint_idx = 1
        while current_node not in p_q_cycle:
            next_node = None
            for node in genome_graphs[genome_graph_idx]:
                for i in range(2):
                    if node[i] == current_node[current_node_other_endpoint_idx]:
                        next_node = node
                        current_node_other_endpoint_idx = (i + 1) % 2
                        break
                    
                if next_node is not None:
                    break
                
            # if current_node == start_node:
                # p_q_cycle.append(current_node)
                # if current_node in genome_visited_nodes[0]:
                #     raise ValueError("Visited node ({0}, {1}), part of graph {2} twice.".format(str(current_node[0]), str(current_node[1]), str(0)))
                # genome_visited_nodes[0][current_node] = 0
            p_q_cycle.append(current_node)
            if current_node in genome_visited_nodes[(genome_graph_idx +1) % 2]:
                raise ValueError("Visited node ({0}, {1}), part of graph {2} twice.".format(str(current_node[0]), str(current_node[1]), str(0)))
            genome_visited_nodes[(genome_graph_idx +1) % 2][current_node] = 0

            if next_node is None:
                raise ValueError("Expected next_node to be populated.")
            # p_q_cycle.append(next_node)
            # if next_node in genome_visited_nodes[genome_graph_idx]:
            #     raise ValueError("Visited node ({0}, {1}), part of graph {2} twice.".format(str(current_node[0]), str(current_node[1]), str(genome_graph_idx)))
            # genome_visited_nodes[genome_graph_idx][next_node] = 0

            if _debug_:
                print(next_node)
            current_node = next_node
            genome_graph_idx = (genome_graph_idx + 1) % 2
        
        if _debug_:
            print("p_q_cycle:")
            print(p_q_cycle)
            print(genome_visited_nodes[0])
            print(genome_visited_nodes[1])

        p_q_cycles.append(p_q_cycle)

        unvisited_nodes = len(genome_visited_nodes[0]) < synteny_count

    if _debug_:
        print(len(p_q_cycles))
    return len(p_q_cycles)






def calc_two_break_distance_from_synteny_list(genomes, synteny_count):
    if _debug_:
        print("synteny count (synteny list)= {0}".format(str(synteny_count)))
    genome_graphs = create_graphs(genomes)
    return calc_two_break_distance_from_node_graph(genome_graphs, synteny_count)

def calc_two_break_distance_from_node_graph(genome_graphs, synteny_count):
    if _debug_:
        print("synteny count (genome graphs) = {0}".format(str(synteny_count)))
    p_q_cycles_count = analyze_p_q_cycles(genome_graphs, synteny_count)
    two_break_distance = synteny_count - p_q_cycles_count
    # two_break_distance = -1
    return two_break_distance

def perform_two_break_transform(genomes, process, synteny_count):
    if len(genomes) != 2 or len(genomes[0]) == 0 or len(genomes[1]) == 0:
        raise ValueError("genomes are set up incorrectly.")

    found_next_target = False
    first_node_to_break = None
    if _debug_:
        print("before find next_target_node")
        print(genomes[0])
        print(genomes[1])
    # for i in range(len(genomes[1])-1,-1,-1):
    #     node = genomes[1][i]
    #     if not found_next_target:
    #         found_match_in_gen_1 = False
    #         for gen_1_node in genomes[0]:
    #             if ((gen_1_node[0] == node[0] and gen_1_node[1] == node[1]) or (gen_1_node[1] == node[0] and gen_1_node[0] == node[1])):
    #                 found_match_in_gen_1 = True
    #                 break

    #         first_node_to_break = node
    #         found_next_target = not found_match_in_gen_1
    for i in range(len(genomes[0])-1,-1, -1):
        end_state_node = genomes[0][i]
        if _debug_:
            print(end_state_node)
        if not found_next_target:
            partial_match_in_gen_2 = False
            partial_match_node = None
            for gen_2_node in genomes[1]:
                if _debug_:
                    print(gen_2_node)
                # there are 4 partial matches, but to game the system I'm specifically only going to look for a second-half partial match
                # if ((gen_2_node[0] == end_state_node[0] and gen_2_node[1] != end_state_node[1]) or (gen_2_node[0] == end_state_node[1] and gen_2_node[1] != end_state_node[0]
                #  or (gen_2_node[0] != end_state_node[0] and gen_2_node[1] == end_state_node[1]) or (gen_2_node[0] != end_state_node[1] and gen_2_node[1] == end_state_node[0]))):
                if ((end_state_node[0] != gen_2_node[0] and end_state_node[1] == gen_2_node[1]) or (end_state_node[0] != gen_2_node[1] and end_state_node[1] == gen_2_node[0])):
                    if _debug_:
                        print("found partial")
                    partial_match_in_gen_2 = True
                    partial_match_node = gen_2_node
                    break
            
            if partial_match_in_gen_2:
                if _debug_:
                    print("found partial")

                first_node_to_break = partial_match_node
                found_next_target = True

    if _debug_:
        print("onward with node")
        print(first_node_to_break)

    second_node_to_break = None
    if  first_node_to_break is not None:
        second_break_first_endpoint = -1

        for gen_1_node in genomes[0]:
            if gen_1_node[0] == first_node_to_break[1]:
                second_break_first_endpoint = gen_1_node[1]
            elif gen_1_node[1] == first_node_to_break[1]:
                second_break_first_endpoint = gen_1_node[0]

        for gen_2_node in genomes[1]:
            if gen_2_node[0] == second_break_first_endpoint:
                second_break_second_endpoint = gen_2_node[1]
                second_node_to_break = gen_2_node
            elif gen_2_node[1] == second_break_first_endpoint:
                second_break_second_endpoint = gen_2_node[0]
                second_node_to_break = gen_2_node
    
    second_node_to_insert = (first_node_to_break[1], second_break_first_endpoint)
    first_node_to_insert = (second_break_second_endpoint, first_node_to_break[0])

    # second_node_to_break = None
    # if  first_node_to_break is not None:
    #     second_break_first_endpoint = -1

    #     for gen_1_node in genomes[0]:
    #         if gen_1_node[0] == first_node_to_break[0]:
    #             second_break_first_endpoint = gen_1_node[1]
    #         elif gen_1_node[1] == first_node_to_break[0]:
    #             second_break_first_endpoint = gen_1_node[0]

    #     for gen_2_node in genomes[1]:
    #         if gen_2_node[0] == second_break_first_endpoint:
    #             second_break_second_endpoint = gen_2_node[1]
    #             second_node_to_break = gen_2_node
    #         elif gen_2_node[1] == second_break_first_endpoint:
    #             second_break_second_endpoint = gen_2_node[0]
    #             second_node_to_break = gen_2_node
    
    # first_node_to_insert = (first_node_to_break[0], second_break_first_endpoint)
    # second_node_to_insert = (second_break_second_endpoint, first_node_to_break[1])


    if _debug_:
        print("rejigger graph to accomidate for breaks")
        print("break targets:")
        print(first_node_to_break)
        print(second_node_to_break)
        print("new nodes:")
        print(first_node_to_insert)
        print(second_node_to_insert)
        print(genomes[1])
    new_genome_2 = []
    prev_node = None
    genome_2_nodes_attached = []
    insert_nodes_attached = []

    prev_node, node_idx, insert_nodes_idx = find_next_unattached_node(genomes[1], genome_2_nodes_attached, insert_nodes_attached, first_node_to_break, second_node_to_break, first_node_to_insert, second_node_to_insert)
    if node_idx >= 0:
        genome_2_nodes_attached.append(node_idx)
    if insert_nodes_idx >= 0:
        insert_nodes_attached.append(insert_nodes_idx)
    
    new_genome_2.append(prev_node)
    gen_2_len = len(genomes[1])
    next_node = None
    # start_cycle = prev_node[0] - 1 if prev_node[0] % 2 == 0 else prev_node[0] + 1
    start_cycle = prev_node[0]
    prev_endpoint = prev_node[1]
    if _debug_:
        print("starting new_genome_2")
        print(new_genome_2)
        print(prev_node)
        print(start_cycle)
        print()

    while len(new_genome_2) < gen_2_len:
        next_node = None
        new_cycle = False
        nn_offset =  1 if prev_endpoint % 2 == 1 else -1

        next_node, node_idx, insert_nodes_idx = find_next_node_target_in_collection(prev_endpoint + nn_offset, genomes[1], genome_2_nodes_attached, insert_nodes_attached, first_node_to_break, second_node_to_break, first_node_to_insert, second_node_to_insert)
        if node_idx >= 0:
            genome_2_nodes_attached.append(node_idx)
        if insert_nodes_idx >= 0:
            insert_nodes_attached.append(insert_nodes_idx)
        # else:
        #     if next_node[1] + 1 if next_node[1] % 2 == 1 else -1 == start_cycle:
        #         new_genome_2.append(next_node)

        #         start_cycle = prev_node[0]
        #         prev_endpoint = prev_node[1]
        #         new_genome_2.append(prev_node)
        #         new_cycle = True
        #         if _debug_:
        #             print("found cycle at {0} to {1}".format(str(prev_endpoint), str(start_cycle)))
        
        # if we Still couldn't find it, maybe we're at the end of a cycle
        if next_node is None:
            if prev_endpoint + nn_offset == start_cycle:
                if _debug_:
                    print("found cycle at {0} to {1}".format(str(prev_endpoint), str(start_cycle)))
                prev_node, node_idx, insert_nodes_idx = find_next_unattached_node(genomes[1], genome_2_nodes_attached, insert_nodes_attached, first_node_to_break, second_node_to_break, first_node_to_insert, second_node_to_insert)
                if node_idx >= 0:
                    genome_2_nodes_attached.append(node_idx)
                if insert_nodes_idx >= 0:
                    insert_nodes_attached.append(insert_nodes_idx)

                start_cycle = prev_node[0]
                prev_endpoint = prev_node[1]
                new_genome_2.append(prev_node)
                new_cycle = True
            else:
                raise ValueError("next_node can't be None, loop constants incorrect")


        if not new_cycle:
            new_genome_2.append(next_node)
            prev_node = next_node
            prev_endpoint = next_node[1]
        
        if _debug_:
            print("new genome 2 building...")
            print(new_genome_2)
            print(next_node)
            print(prev_endpoint)

    # if next_node is None:
    #     raise ValueError("exited loop with next_node = None")

    # new_genome_2.append(next_node)

    genomes.remove(genomes[1])
    genomes.append(new_genome_2)

    process.append(new_genome_2)

    if _debug_:
        print("exiting one transform")
        print(new_genome_2)

    tb_distance = calc_two_break_distance_from_node_graph(genomes, synteny_count)

    return genomes, process, tb_distance


def find_next_node_target_in_collection(target, gen_graph, already_used_graph_nodes, already_used_insert_nodes, first_node_to_break, second_node_to_break, first_node_to_insert, second_node_to_insert):
    gen_graph_node_idx = -1
    insert_node_idx = -1
    if 0 not in already_used_insert_nodes and first_node_to_insert[0] == target:
        return first_node_to_insert, gen_graph_node_idx, 0
    if 0 not in already_used_insert_nodes and first_node_to_insert[1] == target:
        # reverse node
        return (first_node_to_insert[1], first_node_to_insert[0]), gen_graph_node_idx, 0
    if 1 not in already_used_insert_nodes and second_node_to_insert[0] == target:
        return second_node_to_insert, gen_graph_node_idx, 1
    if 1 not in already_used_insert_nodes and second_node_to_insert[1] == target:
        # reverse node
        return (second_node_to_insert[1], second_node_to_insert[0]), gen_graph_node_idx, 1

    for i in range(len(gen_graph)-1,-1,-1):
        if i not in already_used_graph_nodes:
            node = gen_graph[i]
            if node == first_node_to_break or node == second_node_to_break or i in already_used_graph_nodes:
                continue
            if node[0] == target:
                return node, i, insert_node_idx
            if node[1] == target:
                # reverse node
                return (node[1], node[0]), i, insert_node_idx
    
    return None, -1, -1
            
def find_next_unattached_node(gen_graph, already_used_graph_nodes, already_used_insert_nodes, first_node_to_break, second_node_to_break, first_node_to_insert, second_node_to_insert):
    gen_graph_node_idx = -1
    insert_node_idx = -1
    for i in range(len(gen_graph)-1,-1,-1):
        if i not in already_used_graph_nodes:
            node = gen_graph[i]
            if node == first_node_to_break or node == second_node_to_break or i in already_used_graph_nodes:
                continue
            return node, i, insert_node_idx
    
    if 0 not in already_used_insert_nodes:
        return first_node_to_insert, gen_graph_node_idx, 0
    else:
        return second_node_to_insert, gen_graph_node_idx, 1


def transform_second_genome_to_first(genomes, synteny_count):
    gen_graphs = create_graphs(genomes)
    tb_distance = calc_two_break_distance_from_node_graph(gen_graphs, synteny_count)
    tb_distances = [tb_distance]
    process = [gen_graphs[1]]
    if _debug_:
        print("process, tb_distance (start) = {0}".format(str(tb_distance)))
        for g in gen_graphs:
            print(g)
        print()
        for p in process:
            print(p)
        print()
    while tb_distance > 0:
        gen_graphs, process, new_distance = perform_two_break_transform(gen_graphs, process, synteny_count)
        tb_distance = new_distance
        tb_distances.append(tb_distance)
        if _debug_:
            print("process, tb_distance = {0}".format(str(tb_distance)))
            for g in gen_graphs:
                print(g)
            print()
            for p in process:
                print(p)
            print()

    return process, tb_distances

def transform_graph_list_to_synteny_string(graph_list):
    # end_cycle = graph_list[len(graph_list) - 1][1]
    synteny_string = ""
    start_cycle = next = graph_list[0][0]
    prev = graph_list[0][0] - 1 if graph_list[0][0] % 2 == 0 else graph_list[0][0] + 1
    for i in range(0,len(graph_list)):
        next = graph_list[i][0]

        synteny = ""
        if (prev % 2 == 0 and (prev-1 != next)) or (next % 2 == 0 and (prev+1 != next)) or (abs(prev-next) != 1):
        # if abs(prev-next) != 1:
            if (prev % 2 == 0 and (prev-1 != start_cycle)) or (start_cycle % 2 == 0 and (prev+1 != start_cycle)):
            # if abs(prev - start_cycle) != 1:
                raise ValueError("can't correctly transform graph into synteny because prev / next distance != 1 and prev / start != 1. prev: {0} next: {1} start_cycle: {2}".format(str(prev), str(next), str(start_cycle)))
            start_cycle = next
            prev = next -1 if next % 2 == 0 else next + 1
            synteny += ")("
        if prev > next:
            if prev % 2 == 1 and next != start_cycle:
                raise ValueError("value needed to be multiple of 2 (prev) {0}".format(str(prev)))
            synteny += " -" + str(int(prev/2))
        else:
            if next % 2 == 1 and next != start_cycle:
                raise ValueError("value needed to be multiple of 2 (next) {0}".format(str(next)))
            synteny += " +" + str(int(next/2))
        synteny_string += synteny

        prev = graph_list[i][1]

    synteny_string = "(" + synteny_string + ")"
    synteny_string = synteny_string.replace("( ", "(").replace(" )","")

    return synteny_string


def transform_graph_lists_to_synteny_strings(graph_lists):
    results = []
    for g in graph_lists:
        results.append(transform_graph_list_to_synteny_string(g))
    return results

def organize_inputs(input_file):
    with open(input_file) as f:
        genome_2 = f.readline().rstrip()
        genome_1 = f.readline().rstrip()
    
    g1_cycles = genome_1.split(")(")
    g2_cycles = genome_2.split(")(")
    cycles_set= [g1_cycles, g2_cycles]
    for cycles in cycles_set:
        first = 0
        last = len(cycles) - 1
        for i in range(len(cycles)):
            if i != first:
                cycles[i] = "(" + cycles[i]
            if i != last:
                cycles[i] = cycles[i] + ")"
    if _debug_:
        print(cycles_set)
    return cycles_set


def convert_cycles_strings_to_int_lists(cycles_set):
    gen_as_ints = []
    for cycles in cycles_set:
        synteny_count = 0 # assumption, the sum of all synteny counts in the cycles in a cycles_sets is the same for all cycles_sets (aka genomes we're looking at)
        token_group = []
        for cycle in cycles:
            # if _debug_:
            #     print(cycle)
            tokens = cycle.split()
            # if _debug_:
            #     print(tokens)
            if tokens[0][0] == "(":
                tokens[0] = tokens[0][1:]
            if tokens[len(tokens)-1][len(tokens[len(tokens)-1])-1] == ")":
                tokens[len(tokens)-1] = tokens[len(tokens)-1][0:len(tokens[len(tokens)-1])-1]
            ints_set = [int(t) for t in tokens]
            # if _debug_:
            # print(ints_set)
            token_group.append(ints_set)
            synteny_count += len(ints_set)

        gen_as_ints.append(token_group)

    if _debug_:
        print(cycles_set[0])
        print(cycles_set[1])
        print(synteny_count)
        print(gen_as_ints)

    return gen_as_ints, synteny_count


if __name__ == '__main__':
    start = time.process_time()
    if len(sys.argv) < 1:
        print("Expected input:\n[str: filename path]\n\nfile contents:\n[string: nucleotide permutation]\n[string: nucleotide permutation]\n-v = verbose, -vv = regular output -vvv = debug")


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

    cycles_set = organize_inputs(sys.argv[1])
    genomes,synteny_count = convert_cycles_strings_to_int_lists(cycles_set)

    process, tb_distances = transform_second_genome_to_first(genomes,synteny_count)
    results = transform_graph_lists_to_synteny_strings(process)
    for r in results:
        print(r)

    end = time.process_time()
    print("Time: {0}".format(end-start))
