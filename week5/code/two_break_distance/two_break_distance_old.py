#! /usr/bin/python3

import sys
import time
import math

_verbose_ = False
_timed_output_ = False
_debug_ = False


def check_graph(graph):
    if len(graph) > 0:
        start = graph[0]
        prev = graph[0][1]
        for i in range(1, len(graph)):
            if prev[0] == "b":
                if "f" + prev[1:] != graph[i][0]:
                    raise ValueError("idx: {0}. {1} needs to correspond with {2}".format(str(i), prev, graph[i][0]))
            if prev[1] == "f":
                if "b" + prev[1:] != graph[i][0]:
                    raise ValueError("idx: {0}. {1} needs to correspond with {2}".format(str(i), prev, graph[i][0]))
            prev = graph[i][1]
        if prev[0] == "b":
            if "f" + prev[1:] != start[0]:
                raise ValueError("wraparound. {1} needs to correspond with {2}".format(str(i), prev, graph[i][0]))
        if prev[1] == "f":
            if "b" + prev[1:] != start[0]:
                raise ValueError("wraparound. {1} needs to correspond with {2}".format(str(i), prev, graph[i][0]))
        

def create_graph(genome):
    gen_graph = []
    for gen_cycle in genome:
        if len(gen_cycle) > 0:
            gen_cycle_graph = []
            start = gen_cycle[0]
            end = gen_cycle[len(gen_cycle)-1]

            prev = None
            for synteny in gen_cycle:
                if synteny == start:
                    if synteny > 0:
                        prev = "f" + str(synteny)
                    else:
                        prev = "b" + str(-1 * synteny)
                    continue

                if synteny > 0:
                    node = (prev, "b" + str(synteny))
                    prev = "f" + str(synteny)
                else:
                    node = (prev, "f" + str(-1 * synteny))
                    prev = "b" + str(-1 * synteny)
                gen_cycle_graph.append(node)

            if end > 0:
                if start > 0:
                    node = ("f" + str(end), "b" + str(start))
                else:
                    node = ("f" + str(end), "f" + str(-1 * start))
            else:
                if start > 0:
                    node = ("b" + str(-1 * end), "b" + str(start))
                else:
                    node = ("b" + str(-1 * end), "f" + str(-1 * start))
            gen_cycle_graph.append(node)
            check_graph(gen_cycle_graph)
            gen_graph.append(gen_cycle_graph)
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

    p_q_cycles = []
    while unvisited_nodes:
        current_node = None
        p_q_cycle = []
        for graph in genome_graphs[0]:
            for node in graph:
                if node not in genome_visited_nodes[0] and current_node is None:
                    current_node = node
            if current_node is not None:
                continue
        if current_node is None:
            raise ValueError("Should have found node. Some problem with loop constants.")
        
        start_node = current_node
        genome_graph_idx = 1
        current_node_match_idx = 1
        while current_node not in p_q_cycle:
            next_node = None
            for graph in genome_graphs[genome_graph_idx]:
                for node in graph:
                    if node[0] == current_node[current_node_match_idx] and next_node is None:
                        next_node = node
                        current_node_match_idx = 1
                    elif node[1] == current_node[current_node_match_idx] and next_node is None:
                        next_node = node
                        current_node_match_idx = 0

                if next_node is not None:
                    continue
                
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






def calc_two_break_distance(genomes, synteny_count):
    if _debug_:
        print("synteny count = {0}".format(str(synteny_count)))
    genome_graphs = create_graphs(genomes)
    p_q_cycles_count = analyze_p_q_cycles(genome_graphs, synteny_count)
    two_break_distance = synteny_count - p_q_cycles_count
    # two_break_distance = -1
    return two_break_distance

def organize_inputs(input_file):
    with open(input_file) as f:
        genome_1 = f.readline().rstrip()
        genome_2 = f.readline().rstrip()
    
    g1_cycles = genome_1.split(")(")
    g2_cycles = genome_2.split(")(")
    cycles_set= [g1_cycles, g2_cycles]

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
        print(genome_1)
        print(genome_2)
        print(synteny_count)
        print(gen_as_ints)

    return gen_as_ints, synteny_count


if __name__ == '__main__':
    start = time.process_time()
    if len(sys.argv) < 1:
        print("Expected input:\n[str: filename path]\n\nfile contents:\n[string: nucleotide permutation]\n-v = verbose, -vv = regular output -vvv = debug")


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

    genomes,synteny_count = organize_inputs(sys.argv[1])


    results = calc_two_break_distance(genomes,synteny_count)
    print(results)

    end = time.process_time()
    print("Time: {0}".format(end-start))
