'''
Function used to reproduce the results in "On the complexity of the Dantzig-Fulkerson-Johnson TSP formulation for few subtour constraints"
Authors: A. Buchanan, E. Vercesi
'''

import heapq
import gurobipy as gp
from gurobipy import GRB
import math
import networkx as nx
from itertools import combinations
import numpy as np
from scipy.spatial.distance import cdist
import tsplib95
import os

def sample(n, distribution):
    '''
    Sample n points from a given distribution
    :param n: Number of points to sample
    :param distribution: Could be "uniform_points_in_unit_square", "uniform", "exponential"
    :return: nx.Graph with the cost of the edges labelled as "cost"
    '''
    if distribution not in ["uniform_points_in_unit_square", "uniform", "exponential"]:
        raise ValueError("Invalid distribution. Choose between 'uniform_points_in_unit_square', 'uniform', 'exponential'")
    if distribution == "uniform_points_in_unit_square":
        X = np.random.uniform(0, 1, (n, 2))
        # Compute the euclidean distance matrix
        D = cdist(X, X, 'euclidean')
        # Create a graph
        G = nx.Graph()
        for i in range(n):
            for j in range(i + 1, n):
                G.add_edge(i, j, cost=round(1000*D[i, j])) # Multiply every value by 1000 and round
        return G
    if distribution == "uniform":
        D = np.random.uniform(0, 1000, (n, n)) # You sample the whole matrix at once, but you only use the upper triangular part
        G = nx.Graph()
        for i in range(n):
            for j in range(i + 1, n):
                G.add_edge(i, j, cost=round(D[i, j]))  # Round to the nearest integer
        return G
    if distribution == "exponential":
        D = np.random.exponential(1/999, (n, n)) # You sample the whole matrix at once, but you only use the upper triangular part
        G = nx.Graph()
        for i in range(n):
            for j in range(i + 1, n):
                G.add_edge(i, j, cost=round(1000*D[i, j]))  # Multiply every value by 1000 and round just to avoid having 0 costs
        return G


def is_power_of_two(n):
    '''
    Check if n is a power of 2
    :param n: an integer
    :return: True, if n is a power of 2, False otherwise
    '''
    return (n != 0) and (n & (n - 1) == 0)

def E(G,S):
    '''
    Return the edges of G that are in S
    :param G: nx.Graph
    :param S: Could be a set or a list
    :return: A list containing the edges of G that are in S
    '''
    return [(i,j) for (i,j) in G.edges if i in S and j in S]

def ialg(G, verbose=False):
    '''
    Implement the ialg algorithm as described in the paper "On the complexity of the Dantzig-Fulkerson-Johnson TSP formulation for few subtour constraints"
    :param G: nx.Graph, with weight on the edges labelled as "cost"
    :param verbose: If True, print the subtours found every time a solution is found and the size of the problem every 2**n iterations
    :return: List of list containing the subtours found and the number of subtours found
    '''
    m = gp.Model()
    m.Params.OutputFlag = 0
    x = m.addVars(G.edges, vtype=GRB.BINARY)

    m.setObjective(gp.quicksum(G.edges[e]['cost'] * x[e] for e in G.edges), GRB.MINIMIZE)

    # Each city should touch the tour twice (enter-and-leave)
    m.addConstrs(gp.quicksum(x[e] for e in G.edges if e in G.edges(i)) == 2 for i in G.nodes)

    m.update()

    # initialize upper bound to infinity
    incumbent = None
    UB = math.inf

    # store branch-and-bound nodes in a (min) heap
    B = list()

    # (priority=n*len(S_family)+num_comp, size=len(S_family), S_family)
    root = (0, 0, list())
    heapq.heappush(B, root)


    # run branch-and-bound
    num_nodes = 0

    while len(B) > 0:
        (priority, size, S_family) = heapq.heappop(B)
        num_nodes += 1

        if size >= UB:
            continue  # prune by bound

        # add subtour constraints
        c = m.addConstrs(
            gp.quicksum(x[e] for e in E(G, S_family[p])) <= len(S_family[p]) - 1 for p in range(len(S_family)))

        # optimize
        m.optimize()
        tour_edges = [e for e in G.edges if x[e].x > 0.5]

        # remove subtour constraints (so m can be re-used for later solves)
        m.remove(c)
        m.update()

        # produces a tour?
        if nx.is_connected(G.edge_subgraph(tour_edges)):
            assert size == len(S_family)
            if verbose:
                print("Found a solution with this many subtours:", size)
                print("Specifically, they are:")
                for S in S_family:
                    print(S)
                print("S_family = ", S_family)
            UB = size
            incumbent = S_family.copy()
            continue  # prune by feasibility

        # prune *children* by bound
        if size + 1 >= UB:
            continue

        # add subproblems to heap
        sorted_components = list(sorted(nx.connected_components(G.edge_subgraph(tour_edges)), key=len))
        num_comp = len(sorted_components)
        indices = list(range(num_comp))

        # print status update after 1, 2, 4, 8, 16, 32, ... BB nodes
        if is_power_of_two(num_nodes) and verbose:
            print("num_bb_nodes, num_subtour_constrs, num_conn_comp =", num_nodes, size, num_comp)

        # for each subset of subtours (S_1, S_2, ..., S_t ), create a subproblem
        #       that has a new constraint for S = S_1 \cup S_2 \cup ... \cup S_t
        for s in range(1, num_comp):
            comb = combinations(indices, s)
            for subset in comb:

                # create vertex subset S
                S = list()
                for i in subset:
                    S += list(sorted_components[i])

                    # It suffices to impose the constraint for either S or V\S.
                # We choose to pick S with:
                #   1. |S| < |V|/2, or
                #   2. |S| = |V|/2 and 0 \in S
                case1 = (2 * len(S) < G.number_of_nodes())
                case2 = (2 * len(S) == G.number_of_nodes() and 0 in S)

                if case1 or case2:
                    S_family.append(S)
                    priority = G.number_of_nodes() * size + num_comp
                    new_node = (priority, size + 1, S_family.copy())
                    heapq.heappush(B, new_node)
                    S_family.pop()

    S_family = incumbent
    return S_family, len(S_family)

#########################################
# Custom function to parse TSPLIB files #
#########################################
# This functions use the standard library tsplib95 to parse TSPLIB files
def debug_tsplib95(path):
    '''
    Remove the last line of a TSPLIB file if it contains EOF and if rise an error when trying to parse it
    :param path: path to file
    :return: a problem object from the tsplib library
    '''
    try:
        problem = tsplib95.load(path)
    except:
        # Remove EOF
        clean_problem = open(path, "r").readlines()[:-1]
        TEMP = open("tmp.tsp", "w+")
        for line in clean_problem:
            TEMP.write(line)
        TEMP.close()
        problem = tsplib95.load("tmp.tsp")
        os.remove("tmp.tsp")
    return problem


def parse_TSPLIB_file(problem):
    '''
    Parse a TSPLIB file and return the (linearized) matrix of costs and the number of nodes. It's prepared to deal with the most common the edge_weight_types and edge_weight_formats
    :param problem: an instance of problem for the tsplib95 library
    :return: C (list), n (int)
    '''
    ewt = problem.as_name_dict()['edge_weight_type']
    n = problem.dimension
    Cin = problem.edge_weights
    C = []
    if ewt == "EXPLICIT":
        ewf = problem.as_name_dict()['edge_weight_format']
        if ewf == "FULL_MATRIX":
            for i in range(n):
                for j in range(i + 1, n):
                    C.append(Cin[i][j])
        elif ewf == 'UPPER_ROW':
            for x in Cin:
                for y in x:
                    C.append(y)
        elif ewf in ['LOWER_DIAG_ROW', 'UPPER_DIAG_ROW']:
            for x in Cin:
                for y in x:
                    if y != 0:
                        C.append(y)
            if ewf == 'LOWER_DIAG_ROW':
                C = list(reversed(C))
    else:
        X = np.asarray(list(problem.node_coords.values()))
        if ewt in ["EUC_2D", "EUC_3D"]:
            for i in range(n):
                for j in range(i + 1, n):
                    C.append(tsplib95.distances.euclidean(X[i], X[j]))
        if ewt in ["MAN_2D", "MAN_3D"]:
            for i in range(n):
                for j in range(i + 1, n):
                    C.append(tsplib95.distances.manhattan(X[i], X[j]))
        if ewt in ["MAX_2D", "MAX_3D"]:
            for i in range(n):
                for j in range(i + 1, n):
                    C.append(tsplib95.distances.maximum(X[i], X[j]))
        if ewt == 'ATT':
            for i in range(n):
                for j in range(i + 1, n):
                    C.append(tsplib95.distances.pseudo_euclidean(X[i], X[j]))
        if ewt == 'GEO':
            for i in range(n):
                for j in range(i + 1, n):
                    C.append(tsplib95.distances.geographical(X[i], X[j]))
        if ewt == 'CEIL_2D':
            for i in range(n):
                for j in range(i + 1, n):
                    C.append(tsplib95.distances.euclidean(X[i], X[j], round=math.ceil))
    assert len(C) == n * (n - 1) / 2
    return C, n

def make_matrix(C, n):
    '''
    Create a matrix from a list of costs
    :param C: A list of costs of dimension n*(n-1)/2
    :param n: The dimension n
    :return: a np.array of dimension n x n
    '''
    # Create a matrix
    M = np.zeros((n, n))
    # Set a counter
    cont = 0
    for i in range(n):
        for j in range(i + 1, n):
            M[i, j] = C[cont]
            M[j, i] = C[cont]
            cont += 1
    return M

def from_tsplib_file_to_graph(filename):
    '''
    Create a graph from a TSPLIB file
    :param filename: path of the tsplib file
    :return: nx.Graph() : a complete graph; each edge has a cost attribute
    '''
    if "tsp" not in filename:
        filename = filename + ".tsp"
    C, n = parse_TSPLIB_file(debug_tsplib95(filename))
    C = make_matrix(C, n)
    G = nx.Graph()
    for i in range(n):
        for j in range(i + 1, n):
            G.add_edge(i, j, cost=C[i, j])
    return G
