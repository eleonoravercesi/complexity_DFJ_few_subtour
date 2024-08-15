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
import more_itertools as mit
import elkai

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


def tour_cost(tour, cost):
    n = len(tour)
    start = tour[0]
    end = tour[n-1]
    return sum( cost[tour[i]][tour[i+1]] for i in range(n-1) ) + cost[end][start]

# def all_partitions_simple(lst):
#     for k in range(1, len(lst)+1):
#         for partition in mit.set_partitions(lst, k):
#             yield partition

def all_partitions(lst):
    # k=1 case...
    for partition in mit.set_partitions(lst, 1):
        yield partition
    for k in range(2, len(lst)+1):
        # first try partitions like [ [0], [1], [2], [3,4,5,6,7] ]
        # that are primarily singletons
        for subset in combinations(lst, k-1):
            yield [ [lst[i]] for i in subset ] + [ [ i for i in lst if i not in subset ] ]
        
        # now, try the others
        for partition in mit.set_partitions(lst, k):
            if max( len(part) for part in partition ) + (k-1) < len(lst):
                yield partition
                
def minimalize_partition(G, partition, tsp_cost, verbose=True):
    # lkh[J] = lkh ( G[\cup_{j \in J} partition[j] ] )
    lkh = dict()
    for Js in all_partitions(range(len(partition))):
        if 1 < len(Js) and len(Js) < len(partition):
            #if verbose:
            #    print("trying Js =",Js)
            if len(Js) == 2:
                for J in Js:
                    new_part = list()
                    for j in J:
                        new_part += partition[j]
                    
                    cost = [ [ 0 if i==j else G.edges[i,j]['cost'] for j in new_part ] for i in new_part ]
                    tour = elkai.solve_int_matrix(cost)
                    lkh[frozenset(J)] = tour_cost(tour,cost)

            lkh_cost = sum( lkh[frozenset(J)] for J in Js )
            if lkh_cost < tsp_cost:
                new_partition = list()
                for J in Js:
                    new_part = list()
                    for j in J:
                        new_part += partition[j]
                    new_partition.append(new_part)
                if verbose:
                    #print("lkh =",lkh)
                    print("Found smaller partition:",len(partition),"->",len(new_partition))
                    #print("old_partition =",partition)
                    #print("new_partition =",new_partition)
                return new_partition
        
    # partition could not be made smaller
    return partition


def is_two_factor(G):
    return all( G.degree(i) == 2 for i in G.nodes )


def find_2SECs(G, tour):
    
    S_family = set()
    n = len(tour)
    assert n == G.number_of_nodes()
    root = nx.utils.arbitrary_element(G.nodes)
    
    em = dict()
    for i,j in G.edges:
        em[i,j] = (i,j)
        em[j,i] = (i,j)
        
    tour_edges = [ em[tour[i-1],tour[i]] for i in range(n) ]
    
    for p1 in range(n):
        for p2 in range(n):
            e1 = em[tour[p1],tour[(p1+1)%n]]
            e2 = em[tour[p2],tour[(p2+1)%n]]
            old_cost = G.edges[e1]['cost'] + G.edges[e2]['cost']
            
            if p2 == ( (p1+1)%n ) or p1 == ( (p2+1)%n ):
                continue
                
            e3 = em[tour[(p1+1)%n],tour[p2]]
            e4 = em[tour[(p2+1)%n],tour[p1]]
            new_cost = G.edges[e3]['cost'] + G.edges[e4]['cost']
            
            if new_cost >= old_cost:
                continue
                
            new_edges = tour_edges.copy()
            new_edges.remove(e1)
            new_edges.remove(e2)
            new_edges.append(e3)
            new_edges.append(e4)
            
            GN = G.edge_subgraph(new_edges)
            if is_two_factor(GN):
                comp = list(nx.connected_components(GN))
                assert len(comp)==2
                for S in comp:
                    case1 = (len(S) < n/2)
                    case2 = (len(S) == n/2 and root in S)
                    if case1 or case2:
                        S_family.add(frozenset(S))
    return S_family


# can we remote any three edges from tour, add three back, and get 2-factor with less cost than TSP tour?
def find_3SECs(G, tour, twoSECs):
    
    partitions = set()
    n = len(tour)
    assert n == G.number_of_nodes()
    root = nx.utils.arbitrary_element(G.nodes)
    
    em = dict()
    for i,j in G.edges:
        em[i,j] = (i,j)
        em[j,i] = (i,j)
        
    tour_edges = [ em[tour[i-1],tour[i]] for i in range(n) ]
    
    for p1 in range(n):
        for p2 in range(n):
            for p3 in range(n):
                e1 = em[tour[p1],tour[(p1+1)%n]]
                e2 = em[tour[p2],tour[(p2+1)%n]]
                e3 = em[tour[p3],tour[(p3+1)%n]]
                old_cost = G.edges[e1]['cost'] + G.edges[e2]['cost'] + G.edges[e3]['cost']

                if p1 == p2 or p2 == p3 or p3 == p1:
                    continue
                
                if p2 == ( (p1+1)%n ) or p3 == ( (p2+1)%n ) or p1 == ( (p3+1)%n ):
                    continue
                    
                if p1 == ( (p2+1)%n ) or p2 == ( (p3+1)%n ) or p3 == ( (p1+1)%n ):
                    continue

                e4 = em[tour[(p2+1)%n],tour[p1]]
                e5 = em[tour[(p3+1)%n],tour[p2]]
                e6 = em[tour[(p1+1)%n],tour[p3]]
                new_cost = G.edges[e4]['cost'] + G.edges[e5]['cost'] + G.edges[e6]['cost']

                if new_cost >= old_cost:
                    continue

                new_edges = tour_edges.copy()
                new_edges.remove(e1)
                new_edges.remove(e2)
                new_edges.remove(e3)
                new_edges.append(e4)
                new_edges.append(e5)
                new_edges.append(e6)

                GN = G.edge_subgraph(new_edges)
                if is_two_factor(GN):
                    comp = list(nx.connected_components(GN))
                    #print("comp =",comp)
                    assert len(comp)==3
                    partitions.add(frozenset(frozenset(p) for p in comp))
#     print("partitions =")
#     for p in partitions:
#         print(p)
        
    final_partitions = set()
        
    for S_family in partitions:
        found = False
        for S in S_family:
            if S in twoSECs:
                found = True
                break

        for S1 in S_family:
            for S2 in S_family:
                if S1 != S2:
                    S = frozenset(list(S1)+list(S2))
                    if S in twoSECs:
                        found = True
                        break

        if not found:
            final_partitions.add(S_family)
        
#     print("final_partitions =")
#     for p in final_partitions:
#         print(p)
    return final_partitions


def ialg(G, minimalize=True, smart_initialize=True, verbose=False):
    '''
    Implement the ialg algorithm as described in the paper "On the complexity of the Dantzig-Fulkerson-Johnson TSP formulation for few subtour constraints"
    :param G: nx.Graph, with weight on the edges labelled as "cost"
    :param minimalize: If True, attempts to reduce the partition size (# subtours) before branching
    :param verbose: If True, print the subtours found every time a solution is found and the size of the problem every 2**n iterations
    :return: List of list containing the subtours found and the number of subtours found
    '''
    if verbose:
        print("*****************************")
        print("Computing TSP cost")
        print("*****************************")
    (tsp_cost, tour_edges) = mip(G, subtour_callbacks=True, return_edges=True, verbose=verbose) 
    tsp_cost = round(tsp_cost)
    
    m = gp.Model()
    if not verbose:
        m.Params.OutputFlag = 0
    x = m.addVars(G.edges, vtype=GRB.BINARY)

    m.setObjective(gp.quicksum(G.edges[e]['cost'] * x[e] for e in G.edges), GRB.MINIMIZE)

    # Each city should touch the tour twice (enter-and-leave)
    m.addConstrs(gp.quicksum(x[e] for e in G.edges if e in G.edges(i)) == 2 for i in G.nodes)

    m.update()

    # store branch-and-bound nodes in a (min) heap
    B = list()
    partitions = dict()

    # (priority=n*len(S_family)+num_comp, size=len(S_family), S_family)
    if smart_initialize:
        (start, end) = tour_edges[0]
        tour = list( nx.dfs_preorder_nodes(G.edge_subgraph(tour_edges[1:]), source=start) )
        SECs_set = find_2SECs(G, tour)
        SECs = [ list(S) for S in SECs_set ]
        partitions[2] = set()
        for S in SECs:
            VS = [ i for i in G.nodes if i not in S ]
            partition = { frozenset(S), frozenset(VS) }
            partitions[2].add( frozenset(partition) )
        if verbose:
            print("With smart initialization, we begin with #SECs =",len(SECs))
            print("They are:")
            for S in SECs:
                print(S)
        threeps = find_3SECs(G, tour, SECs_set)
        partitions[3] = set()
        for pn in threeps:
            partition = { frozenset(p) for p in pn }
            partitions[3].add(frozenset(partition))
        print("partitions =",partitions)
        root = (0, len(SECs), SECs)
    else:
        root = (0, 0, list())
    heapq.heappush(B, root)

    # run branch-and-bound
    num_nodes = 0

    while len(B) > 0:
        (priority, size, S_family) = heapq.heappop(B)
        num_nodes += 1

        # add subtour constraints
        c = m.addConstrs(
            gp.quicksum(x[e] for e in E(G, S_family[p])) <= len(S_family[p]) - 1 for p in range(len(S_family)))

        # optimize
        m.optimize()
        tour_edges = [e for e in G.edges if x[e].x > 0.5]
        cost = round(m.objVal)

        # remove subtour constraints (so m can be re-used for later solves)
        m.remove(c)
        m.update()

        # cost equals TSP cost?
        if cost == tsp_cost:
            assert size == len(S_family)
            if verbose:
                print("Found a solution with #SECs:", size)
                print("Specifically, they are:")
                for S in S_family:
                    print(S)
                print("S_family = ", S_family)
                print("And finally, we have partitions =",partitions)
            return S_family, len(S_family)

        # check if we can use partition from partition pool
        components = False
        if nx.number_connected_components(G.edge_subgraph(tour_edges)) > 2:
            for partition in partitions[2]:
                feasible = True

                for S in S_family:
                    found_crossing = False
                    for part in partition:
                        ctS = sum( 1 for i in part if i in S )
                        ctVS = sum( 1 for i in part if i not in S )
                        if ctS>0 and ctVS>0:
                            found_crossing = True
                            break

                    if not found_crossing:
                        feasible = False
                        break

                if feasible:
                    print("Found partition from pool:",partition)
                    components = [ list(part) for part in partition ]
                    print("Its components are:",components)
                    break
                    
            if not components:
                for partition in partitions[3]:
                    feasible = True

                    for S in S_family:
                        found_crossing = False
                        for part in partition:
                            ctS = sum( 1 for i in part if i in S )
                            ctVS = sum( 1 for i in part if i not in S )
                            if ctS>0 and ctVS>0:
                                found_crossing = True
                                break

                        if not found_crossing:
                            feasible = False
                            break

                    if feasible:
                        print("Found partition from pool:",partition)
                        components = [ list(part) for part in partition ]
                        print("Its components are:",components)
                        break
        
        # cannot find suitable partition in pool. Use 2-factor instead
        if not components:
            components = list(nx.connected_components(G.edge_subgraph(tour_edges)))
            if minimalize and len(components) > 2:
                components = minimalize_partition(G, components, tsp_cost, verbose=verbose)
        
        num_comp = len(components)
        indices = list(range(num_comp))
        
        # add partition to pool
        try:
            partitions[num_comp]
        except:
            partitions[num_comp] = set()
            
        partition = { frozenset(comp) for comp in components }
        partitions[num_comp].add( frozenset(partition) )

        # print status update after 1, 2, 4, 8, 16, 32, ... BB nodes
        #if is_power_of_two(num_nodes) and verbose:
        if verbose:
            print("num_bb_nodes, num_subtour_constrs, num_conn_comp =", num_nodes, size, num_comp)
            
        # for each subset of subtours (S_1, S_2, ..., S_t ), create a subproblem
        #       that has a new constraint for S = S_1 \cup S_2 \cup ... \cup S_t
        for s in range(1, num_comp):
            comb = combinations(indices, s)
            for subset in comb:

                # create vertex subset S
                S = list()
                for i in subset:
                    S += list(components[i])

                # It suffices to impose the constraint for either S or V\S.
                # We choose to pick S with:
                #   1. |S| < |V|/2, or
                #   2. |S| = |V|/2 and 0 \in S
                case1 = (2 * len(S) < G.number_of_nodes())
                case2 = (2 * len(S) == G.number_of_nodes() and 0 in S)

                if case1 or case2:
                    S_family.append(S)
                    #priority = G.number_of_nodes() * size + num_comp
                    priority = tsp_cost * (size+1) - cost
                    new_node = (priority, size + 1, S_family.copy())
                    heapq.heappush(B, new_node)
                    S_family.pop()


# Generic MIP function that can solve:
# 1. min-weight 2-factor 
#      mip(G)
# 2. min-weight 2-factor subject to select subtour constraints
#      mip(G, initial_subtours=initial_subtours)
# 3. DFJ model with subtour callbacks
#      mip(G, subtour_callbacks=True)
#      option: initial_subtours
#      option: one_cut=True adds one cut per callback
# ...
#
def mip(G, initial_subtours=list(), subtour_callbacks=False, one_cut=False, return_components=False, return_edges=False, verbose=True):
    
    assert not return_components or not return_edges, "Cannot return both. Pick one."
    
    # start with the 2-matching relaxation
    # Create model object
    m = gp.Model()
    if not verbose:
        m.Params.OutputFlag = 0

    # Create variable for each edge
    x = m.addVars( G.edges, vtype=GRB.BINARY )
    
    # Objective function
    m.setObjective( gp.quicksum( G.edges[e]['cost'] * x[e] for e in G.edges), GRB.MINIMIZE )

    # Add degree-two constraint for each vertex
    m.addConstrs( gp.quicksum( x[e] for e in G.edges if e in G.edges(i)  ) == 2 for i in G.nodes )
    
    # Add initial subtour constraints (if any)
    for S in initial_subtours:
        m.addConstr( gp.quicksum( x[e] for e in E(G, S) ) <= len(S) - 1 )

    # add subtour elimination constraints in callback
    m._x = x
    if subtour_callbacks:
        m.Params.LazyConstraints = 1
        m._G = G
        m._one_cut = one_cut
        m._subtours = list()
        m._callback = subtour_elimination 
        m.optimize(m._callback)
        if verbose:
            print("number of subtours:",len(m._subtours))
            print("subtours =",m._subtours)
    else:
        m.optimize()
        
    if return_edges:
        chosen_edges = [ e for e in G.edges if m._x[e].x > 0.5 ]
        return ( m.objVal, chosen_edges )
    elif return_components:
        chosen_edges = [ e for e in G.edges if m._x[e].x > 0.5 ]
        return ( m.objVal, list(nx.connected_components(G.edge_subgraph(chosen_edges))) )
    else:
        return m.objVal

    
# a function to separate subtour elimination constraints
def subtour_elimination(m, where):
    
    # check if LP relaxation at this branch-and-bound node has an integer solution
    if where != GRB.Callback.MIPSOL: 
        return
        
    # retrieve the LP solution
    xval = m.cbGetSolution(m._x)
    G = m._G

    # which edges are selected?
    tour_edges = [ e for e in G.edges if xval[e] > 0.5 ]

    # edges already form a tour; no cut needed
    if nx.is_connected( G.edge_subgraph( tour_edges )):
        return

    # for each subtour, add a cut
    root = nx.utils.arbitrary_element(G.nodes)
    n = G.number_of_nodes()
    for S in sorted(nx.connected_components( G.edge_subgraph( tour_edges ) ), key=len):
        case1 = ( len(S) < n/2 )
        case2 = ( len(S) == n/2 and root in S )
        if case1 or case2:
            m.cbLazy( gp.quicksum( m._x[e] for e in E(G,S) ) <= len(S) - 1 )
            m._subtours.append(S)
            if m._one_cut:
                return

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