'''
Function used to reproduce the results in "On the complexity of the Dantzig-Fulkerson-Johnson TSP formulation for few subtour constraints"
Authors: E. Vercesi, A. Buchanan
'''

import heapq
import gurobipy as gp
from gurobipy import GRB
import networkx as nx
from itertools import combinations
import more_itertools as mit
import elkai
import time


def E(G,S):
    '''
    Return the edges of G that are in S
    :param G: nx.Graph
    :param S: Could be a set or a list
    :return: A list containing the edges of G that are in S
    '''
    return [(i,j) for (i,j) in G.edges if i in S and j in S]


def tour_cost(tour, cost):
    '''
    Computes length/cost of tour
    :param tour: list
    :param cost: cost matrix stored as list of lists
    :return: Cost of tour (int)
    '''
    n = len(tour)
    start = tour[0]
    end = tour[n-1]
    return sum( cost[tour[i]][tour[i+1]] for i in range(n-1) ) + cost[end][start]


def all_partitions(lst):
    '''
    Returns all partitions of lst
    :param lst: list
    :return: python generator, each time yielding a list of lists
    '''
    # k=1 case...
    for partition in mit.set_partitions(lst, 1):
        yield partition
        
    # k>=2 case...
    for k in range(2, len(lst)+1):
        
        # first try partitions like [ [0], [1], [2], [3,4,5,6,7] ] that are primarily singletons
        for subset in combinations(lst, k-1):
            yield [ [lst[i]] for i in subset ] + [ [ i for i in lst if i not in subset ] ]
        
        # now, try the others
        for partition in mit.set_partitions(lst, k):
            if max( len(part) for part in partition ) + (k-1) < len(lst):
                yield partition
                
def minimalize_partition(G, partition, tsp_cost, verbose=True):
    '''
    Merges parts of the partition such that we still have \sum_j tsp(G[V_j]) < tsp(G)
    :param G: nx.Graph
    :param partition: list of lists
    :param tsp_cost: int
    :return: list of lists
    '''
    # caches lkh values for subsets, like lkh[J] = lkh ( G[\cup_{j \in J} partition[j] ] )
    lkh = dict()
    for Js in all_partitions(range(len(partition))):
        if 1 < len(Js) and len(Js) < len(partition):
            
            # compute lkh(G[V_1]) + lkh(G[V_2]) for 2-partitions (V_1, V_2)
            if len(Js) == 2:
                for J in Js:
                    new_part = list()
                    for j in J:
                        new_part += partition[j]
                    
                    cost = [ [ 0 if i==j else G.edges[i,j]['cost'] for j in new_part ] for i in new_part ]
                    tour = elkai.solve_int_matrix(cost)
                    lkh[frozenset(J)] = tour_cost(tour,cost)

            # use cached values, rather than re-computing from scratch
            lkh_cost = sum( lkh[frozenset(J)] for J in Js )
            if lkh_cost < tsp_cost:
                new_partition = list()
                for J in Js:
                    new_part = list()
                    for j in J:
                        new_part += partition[j]
                    new_partition.append(new_part)
                if verbose:
                    print("Found smaller partition:",len(partition),"->",len(new_partition))
                return new_partition
        
    # partition could not be made smaller; return the original one
    return partition


def quick_bipartitions(G, tour):
    '''
    Drops two edges from tour, splice endpoints together, and see if resulting cycles have cost(C1) + cost(C2) < tour_cost
    :param G: nx.Graph
    :param tour: list 
    :return: set of 2-partitions (V(C1), V(C2))
    '''
    partitions = set()
    n = len(tour)
    assert n == G.number_of_nodes()
    root = nx.utils.arbitrary_element(G.nodes)
    
    # edge map stores orientation of each edge. Needed because networkx stores undirected edges as if there are directed...
    em = dict()
    for i,j in G.edges:
        em[i,j] = (i,j)
        em[j,i] = (i,j)
        
    tour_edges = [ em[tour[i-1],tour[i]] for i in range(n) ]
    
    for p1 in range(n):
        for p2 in range(p1+3,n):

            # too close?
            if abs(p2-p1)<3 or abs(p1+n-p2)<3:
                continue
            
            # pick two edges to drop
            e1 = em[tour[p1],tour[(p1+1)%n]]
            e2 = em[tour[p2],tour[(p2+1)%n]]
            old_cost = G.edges[e1]['cost'] + G.edges[e2]['cost']
            
            if p2 == ( (p1+1)%n ) or p1 == ( (p2+1)%n ):
                continue
                
            # which edges splice the path endpoints
            e3 = em[tour[(p1+1)%n],tour[p2]]
            e4 = em[tour[(p2+1)%n],tour[p1]]
            new_cost = G.edges[e3]['cost'] + G.edges[e4]['cost']
            
            if new_cost >= old_cost:
                continue
                
            path_edges = [ e for e in tour_edges if e not in {e1,e2} ]
            partition = { frozenset(comp) for comp in nx.connected_components(G.edge_subgraph(path_edges)) }
            partitions.add( frozenset(partition) )
    return partitions

def get_partition_from_pool(S_family, partitions, max_parts=4, verbose=False):
    '''
    Seeks a partition that satisfies x(E(S))<=|S|-1 for all S in S_family
    :param S_family: list of lists
    :param partitions: set of partitions, with each partition stored as a set of sets
    :param max_parts: when selecting a partition, how many parts do we allow it to have?
    :return: (if any exist) a feasible partition from partitions, else False
    '''
    for pts in range(2,max_parts+1):
        try:
            partitions[pts]
        except:
            partitions[pts] = set()
        for partition in partitions[pts]:
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
                components = [ list(part) for part in partition ]
                if verbose:
                    print("Found partition from pool:",partition)
                    print("Its components are:",components)
                return components
    return False

def ialg(G, minimalize=True, smart_initialize=True, partition_pool=True, verbose=False, time_limit=3600):
    '''
    Implement the ialg algorithm as described in the paper "The Dantzig-Fulkerson-Johnson TSP formulation is easy to solve for few subtour constraints"
    :param G: nx.Graph, with weight on the edges labeled as "cost"
    :param minimalize: If True, attempts to reduce the partition size (# subtours) before branching
    :param smart_initialize: If True, initializes the algorithm with 2-partitions obtained from the TSP tour
    :return: S_family, list of vertex subsets S for the SECs
    :return: len(S_family) if terminates; -1 if time limit exceeded
    :return: partitions, collection of *all* vertex partitions that the algorithm branched on
    :return: max_comp, the max number of components (or parts) across the partitions
    '''
    start_tsp_time = time.time()
    # Here, we set the verbosity of the mip = False to keep the output of a reasonable size
    (tsp_cost, tour_edges) = mip(G, subtour_callbacks=True, return_edges=True, verbose=False)
    tsp_cost = round(tsp_cost)
    end_tsp_time = time.time() - start_tsp_time
    if verbose:
        print("TSP compute in {} seconds. TSP cost = {}".format(end_tsp_time, tsp_cost))
    
    start_time = time.perf_counter()
    m = gp.Model()
    # Try to fix bug
    m.Params.MIPGap = 1e-16
    m.Params.FeasibilityTol = 1e-9
    m.Params.IntFeasTol = 1e-9
    # We refer to the verbosity of the ialg algorithm. All the output of Gurobi have been turned off
    m.Params.OutputFlag = 0
    x = m.addVars(G.edges, vtype=GRB.BINARY)

    # minimize the sum of edge weights
    m.setObjective(gp.quicksum(G.edges[e]['cost'] * x[e] for e in G.edges), GRB.MINIMIZE)

    # Each city should touch the tour twice (enter-and-leave)
    m.addConstrs(gp.quicksum(x[e] for e in G.edges if e in G.edges(i)) == 2 for i in G.nodes)

    m.update()

    # store search tree nodes in a (min) heap
    B = list()
    partitions = dict()

    if smart_initialize:
        (start, end) = tour_edges[0]
        tour = list( nx.dfs_preorder_nodes(G.edge_subgraph(tour_edges[1:]), source=start) )
        partitions[2] = quick_bipartitions(G, tour)
        twoSECs = [ list(partition)[0] for partition in partitions[2] ]
        if verbose:
            print("With smart initialization, we begin with #SECs =",len(twoSECs))
            print("They are:")
            for S in twoSECs:
                print(S)
        max_comp = 2 if len(twoSECs)>0 else 0
        root = (0, len(twoSECs), twoSECs)
    else:
        max_comp = 0
        root = (0, 0, list())
    
    # run the main algorithm
    heapq.heappush(B, root)
    num_nodes = 0
    arb = nx.utils.arbitrary_element(G.nodes) # arbitrary vertex, used to break ties

    while len(B) > 0:
        
        # check time limit
        elapsed = time.perf_counter() - start_time
        #if verbose:
        #    print("elapsed time =", elapsed)
        if elapsed > time_limit:
            if verbose:
                print("Exceeded time limit. Exiting")
            return (list(), -1, partitions, max_comp, 3600)
        
        (priority, size, S_family) = heapq.heappop(B)
        num_nodes += 1

        # add subtour constraints
        c = m.addConstrs(
            gp.quicksum(x[e] for e in E(G, S_family[p])) <= len(S_family[p]) - 1 for p in range(len(S_family)))

        # optimize
        m.optimize()
        tour_edges = [ e for e in G.edges if x[e].x > 0.5 ]
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
            return (S_family, len(S_family), partitions, max_comp, time.perf_counter() - start_time)

        # first, try partition from partition pool
        components = False
        if partition_pool and nx.number_connected_components(G.edge_subgraph(tour_edges)) > 2:
            components = get_partition_from_pool(S_family, partitions, verbose=verbose)
                    
        # cannot find suitable partition in pool. Use 2-factor instead
        if not components:
            components = list(nx.connected_components(G.edge_subgraph(tour_edges)))
            if minimalize and len(components) > 2:
                components = minimalize_partition(G, components, tsp_cost, verbose=verbose)
        
        num_comp = len(components)
        max_comp = max( max_comp, num_comp )
        indices = list(range(num_comp))
        
        # add partition to pool
        try:
            partitions[num_comp]
        except:
            partitions[num_comp] = set()
            
        partition = { frozenset(comp) for comp in components }
        partitions[num_comp].add( frozenset(partition) )

        if verbose:
            print("num_bb_nodes, num_subtour_constrs, num_conn_comp =", num_nodes, size, num_comp)
            
        # for each subset of subtours (S_1, S_2, ..., S_t ), create a subproblem
        #       that has a new constraint for S = S_1 \cup S_2 \cup ... \cup S_t
        for s in range(1, num_comp):
            comb = combinations(indices, s)
            for subset in comb:

                # create vertex subset S
                S = [ j for i in subset for j in components[i] ]

                # It suffices to impose the constraint for either S or V\S.
                # We choose to pick S with:
                #   1. |S| < |V|/2, or
                #   2. |S| = |V|/2 and 0 \in S
                case1 = (2 * len(S) < G.number_of_nodes())
                case2 = (2 * len(S) == G.number_of_nodes() and arb in S)

                if case1 or case2:
                    S_family.append(S)
                    #priority = G.number_of_nodes() * size + num_comp
                    big_M = sum( G.edges[e]['cost'] for e in G.edges )
                    priority = big_M * (size+1) - cost
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