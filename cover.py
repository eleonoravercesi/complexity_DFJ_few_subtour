import gurobipy as gp
from gurobipy import GRB

from itertools import chain, combinations
def powerset(iterable):
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

def partition_subsets(partition, verbose=False):
    collection = set()
    for J in powerset( range( len(partition) ) ):
        if 0 < len(J) and len(J) < len(partition):
            subset = [ i for j in J for i in partition[j] ]
            if verbose:
                print(subset)
            collection.add( frozenset(subset) )
    return collection

def partitions_subsets(partitions):
    collection = set()
    for partition in partitions:
        for subset in partition_subsets(partition):
            collection.add(subset)
    return list(collection)

def set_cover_subroutine(partitions, verbose=True):
    m = gp.Model()
    if not verbose:
        m.Params.OutputFlag = 0
    
    collection = partitions_subsets(partitions)
    s = m.addVars(len(collection), vtype=GRB.BINARY)
    
    m.setObjective( gp.quicksum(s), GRB.MINIMIZE )
    
    for partition in partitions:
        # for each partition, must pick at least one SEC to cut it off
        m.addConstr( gp.quicksum( s[i] for i in range(len(collection)) if collection[i] in partition_subsets(partition) ) >= 1 )
    
    m.optimize()
    
    assert m.solCount > 0, print("partitions =",partitions,"\n collection =",collection)
    
    num_SECs = round(m.objVal)
    if verbose:
        print(f"To cut off these partitions, require at least {num_SECs} SECs.")
    
    # re-optimize to get smallest SECs
    cardinality = round( m.objVal )
    m.addConstr( gp.quicksum(s) == cardinality )
    m.setObjective( gp.quicksum( len(collection[i]) * s[i] for i in range(len(collection)) ), GRB.MINIMIZE )
    m.optimize()
    
    if verbose:
        print(f"Among sufficient collections of {num_SECs} SECS, the smallest total size |S_1|+|S_2|+...+|S_k*| is {m.objVal}.")
    
    return [ collection[i] for i in range(len(collection)) if s[i].x > 0.5 ]