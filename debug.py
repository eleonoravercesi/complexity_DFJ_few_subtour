import networkx as nx

from ialg import ialg, mip
from utils import from_tsplib_file_to_graph
from cover import set_cover_subroutine


# Load the instance
instance = "./data/rat99.tsp"
G = from_tsplib_file_to_graph(instance)

(S_family, size_S_family, partitions_keep, c, runtime) = ialg(G, verbose=True)

partitions_list = [[list(part) for part in partition ] for npts in partitions_keep for partition in partitions_keep[npts]]
smallest_S_family = set_cover_subroutine(partitions_list, verbose=True)

two_factor_cost_smallest, x = mip(G, initial_subtours=smallest_S_family, verbose=False, return_edges=True)
H = G.edge_subgraph(x)
print(list(nx.connected_components(H)))
two_factor_cost_S_family, _ = mip(G, initial_subtours=S_family, verbose=False, return_edges=True)

# And the tsp
tsp_cost, x_tour = mip(G, subtour_callbacks=True, verbose=False, return_edges=True)
print("Two-factor cost of smallest S-family: ", two_factor_cost_smallest)
print("Two-factor cost of S-family: ", two_factor_cost_S_family)
print("TSP cost: ", tsp_cost)


for x in S_family:
    print(sorted(x))

print("------")

for x in smallest_S_family:
    print(sorted(x))

#smallest_S_family.append([72, 73, 74, 75, 76, 77, 78, 81, 82, 83, 84, 85, 86, 87, 90, 91, 92, 93, 94])
smallest_S_family.append( [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 40, 41, 42, 43])
two_factor_cost_S_family_pushed, _ = mip(G, initial_subtours=smallest_S_family, verbose=False, return_edges=True)
print("Two-factor cost of smallest S-family: ", two_factor_cost_smallest)
print("Two-factor cost of S-family: ", two_factor_cost_S_family)
print("TSP cost: ", tsp_cost)
print("Two-factor cost of S-family pushed: ", two_factor_cost_S_family_pushed)
