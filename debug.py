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
    print(x)

print("------")

for x in smallest_S_family:
    print(x)