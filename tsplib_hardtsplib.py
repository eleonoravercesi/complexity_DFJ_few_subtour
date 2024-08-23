#!/usr/bin/env python
# coding: utf-8
import networkx as nx

# # Running the `ialg` algorithm on the TSPLIB and the HardTSPLIB instances
# 
# Here, we want to compute the minimum number of SECs needed to prove optimality for some small famous instances in the TSPLIB [1] and some hard-to-solve instances in the HardTSPLIB [2]. We will use the `ialg` algorithm to compute the minimum number of SECs needed to prove optimality for these instances.
# 
# ### References
# [1] Reinelt, Gerhard. "TSPLIB—A traveling salesman problem library." ORSA journal on computing 3.4 (1991): 376-384.
# 
# [2] Vercesi, Eleonora, et al. "On the generation of metric TSP instances with a large integrality gap by branch-and-cut." Mathematical Programming Computation 15.2 (2023): 389-416.
# 

# In[40]:


from ialg import ialg, mip
from cover import set_cover_subroutine
from utils import from_tsplib_file_to_graph
import pandas as pd
import pickle


# tsplib_instances = [("burma14", 14), ("ulysses16", 16), ("gr17", 17), ("gr21", 21), ("ulysses22", 22),  ("gr24", 24),
#                     ("fri26", 29), ("bayg29", 29), ("bays29", 29), ("dantzig42", 42), ("swiss42", 42), ("att48", 48),
#                      ("gr48", 48), ("hk48", 48), ("eil51", 51), ("berlin52", 52), ("brazil58", 58), ("st70", 70),
#                      ("eil76", 76), ("pr76", 76), ("gr96", 96), ("rat99", 99), ("kroA100", 100), ("kroB100", 100),
#                     ("kroC100", 100), ("kroD100", 100), ("kroE100", 100), ("rd100", 100)]



hard_tsplib_instances_random = [("10001_hard", 10), ("10007_hard", 10), ("10008_hard", 10), ("10010_hard", 10),
                         ("11675_hard", 11), ("12290_hard", 12), ("14850_hard", 14), ("15002_hard", 15),
                         ("15005_hard", 15), ("15007_hard", 15), ("16038_hard", 16), ("20007_hard", 20),
                            ("20009_hard", 20), ("20181_hard", 20), ("25001_hard", 25), ("25004_hard", 25),
                            ("25006_hard", 25), ("30001_hard", 30), ("30003_hard", 30), ("30005_hard", 30),
                            ("33001_hard", 33), ("35002_hard", 35), ("35003_hard", 35), ("35009_hard", 35),
                            ("40003_hard", 40), ("40004_hard", 40), ("40008_hard", 40)]

hard_tsplib_instances_from_tsplib = [("gr24_hard", 24), ("bayg29_hard", 29), ("bays29_hard", 29),
                                     ("dantzig42_hard", 42), ("swiss42_hard", 42), ("gr48_hard", 48), ("hk48_hard", 48),
                                     ("att48_hard", 48), ("eil51_hard", 51), ("berlin52_hard", 52), ("brazil58_hard", 58),
                                     ("st70_hard", 70), ("pr76_hard", 76), ("eil76_hard", 76)]



# Store the values in a dictionary
out = {}


'''
Random instances
'''
for minimalize in [True, False]:
    for instance_name, n in hard_tsplib_instances_random:
        # Parse the instance
        G = from_tsplib_file_to_graph("./data/" + instance_name)
        print("******* Instance:", instance_name, "*******")
        (S_family, size_S_family, partitions, c, runtime, num_nodes) = ialg(G, verbose=False, minimalize=minimalize)
        out[instance_name] = [S_family, size_S_family, partitions, c, runtime, num_nodes]
        # At each iteration save the out dictionary in pickle
        with open("OUT_HardTSPLIB_random_minimalize_{}.pickle".format(minimalize), "wb") as f:
            pickle.dump(out, f)

        if size_S_family == -1:
            print("Ran into time limit.")
            continue

        # check that k* >= size_S_family
        partitions_list = [ [ list(part) for part in partition ] for npts in partitions for partition in partitions[npts] ]
        smallest_S_family = set_cover_subroutine(partitions_list, verbose=True)
        print("k* =",size_S_family,"for S_family =", smallest_S_family)
        assert len(smallest_S_family) == size_S_family

        # Compute the TSP
        tsp_cost = mip(G, subtour_callbacks=True, verbose=False)

        # check that k* <= size_S_family
        two_factor_cost_with_S_family, _ = mip(G, initial_subtours=S_family, verbose=False, return_edges=True)
        assert round(two_factor_cost_with_S_family) == round(tsp_cost)
        two_factor_cost, x = mip(G, initial_subtours=smallest_S_family, verbose=False, return_edges=True)
        assert round(two_factor_cost) == round(tsp_cost)


for minimalize in [True, False]:
    for instance_name, n in hard_tsplib_instances_from_tsplib:
        # Parse the instance
        G = from_tsplib_file_to_graph("./data/" + instance_name)
        print("******* Instance:", instance_name, "*******")
        (S_family, size_S_family, partitions, c, runtime, num_nodes) = ialg(G, verbose=False, minimalize=minimalize)
        out[instance_name] = [S_family, size_S_family, partitions, c, runtime, num_nodes]
        # At each iteration save the out dictionary in pickle
        with open("OUT_HardTSPLIB_TSPLIB_minimalize_{}.pickle".format(minimalize), "wb") as f:
            pickle.dump(out, f)

        if size_S_family == -1:
            print("Ran into time limit.")
            continue

        # check that k* >= size_S_family
        partitions_list = [ [ list(part) for part in partition ] for npts in partitions for partition in partitions[npts] ]
        smallest_S_family = set_cover_subroutine(partitions_list, verbose=True)
        print("k* =",size_S_family,"for S_family =", smallest_S_family)
        assert len(smallest_S_family) == size_S_family

        # Compute the TSP
        tsp_cost = mip(G, subtour_callbacks=True, verbose=False)

        # check that k* <= size_S_family
        two_factor_cost_with_S_family, _ = mip(G, initial_subtours=S_family, verbose=False, return_edges=True)
        assert round(two_factor_cost_with_S_family) == round(tsp_cost)
        two_factor_cost, x = mip(G, initial_subtours=smallest_S_family, verbose=False, return_edges=True)
        assert round(two_factor_cost) == round(tsp_cost)

        print(" ") # Leave some space







# # ## On the impact of the minimalization procedure
#
# # ## HardTSPLIB instances
# #
# # HardTSPLIB is made of instances generated both at random and starting from instances of the TSPLIB.  We divide such instances into to, to make direct comparison with TSPLIB. `hardtsplib_instances_random` and `hardtsplib_instances_tsplib`. We will only run the algorithm on the instances that are feasible to run on our machine.
#
# # In[46]:
#
#
# hardtsplib_instances_random = [("10001_hard", 10), ("10007_hard", 10), ("10008_hard", 10), ("10010_hard", 10), ("11675_hard", 11), ("12290_hard", 12), ("14850_hard", 14), ("15002_hard", 15), ("15005_hard", 15), ("15007_hard", 15), ("16038_hard", 16), ("20004_hard", 20), ("20007_hard", 20), ("20009_hard", 20), ("20181_hard", 20), ("25001_hard", 25), ("25004_hard", 25), ("25006_hard", 25), ("30001_hard", 30), ("30003_hard", 30), ("30005_hard", 30), ("33001_hard", 33), ("35002_hard", 35), ("35003_hard", 35), ("35009_hard", 35), ("40003_hard", 40), ("40004_hard", 40), ("40008_hard", 40)]
#
# hardtsplib_instances_tsplib = [("gr24_hard", 24), ("bayg29_hard", 29), ("bays29_hard", 29), ("dantzig42_hard", 42),  ("gr48_hard", 48), ("hk48_hard", 48), ("att48_hard", 48), ("eil51_hard", 51), ("brazil58_hard", 58), ("st70_hard", 70), ("pr76_hard", 76)]
#
#
# # In[47]:
#
#
# max_instance_random = 10
# max_instance_from_tsplib = 1
#
#
# #
# # ## Hard instances derived from TSPLIB instances
#
# # In[48]:
#
#
# # Store the values in a dictionary
# out_tsplib = {}
#
#
# # In[49]:
#
#
# for instance_name, n in hardtsplib_instances_tsplib[:max_instance_from_tsplib]:
#     # Parse the instance
#     G = from_tsplib_file_to_graph("./data/" + instance_name)
#     print("******* Instance:", instance_name, "*******")
#     (S_family, size_S_family, partitions, c, runtime) = ialg(G, verbose=True)
#     out_tsplib[instance_name] = (S_family, size_S_family, partitions, c, runtime)
#
#     if size_S_family == -1:
#         print("Ran into time limit.")
#         continue
#
#     # check that k* >= size_S_family
#     partitions_list = [ [ list(part) for part in partition ] for npts in partitions for partition in partitions[npts] ]
#     smallest_S_family = set_cover_subroutine(partitions_list, verbose=False)
#     print("k* =",size_S_family,"for S_family =",smallest_S_family)
#     assert len(smallest_S_family) == size_S_family
#
#     # check that k* <= size_S_family
#     two_factor_cost = mip(G, initial_subtours=smallest_S_family, verbose=False)
#     tsp_cost = mip(G, subtour_callbacks=True, verbose=False)
#     assert round(two_factor_cost) == round(tsp_cost)
#
#     print(" ") # Leave some space
#
#
# # In this case, we can see a line-by-line comparison between the values of TSPLIB and HardTSPLIB
# #
#
# # In[61]:
#
#
# # Create a dataframe out of the dictionary out and out_tsplib
# df_list = []
# for i in range(len(out_tsplib)):
#     tsplib = list(out_tsplib.items())[i]
#     hardtsplib = list(out_tsplib.items())[i]
# df_list.append((tsplib[0], tsplib[1][1], hardtsplib[1][1], tsplib[1][3], hardtsplib[1][3]))
# df = pd.DataFrame(df_list, columns=["instance", "S_min_TSPLIB", "S_min_HardTSPLIB", "b_TSPLIB", "b_HardTSPLIB"])
# # Print the dataframe
# df
#
#
# #
# # ## Hard instances derived from random instances
#
# # In[62]:
#
#
# # Store the values in a dictionary
# out_random = {}
#
#
# # In[63]:
#
#
# for instance_name, n in hardtsplib_instances_random[:max_instance_random]:
#     # Parse the instance
#     G = from_tsplib_file_to_graph("./data/" + instance_name)
#     print("******* Instance:", instance_name, "*******")
#     (S_family, size_S_family, partitions, c, runtime) = ialg(G, verbose=True)
#     out_random[instance_name] = (S_family, size_S_family, partitions, c, runtime)
#
#     if size_S_family == -1:
#         print("Ran into time limit.")
#         continue
#
#     # check that k* >= size_S_family
#     partitions_list = [ [ list(part) for part in partition ] for npts in partitions for partition in partitions[npts] ]
#     smallest_S_family = set_cover_subroutine(partitions_list, verbose=False)
#     print("k* =",size_S_family,"for S_family =",smallest_S_family)
#     assert len(smallest_S_family) == size_S_family
#
#     # check that k* <= size_S_family
#     two_factor_cost = mip(G, initial_subtours=smallest_S_family, verbose=False)
#     tsp_cost = mip(G, subtour_callbacks=True, verbose=False)
#     assert round(two_factor_cost) == round(tsp_cost)
#
#     print(" ") # Leave some space
#
#
# # In[65]:
#
#
# # Create a dataframe out of the dictionary out
# df = pd.DataFrame([(x[0], x[1][1], x[1][3], x[1][4]) for x in out_random.items()],
#                   columns=["instance", "S_min", "b_max", "runtime"])
# # Print the dataframe
# df


# In[ ]:




