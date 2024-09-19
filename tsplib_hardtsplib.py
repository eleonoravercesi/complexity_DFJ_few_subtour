'''
For large instances, it is better to have a command line interface to run the experiments.
'''
# Suppress future warning of pandas
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
# Things that you actually need
from ialg import ialg, mip
from cover import set_cover_subroutine
from utils import from_tsplib_file_to_graph
import pandas as pd
from tqdm import tqdm

# Redirect the standard output to a file
import sys
sys.stdout = open("./output/tsplib_hardtsplib.log", "w")

'''
TSPLIB
'''
#tsplib_instances = [("burma14", 14), ("ulysses16", 16), ("gr17", 17), ("gr21", 21), ("ulysses22", 22),  ("gr24", 24),
#                    ("fri26", 29), ("bayg29", 29), ("bays29", 29), ("dantzig42", 42), ("swiss42", 42)]
# In case you want to test also these instances
large_tsplib_instance = [("att48", 48), ("gr48", 48), ("hk48", 48), ("eil51", 51), ("berlin52", 52), ("brazil58", 58), ("st70", 70),
                     ("eil76", 76), ("pr76", 76), ("gr96", 96), ("rat99", 99), ("kroA100", 100), ("kroB100", 100),
                    ("kroC100", 100), ("kroD100", 100), ("kroE100", 100), ("rd100", 100)]


# ### With minimalization
minim = True

print(" *** Running the ialg algorithm with minimalization on the TSPLIB instances ***")
# Create a dataframe with keys tsplib instances
df = pd.DataFrame(index = [x[0] for x in large_tsplib_instance], columns=["S_min", "b_prime", "runtime", "BeB_n_nodes"])
for instance_name, n in tqdm(large_tsplib_instance):
    # Parse the instance
    G = from_tsplib_file_to_graph("./data/" + instance_name)
    print("Instance:", instance_name)
    (S_family, size_S_family, partitions, c, runtime, num_nodes) = ialg(G, verbose=False, minimalize=minim)
    df.loc[instance_name] = [size_S_family, c, runtime, num_nodes]

    if size_S_family == -1:
        print("\tRan into time limit.")
        continue
    else:
        print(
            "\t After {} seconds and {} BeB nodes, we found {} SECs. Specifically, they are:".format(runtime, num_nodes,
                                                                                                 size_S_family))
        for S in S_family:
            print("\t\t", S)

    # check that k* >= size_S_family
    partitions_list = [ [ list(part) for part in partition ] for npts in partitions for partition in partitions[npts] ]
    smallest_S_family = set_cover_subroutine(partitions_list, verbose=False)
    assert len(smallest_S_family) == size_S_family

    # Compute the TSP
    tsp_cost = mip(G, subtour_callbacks=True, verbose=False)

    # check that k* <= size_S_family
    two_factor_cost_with_S_family, _ = mip(G, initial_subtours=S_family, verbose=False, return_edges=True)
    assert round(two_factor_cost_with_S_family) == round(tsp_cost)

# Save the results as csv file
df.to_csv("./output/tsplib_instances_minimalization.csv")


# ### Without minimalization

minim = False

print(" *** Running the ialg algorithm with NO minimalization on the TSPLIB instances ***")

# Create a dataframe with keys tsplib instances
df = pd.DataFrame(index = [x[0] for x in large_tsplib_instance], columns=["S_min", "b_prime", "runtime", "BeB_n_nodes"])
for instance_name, n in tqdm(large_tsplib_instance):
    # Parse the instance
    G = from_tsplib_file_to_graph("./data/" + instance_name)
    print("Instance:", instance_name)
    (S_family, size_S_family, partitions, c, runtime, num_nodes) = ialg(G, verbose=False, minimalize=minim)
    df.loc[instance_name] = [size_S_family, c, runtime, num_nodes]

    if size_S_family == -1:
        print("Ran into time limit.")
        continue
    else:
        print(
            "\t After {} seconds and {} BeB nodes, we found {} SECs. Specifically, they are:".format(runtime, num_nodes,
                                                                                                     size_S_family))
        for S in S_family:
            print("\t\t", S)

    # check that k* >= size_S_family
    partitions_list = [ [ list(part) for part in partition ] for npts in partitions for partition in partitions[npts] ]
    smallest_S_family = set_cover_subroutine(partitions_list, verbose=False)
    assert len(smallest_S_family) == size_S_family

    # Compute the TSP
    tsp_cost = mip(G, subtour_callbacks=True, verbose=False)

    # check that k* <= size_S_family
    two_factor_cost_with_S_family, _ = mip(G, initial_subtours=S_family, verbose=False, return_edges=True)
    assert round(two_factor_cost_with_S_family) == round(tsp_cost)


# Save the results as csv file
df.to_csv("./output/tsplib_instances_no_minimalization.csv")

'''
HardTSPLIB derived from random instances
'''
# In[ ]:


#hardtsplib_instances_random = [("10001_hard", 10), ("10007_hard", 10), ("10008_hard", 10), ("10010_hard", 10), ("11675_hard", 11), ("12290_hard", 12), ("14850_hard", 14), ("15002_hard", 15), ("15005_hard", 15), ("15007_hard", 15), ("16038_hard", 16), ("20004_hard", 20), ("20007_hard", 20)]

hardtsplib_instances_random_large = [ ("20009_hard", 20), ("20181_hard", 20), ("25001_hard", 25), ("25004_hard", 25), ("25006_hard", 25), ("30001_hard", 30), ("30003_hard", 30), ("30005_hard", 30), ("33001_hard", 33), ("35002_hard", 35), ("35003_hard", 35), ("35009_hard", 35)]


minim = True


# Create a dataframe with keys tsplib instances
df = pd.DataFrame(index = [x[0] for x in hardtsplib_instances_random_large], columns=["S_min", "b_prime", "runtime", "BeB_n_nodes"])
for instance_name, n in tqdm(hardtsplib_instances_random_large):
    # Parse the instance
    G = from_tsplib_file_to_graph("./data/" + instance_name)
    print("******* Instance:", instance_name, "*******")
    (S_family, size_S_family, partitions, c, runtime, num_nodes) = ialg(G, verbose=False, minimalize=minim, smart_initialize=True)
    df.loc[instance_name] = [size_S_family, c, runtime, num_nodes]

    if size_S_family == -1:
        print("Ran into time limit.")
        continue
    else:
        print(
            "\t After {} seconds and {} BeB nodes, we found {} SECs. Specifically, they are:".format(runtime, num_nodes,
                                                                                                     size_S_family))
        for S in S_family:
            print("\t\t", S)

    # check that k* >= size_S_family
    partitions_list = [ [ list(part) for part in partition ] for npts in partitions for partition in partitions[npts] ]
    smallest_S_family = set_cover_subroutine(partitions_list, verbose=False)
    assert len(smallest_S_family) == size_S_family

    # Compute the TSP
    tsp_cost = mip(G, subtour_callbacks=True, verbose=False)

    # check that k* <= size_S_family
    two_factor_cost_with_S_family, _ = mip(G, initial_subtours=S_family, verbose=False, return_edges=True)
    assert round(two_factor_cost_with_S_family) == round(tsp_cost)

# Save the results as csv file
df.to_csv("./output/hardtsplib_instances_random_large_minimalization.csv")

