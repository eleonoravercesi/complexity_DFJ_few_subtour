from utils import from_tsplib_file_to_graph
from ialg import ialg
import pandas as pd
import numpy as np

instance_name = "dantzig42.tsp"
#instance_name = "brazil58.tsp"

G = from_tsplib_file_to_graph("./data/" + instance_name)
print("******* Instance:", instance_name, "*******")

df = pd.DataFrame(0, index=np.arange(4), columns=["minimalize", "smart_init", "S_min", "b", "runtime", "bb_nodes"])

# Redirect the standard output on a log file
import sys
sys.stdout = open("log_" + instance_name + ".txt", "w+")
cont = 0
for minimalize in [True, False]:
    for smart_init in [True, False]:
        (S_family, size_S_family, partitions, c, runtime, num_nodes) = ialg(G, verbose=True, minimalize=minimalize, smart_initialize=smart_init)
        # Add a row to df
        df.iloc[cont] = [minimalize,  smart_init,  size_S_family, c, runtime, num_nodes]
        cont += 1
        print(" ")
        print("----------------------------")
        print(" ")
        # Save df to csv
        df.to_csv("../why_tsp_easy/csv_files/results_" + instance_name + ".csv")

