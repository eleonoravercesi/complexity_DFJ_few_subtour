# The Dantzig-Fulkerson-Johnson TSP formulation is easy to solve for few subtour constraints

This is GitHub repository accompaining the paper "The Dantzig-Fulkerson-Johnson TSP formulation is easy to solve for few subtour constraints".

The folder contains the following files:

- `tsplib_hardtsplib.ipynb`  for reproducing the results of Table 1, 3 and 4.
- `focus_dantzig42_brazil58.ipynb` for the results of Table 2 and Figure 3.
- `random.ipynb` for the results of Figure 4 and 5. 
- `plot.ipynb` for the results of Figure 6.
- `rectilinear_3D_instances.ipynb` for the results of Section 5.4.
- `utils.py, ialg.py, cover.py` containing the functions used in the notebooks, that can be of independent interest.

To run the notebooks, first put in the folder `data` and put there the instances of the [TSPLIB](http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/) [1] and the [HardTSPLIB](https://github.com/eleonoravercesi/HardTSPLIB/) [2].

[1] Reinelt, Gerhard. "TSPLIB—A traveling salesman problem library." ORSA journal on computing 3.4 (1991): 376-384.


[2] Vercesi, Eleonora, et al. "On the generation of metric TSP instances with a large integrality gap by branch-and-cut." Mathematical Programming Computation 15.2 (2023): 389-416.