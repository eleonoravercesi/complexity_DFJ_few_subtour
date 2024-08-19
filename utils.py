'''
Function used to reproduce the results in "On the complexity of the Dantzig-Fulkerson-Johnson TSP formulation for few subtour constraints"
Authors: E. Vercesi, A. Buchanan
'''

import math
import networkx as nx
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
        D = np.random.exponential(999, (n, n)) # You sample the whole matrix at once, but you only use the upper triangular part
        G = nx.Graph()
        for i in range(n):
            for j in range(i + 1, n):
                G.add_edge(i, j, cost=round(D[i, j]))
        return G


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
            G.add_edge(i, j, cost=int(C[i, j])) # Ensure that the cost is an integer
    return G