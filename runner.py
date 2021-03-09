#!/usr/bin/env python3

from utils import read_matrix, read_data
from optimiser import Optimiser

# Read the data,  one number (float) per line, so line 0 has the data for vertex 0 and so-on
data = read_data("data.txt")

# Read an n x n adjacency matrix, each row on a new line and separated by commas
W = read_matrix("adj_matrix.txt")

# Create an optimser object
opt = Optimiser(W, data)

# Call the optimiser with the selected parameters, returning the optimal graph
optimal_graph = opt.iterative_opt(remove=True, add=True, remove_first=False)
