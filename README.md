# Subgraph optimiser

This is the code that implements a heuristic subgraph optimising algorithm as
described in &lt;yet to be published paper by Lee, Meeks and Pettersson&gt;.

This work was (and is) funded by [EPSRC Grant
EP/T004878/1](https://gow.epsrc.ukri.org/NGBOViewGrant.aspx?GrantRef=EP/T004878/1)


## Background

The goal is to find a subgraph of a given weighted graph that maximises a known
function. Roughly speaking, the value of the function benefits from the graph
containing additional edges, but there is also a penalty associated with
vertices that is linked to the difference between the weight of the vertex and
the average weight of the vertex's neighbour.

This particular code finds a approximate solution by either adding or
removing suitable edges that would increase the value of the function.

The code takes as input a graph (as an n x n adjacency matrix, with each row
beginning on a new line and each line having n entries separated by commas),
and a data file (containing n numbers, one per line). Utility functions to read
such files are given in the `util` package as `read_matrix` and `read_data`
respectively.

The optimiser is created by passing the adjacency matrix and data to
`Optimiser()`. The actual optimisation is then performed by calling
`iterative_opt` on this object. As the name suggests, this function will
iteratively try to either add or remove edges until it cannot find new changes
to make. There are three particular parameters which control whether the code
tries to add or remove edges (or both), as well as a starting point.

* `remove` - If `True`, the algorithm will attempt to remove edges from the
	graph at each iteration.
* `add` - If `True`, the algorithm will attempt to add edges to the graph at
	each iteration. Note that only edges that existed in the original graph will
	ever be added back.
* `remove_first` - If `True`, the algorithm will begin with the original graph,
	and first try to remove edges. Otherwise, the algorithm will begin with an
	empty graph and first try to add edges.

Note that as the function is not well defined for any graph with isolated
vertices, if `remove_first` is False then a further heuristic is run that
finds edges in the original graph whose endpoints have a small difference in
weights to use as a starting point.

## Sample usage


```python
# Read the data,  one number (float) per line, so line 0 has the data for vertex 0 and so-on
data = read_data("data.txt")

# Read an n x n adjacency matrix, each row on a new line and separated by commas
W = read_matrix("adj_matrix.txt")

# Create an optimser object
opt = Optimiser(W, data)

# Call the optimiser with the selected parameters, returning the optimal graph
optimal_graph = opt.iterative_opt(remove=True, add=True, remove_first=False)
```
