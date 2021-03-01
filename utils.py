"""Util functions."""

try:
    from math import log as ln, isqrt
except ImportError:  # python3.8 is needed for isqrt
    from math import log as ln
    import operator
    ## Following taken from comments of python3.8 code.
    def isqrt(n):
        """
        Return the integer part of the square root of the input.
        """
        n = operator.index(n)
        if n < 0:
            raise ValueError("isqrt() argument must be nonnegative")
        if n == 0:
            return 0
        c = (n.bit_length() - 1) // 2
        a = 1
        d = 0
        for s in reversed(range(c.bit_length())):
            # Loop invariant: (a-1)**2 < (n >> 2*(c - d)) < (a+1)**2
            e = d
            d = c >> s
            a = (a << d - e - 1) + (n >> 2*c - e - d + 1) // a
        return a - (a*a > n)


def powerset(seq):
    """ generator returning the powerset of a set """
    if len(seq) == 0:
        yield []
    else:
        for item in powerset(seq[1:]):
            yield [seq[0]]+item
            yield item




class Graph:

    def __init__(self, matrix):
        self._adj = matrix
        self._degrees = [sum(row) for row in matrix]
        self._nbs = [set(w for w, val in enumerate(row) if val == 1) for row in self._adj]

    def size(self):
        return len(self._adj)

    def vertices(self):
        return range(self.size())

    def add_edges(self, edges):
        for e in edges:
            self.add_edge(e)

    def add_edge(self, edge):
        self._adj[edge[0]][edge[1]] = 1
        self._adj[edge[1]][edge[0]] = 1
        self._degrees[edge[0]] += 1
        self._degrees[edge[1]] += 1
        self._nbs[edge[0]].add(edge[1])
        self._nbs[edge[1]].add(edge[0])

    def remove_edges(self, edges):
        for e in edges:
            self.remove_edge(e)

    def remove_edge(self, edge):
        self._adj[edge[0]][edge[1]] = 0
        self._adj[edge[1]][edge[0]] = 0
        self._degrees[edge[0]] -= 1
        self._degrees[edge[1]] -= 1
        self._nbs[edge[0]].remove(edge[1])
        self._nbs[edge[1]].remove(edge[0])

    def degree(self, v):
        return self._degrees[v]

    def nbs(self, v):
        return self._nbs[v]

    def edges(self):
        res = []
        for v, row in enumerate(self._adj):
            for w, entry in enumerate(row):
                # Only consider edges
                if entry == 0:
                    continue
                # Only consider (v,w) if v < w
                if w <= v:
                    continue
                res.append([v,w])
        return res

    def hasEdge(self, edge):
        return self._adj[edge[0]][edge[1]] == 1

    def isSubgraphOf(self, other):
        for edge in self.edges():
            if not other.hasEdge(edge):
                return False
        return True


class EmptyGraph(Graph):
    def __init__(self, adj_matrix):
        self._possible = Graph(adj_matrix)
        empty = []
        for i in range(self._possible.size()):
            empty.append( [0]*self._possible.size())
        super().__init__(empty)

    def possNbs(self, v):
        return self._possible.nbs(v) - self.nbs(v)

    def possibleDegree(self, v):
        return self._possible.degree(v)

    def allPossNbs(self, v):
        return self._possible.nbs(v)

    def origEdges(self):
        return self._possible.edges()


def read_data(filename):
    """Read data from a file. One line per value, all as floats."""
    with open(filename, "r") as infile:
        return [float(line.rstrip()) for line in infile]


def read_matrix(filename):
    """Read a square matrix from comma-separated file.
    Note that we ignore markers like line-endings.
    The matrix is returned as a list of lists.
    """
    entries = []
    with open(filename, "r") as infile:
        for line in infile:
            entries.extend( int(x.rstrip()) for x in line.rstrip().split(",") if x.rstrip())
    size = isqrt(len(entries))
    assert(size ** 2 == len(entries))
    matrix = []
    while entries:
        matrix.append(entries[:size])
        entries = entries[size:]
    return matrix


def write_sol(filename, graph):
    with open(filename, "w") as outfile:
        for row in graph._adj:
            outfile.write(",".join(str(e) for e in row))
            outfile.write("\n")


def adj_cont(v, matrix, phi, other_matrix):
    first = ln(degree(v, matrix))/2
    inner = 1 + (degree(v, matrix)*neighbour_discrepancy(v, matrix, phi)) / (sum(degree(w, other_matrix)*neighbour_discrepancy(w, other_matrix, phi) for w in range(len(matrix))) - degree(v, matrix)*neighbour_discrepancy(v, matrix, phi))
    return first - len(matrix)/2 * ln (inner)


def neighbour_discrepancy(vertex, graph, phi):
    """Calculate the neighbourhood discrepancy of a matrix in a graph.
    """
    return (phi[vertex] - sum( phi[n] for n in graph.nbs(vertex))/graph.degree(vertex))**2


def func_value(graph, phi):
    """Given a Graph G, and a list phi representing the function \phi : V(G) \mapsto \mathbb{R}
    calculate J_H(\phi)
    """
    K = graph.size()
    val = 1/2 * sum(ln(graph.degree(v)) for v in range(0, K))
    sum_nd = sum(graph.degree(v) * neighbour_discrepancy(v, graph, phi) for v in range(0, K))
    assert(sum_nd != 0)
    val -= K/2 * ln(sum_nd)
    return val


def ordering(matrix):
    """Given a full adjacency matrix, calculates a degeneracy ordering using Matula and Beck (1983)
    """
    L = []
    D =[[] for _ in range(len(matrix))]
    nbrs = []
    dv =[sum(val for val in row) for row in matrix]
    for ind, row in enumerate(matrix):
        deg = sum(val for val in row)
        D[deg].append(ind)
        nbrs.append([other for other, val in enumerate(row) if val])
    k = 0
    for rep in range(len(matrix)):
        ind = 0
        for Dcell in D:
            if Dcell:
                break
            ind += 1
        k = max(k, ind)
        L.insert(0, D[ind].pop(0))
        for w in nbrs[L[0]]:
            if w not in L:
                assert(w in D[dv[w]])
                D[dv[w]].remove(w)
                dv[w] -= 1
                D[dv[w]].append(w)
    return L, k

