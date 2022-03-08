import os
# authors @ abhinav tamaskar and james bannon
import cvxopt as cx
import networkx as nx
import numpy as np
import random
import math
import sys

def compute_tv_vector(model:BayesianNetwork)->np.array:
    """ Extended BN to Extended BN"""
    # print(model.nodes)
    if len(model.nodes) > 20:
        all_events = sample_binary_events(model.nodes)
    else:
        all_events = pd.DataFrame(data=np.array(list(map(list, itertools.product([0, 1],
        repeat=len(model.nodes))))),columns = model.nodes)

    results = []
    for i,r in all_events.iterrows():
        res = get_event_probability(model,r.to_dict())
        results.append(res)
    return np.array(results)

def compute_union_graph(G1:nx.DiGraph,G2:nx.DiGraph)->nx.DiGraph:
    U = nx.DiGraph()
    for node in G1.nodes():
        U.add_node(node)
    for edge in set(G1.edges()).union(set(G2.edges())):
        U.add_edge(edge[0],edge[1])
    return U

def graph_agony(graph:nx.DiGraph)->int:

    # use this as example and to understand the code
    # https://cvxopt.org/userguide/coneprog.html#cvxopt.solvers.lp
    # UnionGraph = ConstructUnionGraph(G1,G2)
    n = graph.number_of_nodes()
    m = graph.number_of_edges()

    # the variables in our case are p(u, v), r(u), r(v) (as in the paper)
    # the first m variables are for p(u,v)
    # the next n variables are for r(u)

    # we now want to minimize the sum of agony of edges only
    c = [1. for i in range(m)] + [0. for i in range(n)]

    # we have 2*(m + n) equations with n+m variables giving a (n+m)x(2*(m+n)) matrix
    G = [[0. for i in range(2 * (m + n))] for j in range(n + m)]

    # the first m equations
    edge_list = list(graph.edges())
    for i, (u, v) in enumerate(edge_list):
        G[m + v][i] = 1.
        G[m + u][i] = -1.
        G[i][i] = -1.
    for i in range(m, 2 * m):
        G[i - m][i] = -1.
    for i in range(2 * m, 2 * m + n):
        G[i - m][i] = 1.
    for i in range(2 * m + n, 2 * m + 2 * n):
        G[i - m - n][i] = -1.

    c = cx.matrix(c)
    G = cx.matrix(G)
    h = cx.matrix([-1. for i in range(m)] + [0. for i in range(m)] +
                  [0. + n for i in range(n)] + [0. for i in range(n)])
    solution = cx.solvers.lp(c, G, h)
    raw_ag = int(round(solution['primal objective']))
    return raw_ag, (1.0/m)*raw_ag
