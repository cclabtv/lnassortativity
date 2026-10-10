import networkx as nx
import random
import numpy as np
import matplotlib.pyplot as plt
import time



def generate(n, d, cap = None):
    '''Generate a classical barabasi albert in linear time in the number of edges. Note that the linear time algoritm may generate multiedges or self loops'''
    
    #The array of size 2nd where d is the number of neighbours chosen by the node. Used in the linear time algorithm
    M = np.empty(2*n*d, dtype=int)
    degree = np.zeros(n, dtype=int)
    graph = nx.MultiGraph()

    #create a clique with duplicated edges of d+1 nodes as a start
    index = 0
    for v in range(d+1):
        degree[v] += 2*d
        for i in range(d+1):
            if(v!=i):
                M[2*index] = v
                M[2*index + 1] = i
                index += 1

    for v in range(d+1, n):
        degree[v] += d
        for i in range(d):
            M[2*(v*d + i)] = v
            r = random.randint(0, 2*(v*d + i))

            while(M[r] == v or (cap != None and degree[M[r]] >= cap)):
                r = random.randint(0, 2*(v*d + i))

            degree[M[r]] += 1
            M[2*(v*d + i) + 1] = M[r]

    for i in range(n*d):
        graph.add_edge(M[2*i], M[2*i + 1])

    return graph

def choice(graph, u, d, p_vector):
    '''choose new neighbours for the given node according to a probability vector'''
    choices = random.choices(graph.nodes, weights=p_vector, k=d)
    for v in choices:
        graph.add_edge(u, v)


def generate_non_lin(n, d, exp=1):
    '''generate a variation of the barabasi-albert graph where the preferential attachment rule is non linear on the degree where a is the exponent on the degree'''

    graph = nx.MultiGraph()
    degree = np.zeros(n, dtype=int)

    #create a clique with duplicated edges of d+1 nodes as a start
    index = 0
    for v in range(d+1):
        degree[v] += 2*d
        for i in range(d+1):
            if(v!=i):
                graph.add_edge(v, i)
                index += 1
    
    for v in range(d+1, n):
        
        p = degree[:v]**exp
        p = p / np.sum(p)

        choices = np.random.choice(graph.nodes, size=d, p=p)
        degree[v] += d
        for u in choices:
            graph.add_edge(v, u)
            degree[u] += 1

    return graph


if __name__ == '__main__':

    graph = generate_non_lin(100, 2, exp=2)
    nx.draw(graph)
    plt.show()