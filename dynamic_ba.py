import networkx as nx
import random
import numpy as np
import matplotlib.pyplot as plt
import time
from datetime import datetime
import math


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


def update_nodes(graph, nodes, d):
    ''' Add new d edges to nodes 'node' in graph 'graph' according to the BA model'''
    edge_lst = list(graph.edges())
    for node in nodes:
        new_neigh = []
        for _ in range(d):
           edge = random.choice(edge_lst)
           new_neigh.append(random.choice(edge))
        for u in new_neigh:
            graph.add_edge(node,u)

def get_updating_nodes(graph, rate_act_vec):
    ''' The list of nodes chosen for the update'''
    return [node for node in graph.nodes if random.random() < rate_act_vec[node]]

def remove_edges(graph, death_rate_dict):
    ''' Removes from graph edges according to the death rate'''
    for edge in death_rate_dict.keys():
        if random.random() < death_rate_dict[edge]:
            graph.remove_edge(*edge)

def get_rate_act_vec(graph):
    ''' Returns a rate activation vector for nodes'''
    n = len(graph.nodes)
    ne = len(graph.edges())
    p = 50/ne # one pecent of total number of edges are new
    return [p for _ in range(n)]

def get_death_rate_dict(graph, sel = 'unif', beta = 0):
    ''' Returns a death rate dictionary for edges'''
    if sel == 'unif':
        ne = len(graph.edges())
        p = 0.0015
        out = {edge: p for edge in graph.edges()}
    if sel == 'deg_dist':
        out = dict()
        for edge in graph.edges():
            edge_degrees = (graph.degree[edge[0]], graph.degree[edge[1]])
            m, M = min(edge_degrees), max(edge_degrees)
            out[edge] = 1.0/(((M-m)+1))
    if sel == 'exp_deg':
        out = dict()
        for edge in graph.edges():
            edge_degrees = (graph.degree[edge[0]], graph.degree[edge[1]])
            m, M = min(edge_degrees), max(edge_degrees)
            out[edge] = math.exp(-beta*(M-m))
    return out


if __name__ == '__main__':
    
    d = 5
    rounds = 10000
    beta = 2.5
    n = 100
    graph = generate(n, d)
    text = []
    for i in range(rounds):
        death_rate_dict = get_death_rate_dict(graph, sel = 'exp_deg', beta = beta)
        rate_act_vec = get_rate_act_vec(graph)
        updating_nodes = get_updating_nodes(graph, rate_act_vec)
        updating_nodes.append(max(graph.nodes) + 1)
        update_nodes(graph, updating_nodes, d)
        remove_edges(graph, death_rate_dict)
        if i % 100 == 0:
            text.append(str(len(graph.nodes))+','+str(len(graph.edges()))+','+str(nx.degree_assortativity_coefficient(graph))+'\n')
    unix_timestamp = int((datetime.utcnow() - datetime(1970, 1, 1)).total_seconds())
    header = '# '+ str(rate_act_vec) + '\n# '+ str(death_rate_dict) + '\n'
    fname = 'assortativity_' + 'n_'+str(n) + '_d_' + str(d)+'_rounds_' + str(rounds) + '_beta_'+ str(beta)+'_' + str(unix_timestamp)+ '.csv'
    with open(fname, 'a') as f:
        f.write(header + ''.join(text))

    # nx.draw(graph)
    # plt.show()
