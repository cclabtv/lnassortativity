import networkx as nx
import helpers
from math import exp

class Net(nx.Graph):

    def __init__(self, n0 = 5):
        ''' Inheriths all methos of the class networkx.Graph and initializes a
        clique with n0 nodes '''
        super().__init__()
        for u in range(n0):
            self.add_node(u)
            for v in range(u):
                self.add_edge(u,v, capacity = 1)
        self.largest_node_id = n0-1

    def degree_seq(self):
        ''' Returns a pair of arrays (node, degree): "degree[i]" is the
        degree of node "node[i]" '''
        node = list(self.nodes)
        degree = [self.degree[node[i]] for i in range(len(node))]
        return node, degree

    def capacity_seq(self):
        ''' Returns a pair of arrays (node, capacity): capacity[i] is the sum
        of the capacities of the edges incident on "node[i]" '''
        node = list(self.nodes)
        capacity = [sum([self[edge[0]][edge[1]]['capacity'] for edge in \
                      self.edges(node[i])]) for i in range(len(node))]
        return node, capacity

    def node_capacity(self, node):
        ''' Returns the sum of the capacities of edges incident on "node"''' 
        return sum([self[edge[0]][edge[1]]['capacity'] for edge in \
                      self.edges(node)])
    
    def edge_endpoints_capacity(self, edge):
        ''' Returns the sum of the capacities of edges incident on the endpoints
        of "edge" (capacity of "edge" counted twice)''' 
        return self.node_capacity(edge[0]) + self.node_capacity(edge[1]) 
 
    def assortativity(self, weight = None):
        ''' Returns the (possibly weighted) assortativity '''
        return nx.degree_assortativity_coefficient(self, weight = weight)

    def select_node_weight(self, with_capacity = False):
        ''' Chooses a node with probability proportional to the node degree
        or to the node capacity (sum of the capacities of incident edges) '''
        if with_capacity:
            return helpers.rnd_choice_weight(*self.capacity_seq())
        return helpers.rnd_choice_weight(*self.degree_seq())

    def select_nodes_fixed_rate(self, rate):
        ''' Chooses each node with probability "rate" '''
        node = list(self.nodes)
        rates = [rate for _ in range(len(node))]
        return helpers.rnd_choice_rate(node, rates)

    def select_nodes_degree_prop_rate(self, exponent = 2):
        ''' Chooses each node "u" with probability proportional to
        1 / deg(u)^exponent '''
        node = list(self.nodes)
        rates = [1 if self.degree(node[i]) == 0 \
                   else 1.0 / self.degree(node[i])**exponent for i in range(len(node))]
        return helpers.rnd_choice_rate(node, rates)

    def select_nodes_degree_prop_rate_beta(self, beta = 1):
        ''' Chooses each node "u" with probability proportional to
        exp{-beta*deg(u)} '''
        node = list(self.nodes)
        rates = [exp(-beta*self.degree(node[i])) for i in range(len(node))]
        return helpers.rnd_choice_rate(node, rates)

    def select_edges_fixed_rate(self, rate):
        ''' Chooses each edge with probability "rate" '''
        edge = list(self.edges)
        rates = [rate for _ in range(len(edge))]
        return helpers.rnd_choice_rate(edge, rates)

    def select_edges_capacity_prop_rate(self, exponent, local = True):
        ''' Chooses each edge "e" with probability proportional to
        capacity of e divided by the sum of the capacities of endponts'''
        edge = list(self.edges)
        rates = [(self[edge[i][0]][edge[i][1]]['capacity'] / \
                 self.edge_endpoints_capacity(edge[i]))**exponent \
                 for i in range(len(edge))]
        return helpers.rnd_choice_rate(edge, rates)

    def select_edges_capacity_prop_rate_beta(self, beta, local = True):
        ''' Chooses each edge "e = {u,v}" with probability 
        exp(beta * cap(e)) / [exp(beta cap(u)) + exp(beta cap(v))] '''
        edge = list(self.edges)
        rates = [1 / (exp(beta*(self.node_capacity(edge[i][0]) - \
                                self[edge[i][0]][edge[i][1]]['capacity'])) + \
                      exp(beta*(self.node_capacity(edge[i][1]) - \
                                self[edge[i][0]][edge[i][1]]['capacity']))) \
                 for i in range(len(edge))]
        return helpers.rnd_choice_rate(edge, rates)

    def add_one_node(self):
        node = self.largest_node_id + 1
        self.add_node(node)
        self.largest_node_id = node
        return self.largest_node_id

    def add_one_edge(self, edge, capacity):
        u, v = edge
        self.add_edge(u,v, capacity = capacity)

