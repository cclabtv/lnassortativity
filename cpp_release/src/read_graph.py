#import pandas as pd
import networkx as nx


def read_csv_graph(filename):
    '''
    Simple python scripts to convert the generated graph through the c++ code in a networkx multigraph

    '''
    
    f = open(filename, 'r')
    lines = f.readlines()
    g = nx.MultiGraph()

    for line in lines[1:]:
        if line.startswith('n'):
            continue
        #divide the two lists
        x = line.split(']')
        l1 = x[0]+']'
        l2 = x[1][1:]+']'

        #parse the string into two lists of integers
        res1 = [st + ')' for st in l1.strip('][').split('),')]
        res1[-1] = res1[-1][:-1]
        
        newEdges = [(int(st.strip(')(').split(',')[0]), int(st.strip(')(').split(',')[1])) for st in res1]
        
        res2 = [st + ')' for st in l2.strip('][').split('),')]
        res2[-1] = res2[-1][:-1]
        
        removedEdges = [(int(st.strip(')(').split(',')[0]), int(st.strip(')(').split(',')[1])) for st in res2 if len(st)!=0]
        
        for (u,v) in newEdges:
            g.add_edge(u, v)

        
        for (u,v) in removedEdges:
            if (u, v) in g.edges():
                g.remove_edge(u, v)
    
    return g

if __name__ == "__main__":

    '''
    Example on how to read the graph using this script
    '''
    g = read_csv_graph('historycalGraph1712397468.csv')

    print(nx.degree_assortativity_coefficient(g))






#print(nx.degree_assortativity_coefficient(g))