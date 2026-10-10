import barabasi
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import pickle
import os
import time
from tqdm import tqdm

def experiment(n_of_nodes, n_of_trials, graph_model, *args, **kwargs):
    '''Return the average degree assortativity of a graph of size n_of_nodes on n_of_trials trials'''
    assortativity = []
    for i in range(n_of_trials):

        graph = graph_model(n_of_nodes, *args, **kwargs)
        assortativity.append(nx.degree_assortativity_coefficient(graph))
    
    return np.mean(assortativity)


def simulation(n_of_trials, graph_model, start, end, step, save_every_n_minutes=5, progress_file='progress.pkl', *args, **kwargs):
    '''draw a graph on the evolution of the degree assortativity for increasing value of n given a specific graph model'''
    x = list(range(start, end, step))
    y = []
    params = (n_of_trials, graph_model, start, end, step, save_every_n_minutes, args, kwargs)
    
    last_save_time = time.time()

    for i in tqdm(x):  
        y.append(experiment(i, n_of_trials, graph_model, *args, **kwargs))

        # Save progress to a file every n minutes
        current_time = time.time()
        if current_time - last_save_time >= save_every_n_minutes*60:
            with open(progress_file, 'wb') as f:
                pickle.dump((params, x, y), f)
            last_save_time = current_time

    plt.suptitle(str(graph_model.__name__) + "parameters: " +  str(args) + " " + str(kwargs))
    plt.xlabel("Number of nodes")
    plt.ylabel("Degree assortativity")
    plt.plot(x, y)
    plt.savefig(str(graph_model.__name__) + "parameters:" +  str(args) + str(kwargs))
    plt.show()

def restore_simulation(progress_file='progress.pkl'):
    '''Restore the simulation from a pickle file'''

    # Check if a saved progress file exists
    if os.path.exists(progress_file):
        # Load progress from file
        with open(progress_file, 'rb') as f:
            params, x, y = pickle.load(f)
        
        (n_of_trials, graph_model, start, end, step, save_every_n_minutes, args, kwargs) = params

        last_save_time = time.time()

        for i in tqdm(x[len(y):]):  # Continue from where we left off
            y.append(experiment(i, n_of_trials, graph_model, *args, **kwargs))

            # Save progress to a file every n minutes
            current_time = time.time()
            if current_time - last_save_time >= save_every_n_minutes*60:
                with open(progress_file, 'wb') as f:
                    pickle.dump((params, x, y), f)
                last_save_time = current_time
            
        plt.suptitle(str(graph_model.__name__) + "parameters: " +  str(args) + " " + str(kwargs))
        plt.xlabel("Number of nodes")
        plt.ylabel("Degree assortativity")
        plt.plot(x, y)
        plt.savefig(str(graph_model.__name__) + "parameters:" +  str(args) + str(kwargs))
        plt.show()
    else:
        print("No progress file found.")



if __name__ == '__main__':
    simulation(20, barabasi.generate_non_lin, 10000, 100000, 10000, 5, 'progress.pkl', 1, 1.5)
    #restore_simulation('progress.pkl')