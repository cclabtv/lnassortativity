from net import Net
import helpers
from configparser import ConfigParser
import argparse
from time import time
from os import path,makedirs

def args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config', help = 'Configuration File', \
                        default = 'model.conf')
    parser.add_argument('-o', '--output', help = 'Output File', default = \
                        'model_data.csv')
    parser.add_argument('-e', '--exponent', help = 'Exponent for edge disappearing rate', \
                        default = 2)
    parser.add_argument('-b', '--beta', help = 'Numerator for Beta value', \
                        default = 1.0)
    return parser.parse_args()

config = ConfigParser()
fconf = args().config
config.read(fconf)

DAYS = int(config['Default']['days'])
AVG_NODES_PER_DAY = int(config['Default']['avg_nodes_per_day'])
INITIAL_CLIQUE_SIZE = int(config['Default']['initial_clique_size'])
NEW_CHANNELS_NUM = int(config['Default']['new_channels_num'])
CAPACITY_DISTRIBUTION = config['Default']['capacity_distribution']
CAPACITY_A = int(config['Default']['capacity_a'])
CAPACITY_LOC = int(config['Default']['capacity_loc'])
CAPACITY_SCALE = int(config['Default']['capacity_scale'])
CAPACITY_N = int(config['Default']['capacity_n'])
CAPACITY_P = float(config['Default']['capacity_p'])

fout = args().output
def main(fname = fout, runs = 10):
    makedirs(datadir, exist_ok=True)
    datadir = path.join(path.dirname(__file__),'datadir')
    fname = 'e' + str(args().exponent) + '_' + fname
    for _ in range(runs):
        G = Net(n0 = INITIAL_CLIQUE_SIZE)
        text = ''
        ts_fname = 'ts' + str(int(time())) + '_' + fname 
        ts_fname = path.join(datadir, ts_fname)
        for day in range(DAYS):
            new_nodes_num = helpers.rnd_geometric(1.0 / AVG_NODES_PER_DAY)
            for _ in range(new_nodes_num):
                new_node = G.add_one_node()
                for _ in range(NEW_CHANNELS_NUM):
                    neighbor = G.select_node_weight(with_capacity = True)
                    if CAPACITY_DISTRIBUTION == 'binomial': 
                        capacity = helpers.rnd_binomial(CAPACITY_N, CAPACITY_P)
                    elif CAPACITY_DISTRIBUTION == 'powerlaw':
                        capacity = helpers.rnd_powerlaw(CAPACITY_A, \
                                     loc = CAPACITY_LOC, scale = CAPACITY_SCALE)
                    G.add_one_edge((new_node, neighbor), capacity)
            # num_beta = float(args().beta)
            # disappearing_edges = G.select_edges_capacity_prop_rate_beta(beta = num_beta / CAPACITY_N)
            exponent = float(args().exponent)
            disappearing_edges = G.select_edges_capacity_prop_rate(exponent = exponent)
            for e in disappearing_edges:
                G.remove_edge(*e)
            # disappearing_nodes = G.select_nodes_fixed_rate(rate = 1.0/20000)
            # for u in disappearing_nodes:
            #    G.remove_node(u)
            text += str(day)+','+str(G.number_of_nodes()) + ',' \
                    + str(G.number_of_edges()) + ',' \
                    + str(round(G.assortativity(), 3)) + '\n'
            if day % 100 == 0:
                with open(ts_fname, 'a') as f:
                    f.write(text)
                text = ''

if __name__ == '__main__':
    main()
