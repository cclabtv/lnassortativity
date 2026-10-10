from net import Net
import helpers

def main(days = 2305, \
         avg_nodes_per_day = 10, \
         new_channels_num = 5, \
         fname = 'model_data.csv'):
    G = Net()
    text = ''
    for day in range(days):
        new_nodes_num = helpers.rnd_geometric(1.0 / avg_nodes_per_day)
        for _ in range(new_nodes_num):
            new_node = G.add_one_node()
            for _ in range(new_channels_num):
                neighbor = G.select_node_weight(with_capacity = True)
                capacity = helpers.rnd_binomial(10**8, .01)
                G.add_one_edge((new_node, neighbor), capacity)
        disappearing_edges = G.select_edges_capacity_prop_rate(exponent = 2)
        for e in disappearing_edges:
            G.remove_edge(*e)
        disappearing_nodes = G.select_nodes_degree_prop_rate(exponent = 4)
        for u in disappearing_nodes:
           G.remove_node(u)
        text += str(day)+','+str(G.number_of_nodes()) + ',' \
                + str(G.number_of_edges()) + ',' \
                + str(round(G.assortativity(), 3)) + '\n'
        if day % 100 == 0:
            with open(fname, 'a') as f:
                f.write(text)
            text = ''

if __name__ == '__main__':
    main()
