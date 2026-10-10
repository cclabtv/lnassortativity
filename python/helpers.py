import random
import numpy as np
from scipy.stats import powerlaw

def rnd_choice_weight(elements, weights):
    tot = sum(weights)
    probs = [w / tot for w in weights] if tot != 0\
            else [1.0 / len(weights) for w in weights]
    return np.random.choice(elements, p = probs)

def rnd_choice_rate(elements, rates):
    return [elements[i] for i in range(len(elements)) \
            if random.random() < rates[i]]

def rnd_binomial(n, p, size = None):
    return np.random.binomial(n, p, size)

def rnd_geometric(p, size = None):
    return np.random.geometric(p, size)

def rnd_powerlaw(a, loc = None, scale = None, size = None):
    return int(powerlaw.rvs(a, loc = loc, scale = scale, size = size))


if __name__ == '__main__':
    print(rnd_binomial(10, .7))
    print(rnd_geometric(0.2))
    elements = [1,2,3]
    weights = [0,0,0]
    print(rnd_choice_weight(elements, weights))
