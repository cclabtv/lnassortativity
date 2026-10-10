## Notes on the simulations on the dynamic network model inspired by the Barabasi-Albert graph

We experimented with various possible function for the probability of an edge on the dynamic network to disappear.
This is what we have tried so far:

Given an edge (u, v) we tried:

- $\frac{1}{|d_u - d_v|}$ even normalized with the minimum degree $\frac{1}{|d_u - d_v| \cdot min\{d_u, d_v\}}$ both of these approach (especially the first one) resulted in the network assortativity converging rapidly to -1. 

- $e^{-\beta \cdot |d_u - d_v|}$ this resulted in a network that seems to stabilize around a negative constant value depending on the value of $\Beta$. In particulare for $\beta \sim 1$ the network seems to stabilize around an assortativity value of $r \sim 0.2$.
#
However both of these approach, despite being able to skew the assortativity towards negative values, they do not faithfully represent the behaviour of the lightning network. In fact from the data it seems that degree difference is not correlated (it's even mildly negative correlated) with the lifespan of a channel.
#
Some other idea that may depends on the degree could be one that is inversely proportional to the degree only one of the endpoint. 
Preliminary experiment using the minimum of the two gave negative result, while experiment using the maximum seems to show positive outcomes. This latter experiment requires to be validated by actual statistical data on the LN.

