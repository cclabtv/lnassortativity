import networkit as nk
from networkit.correlation import Assortativity
from networkit.generators import BarabasiAlbertGenerator
from tqdm import tqdm
import matplotlib.pyplot as plt
from numpy import mean

ass1 = []
ass2 = []
ass3 = []

for i in tqdm(range(1000, 1000000, 1000)):
    temp1 = []
    temp2 = []
    temp3 = []
    gen1 = BarabasiAlbertGenerator(1, i)
    gen2 = BarabasiAlbertGenerator(2, i)
    gen3 = BarabasiAlbertGenerator(3, i)
    for j in range(10):
        G1 = gen1.generate()
        deg = [G1.degree(u) for u in G1.iterNodes()]
        temp1.append(Assortativity(G1, deg).getCoefficient())
        G2 = gen2.generate()
        deg = [G2.degree(u) for u in G2.iterNodes()]
        temp2.append(Assortativity(G2, deg).getCoefficient())
        G3 = gen3.generate()
        deg = [G3.degree(u) for u in G3.iterNodes()]
        temp3.append(Assortativity(G3, deg).getCoefficient())
    
    ass1.append(mean(temp1))
    ass2.append(mean(temp2))
    ass3.append(mean(temp3))

with open('ass1.txt', 'w') as f:
    for item in ass1:
        f.write(f'{item}\n')

with open('ass2.txt', 'w') as f:
    for item in ass2:
        f.write(f'{item}\n')

with open('ass3.txt', 'w') as f:
    for item in ass3:
        f.write(f'{item}\n')

# Create a list of x values (number of nodes)
y = list(range(1000, 1000000, 1000))

# Create a new figure
plt.figure()

# Plot each list with a label for the corresponding alpha value
plt.plot(ass1, y, label="m=1")
plt.plot(ass2, y, label="m=2")
plt.plot(ass3, y, label="m=3")

# Add a legend
plt.legend()

# Add title and labels for the axes
plt.title('Assortativity Coefficient of Barabasi-Albert')
plt.ylabel('Number of Nodes')
plt.xlabel('Assortativity Coefficient')

# Display the plot
plt.show()
plt.savefig("Assortativity_BA.png")
