import settings
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist
from pointGenerators import generateNodesRect

dim = 2
maxlim = np.array([1., 1.])
minDist = 0.07
trials = 5000
node_coords = [[5, 1], [0, 0], [10, 0]]
generateNodesRect(maxlim, minDist, dim, trials, node_coords)

node_coords = np.array(node_coords)
distances = pdist(node_coords, 'euclidean')
print(distances.min(), distances.max())
print(len(node_coords))

fig, ax = plt.subplots()
ax.plot(node_coords[:, 0], node_coords[:, 1], 'bx')
plt.show()
