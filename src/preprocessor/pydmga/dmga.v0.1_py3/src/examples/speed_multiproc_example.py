import sys
sys.path.append('../')

from pydmga.geometry import OrthogonalGeometry
from pydmga.container import Container
from pydmga.diagram import Diagram
from random import random
from random import seed
import math
import pp
from multiprocessing import Pool
import time
#import sys
#from dmga.io.pdb import PDBReader

(box_x, box_y, box_z) = (300, 300, 300)
print("[TEST] Classic Voronoi test")
g = OrthogonalGeometry(box_x, box_y, box_z, True, True, True)
c = Container(g)
if (len(sys.argv) > 1):
    count = int(sys.argv[1])
else:
    count = 1000

for i in range(count):
    c.add(i, box_x * random(), box_y * random(), box_z * random(), 1.0)
d = Diagram(c)

# make parallel
ppservers = ()
if (len(sys.argv) > 2):
    ncpus = int(sys.argv[2])
else:
    ncpus = 1


def compute_part(params):
    (diagram, start, step, last) = params
    i = start
    v = 0
    while i < last:
        #print "jestem w", start, ",", i
        c = diagram.get_cell(i)
        v += c.volume()
        i += step
    return v


pool = Pool(processes=ncpus)
parts = pool.map(compute_part, [(d, i, ncpus, count) for i in range(ncpus)])
v = sum(parts)

print(v, 'vs', box_x*box_y*box_z)
