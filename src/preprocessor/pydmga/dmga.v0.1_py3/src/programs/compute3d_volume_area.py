from pydmga.geometry import OrthogonalGeometry
from pydmga.container import Container
from pydmga.diagram import Diagram
from pydmga.io import pdb
from math import sqrt

# this is basic filename of the dataset, withou .pdb extension
filename = "drec-1_590-600ns"
# use heavy atoms only - if True then we skip all hydrogen atoms
hao = True
# old lipid stops at what residue?
end_lipid_at = "UDC"
# new lipid starts at what residue?
start_lipid_at = "KDO"

# we will read from <filename>.pdb
in_file = file("{0}.pdb".format(filename), "r")
# we will write to <filename>.dat
out_file = file("{0}.dat".format(filename), "w")

try:
    # this counts current frame, for output
    frame = 0
    # we collect volume and area per lipid per frame
    out_file.write("frame lipid volume area\n")
    while True:  # iterate over all frames
        # here we remember data from pdb for future use
        particles = []
        while True:  # iterate over elements in frame
            item = pdb.PDBRecord(pdb.readnext(in_file))
            if item.type() == "CRYST1":
                # create box geometry periodic in all directions
                geometry = OrthogonalGeometry(
                    item["a"], item["b"], item["c"], True, True, True)
                # create container for atoms
                container = Container(geometry)
            if item.type() == "ATOM":
                # add new item, if we are using only hao, then skip it if its hydrogen
                if (item["atom"][0] != 'H' and item["atom"][1] != 'H') or not hao:
                    # add perticle to the container
                    container.add(item.as_particle())
                    # remember all data for future use
                    particles.append(item)
            if item.type() == "ENDMDL" or item.type() == "TER":
                # stop reading this frame
                break
        # create new Power Diagram
        diagram = Diagram(container)
        # here we will count which lipid molecule we have
        molecule_no = 0
        # this is to allow triggering start_lipid_at test for the first iteration
        residue = end_lipid_at
        # iterate over all cells in this diagram
        for i, cell in enumerate(diagram):
            # if we found one solvent then skip the rest of atoms - they are solvent
            if particles[i]["residue"] == "SOL" or particles[i]["residue"] == "ION":
                break
            # we change lipid molecule only if we pass from <end_lipid_at> residue to <start_lipid_at> residue
            if residue == end_lipid_at and particles[i]["residue"] == start_lipid_at:
                # write data to a file if there is a data
                if molecule_no:
                    out_file.write("{0} {1} {2} {3}\n".format(
                        frame, molecule_no, volume, area))
                molecule_no += 1
                area = 0.0
                volume = 0.0
            # to track when to change residue
            residue = particles[i]["residue"]
            # find cell volume
            volume += cell.volume()
            # iterate over all sides to find contact area for this atom
            for side in cell.sides:
                neighbour = side.neighbour
                if particles[neighbour]['residue'] == "SOL" or particles[neighbour]['residue'] == "ION":
                    # we have solvent <-> particle contact
                    area += side.area()
        # need to write last lipid (we changed frames, not lpipid inside the frame)
        out_file.write("{0} {1} {2} {3}\n".format(
            frame, molecule_no, volume, area))
        # count new frame
        frame += 1
except StopIteration:
    # pdb file ends
    pass
