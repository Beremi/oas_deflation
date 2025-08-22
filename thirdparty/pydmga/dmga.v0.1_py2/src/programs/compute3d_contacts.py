from pydmga.geometry import OrthogonalGeometry
from pydmga.container import Container
from pydmga.diagram import Diagram
from pydmga.io import pdb
from math import sqrt

# this is basic filename of the dataset, withou .pdb extension
filename = "drec-1_590-600ns"
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
	# we collect number of water in contact and number of ions in contact
	out_file.write("frame lipid nb_h2o nb_ion\n")
	while True: # iterate over all frames		
		# here we remember data from pdb for future use
		particles = []		
		while True:	# iterate over elements in frame
			item = pdb.PDBRecord(pdb.readnext(in_file))
			if item.type() == "CRYST1":	
				# create box geometry periodic in all directions
				geometry = OrthogonalGeometry(item["a"], item["b"], item["c"], True, True, True)				
				# create container for atoms
				container = Container(geometry)				
			if item.type() == "ATOM":
				# add new item, for simplicity we use only heavy atoms
				if (item["atom"][0] != 'H' and item["atom"][1] != 'H'):
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
		#iterate over all cells in this diagram				
		for i, cell in enumerate(diagram):
			# break from iteration when SOL/ION part of file reached
			if particles[i]["residue"] == "SOL" or particles[i]["residue"] == "ION": break
			# we change lipid molecule only if we pass from <end_lipid_at> residue to <start_lipid_at> residue
			if residue == end_lipid_at and particles[i]["residue"] == start_lipid_at: 
				# write data to a file if there is a data
				if molecule_no: out_file.write("{0} {1} {2} {3}\n".format(frame, molecule_no, len(nb_h2o), len(nb_ion)))
				molecule_no +=1
				# here we hold water particles in contact
				nb_h2o = []
				# here we hold ions in contact
				nb_ion = []

			# to track when to change residue
			residue = particles[i]["residue"]		
			# iterate over all sides to find neighbours
			for side in cell.sides:
				n = side.neighbour
				if particles[n]['residue'] == "SOL" and not n in nb_h2o: nb_h2o.append(n)
				if particles[n]['residue'] == "ION" and not n in nb_ion: nb_ion.append(n)
					
		# need to write last lipid (we changed frames, not lpipid inside the frame)
		out_file.write("{0} {1} {2} {3}\n".format(frame, molecule_no, len(nb_h2o), len(nb_ion)))
		# count new frame
		frame += 1
except StopIteration:
	# pdb file ends
	pass