from pydmga.geometry import OrthogonalGeometry
from pydmga.geometry import CastOrthoPlaneGeometry
from pydmga.container import Container
from pydmga.diagram import Diagram
from pydmga.io import pdb
filename = "drec-1_590-600ns"
groups = ["GN2"]
end_lipid_at = "UDC"
start_lipid_at = "KDO"
in_file = file("{0}.pdb".format(filename), "r")
out_file = file("{0}.dat".format(filename), "w")
try:
	frame = 0
	out_file.write("frame lipid sa\n")
	while True:
		lipid_by_atom = []
		residue = ''
		molecule_id = 0
		while True:
			item = pdb.PDBRecord(pdb.readnext(in_file))
			if item.type() == "CRYST1":
				(a, b, c) = (item["a"], item["b"], item["c"])
				box_geometry = OrthogonalGeometry(a, b, c, True, True, True)
				geometry_up = CastOrthoPlaneGeometry(box_geometry, (0, 0, 1))
				container_up = Container(geometry_up)
			if item.type() == "ATOM":
				if residue == end_lipid_at and item['residue'] == start_lipid_at:
					molecule_id += 1
				if item["residue"] in groups and item["atom"][0] != 'H' and item["atom"][1] != 'H':										
					if item['z'] > c / 2:
						container_up.add(item.as_particle())
						lipid_by_atom.append(molecule_id)
				residue = item["residue"]
			if item.type() == "ENDMDL" or item.type() == "TER":
				break
		diagram_up = Diagram(container_up)		
		all_area = 0.0
		molecule_id = lipid_by_atom[0]
		sa = 0.0
		for i, cell in enumerate(diagram_up):
			if molecule_id != lipid_by_atom[i]:
				out_file.write("{0} {1} {2}\n".format(frame, molecule_id, sa))
				sa = 0.0
				molecule_id = lipid_by_atom[i]
			for side in cell.sides:
				neighbour = side.neighbour
				if geometry_up.on_boundary(neighbour):
					sa += side.area()
		out_file.write("{0} {1} {2}\n".format(frame, molecule_id, sa))
		frame += 1
except StopIteration:
	pass