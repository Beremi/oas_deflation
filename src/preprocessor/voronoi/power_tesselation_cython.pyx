# distutils: language=c++
cimport cython
import numpy as np
cimport numpy as np
import random
from libc.stdio cimport printf
from libc.stdlib cimport rand, RAND_MAX
from libcpp.vector cimport vector
#from cpprandom cimport mt19937_64, uniform_real_distribution
from libcpp cimport bool
from libc.math cimport sin, cos, abs, fmin, sqrt
import sys

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.cdivision(True)
def _generate_voronoi_parts_2d_cython(diagram, container, list connection_list, list vertices,
                                        list regions, list ridge_vertices, list point_region):
    cdef:
        double _total_volume = 0.0
        int i, j, vert_id, pointA_id, pointB_id
        int id_, coni, vertex_start_num
        double x, y, z, r
        int diagram_len = diagram.size()
        double relative_difference = 1e-10
        vector[int] point_region_tmp
        vector[double] vertices_tmp
        vector[int] regions_tmp
        vector[int] point_con_tmp
        vector[int] connection_list_tmp
        vector[int] ridge_vertices_tmp
        bool in_list
        # list side_indices
        # list side_coords
    point_con_tmp.push_back(0)
    point_con_tmp.push_back(0)
    vertex_start_num = 0
    #for i, cell in enumerate(diagram):
    for i in range(diagram_len):
        cell = diagram.get_cell(i)
        point_region_tmp.push_back(i+1)
        (id_, x, y, z, r) = container.get(i)
        _total_volume += cell.volume()
        for side in cell.sides:
            # TODO: why not - if not np.allclose(np.array(side.as_coords())[:, -1], 0):
            side_coords = side.as_coords()
            side_indices = side.as_list()
            if not ((side_coords[0][2] == 0) and (side_coords[1][2] == 0)):
                pass
            else:
                for j in range(2):
                    vert_id = side_indices[j] + vertex_start_num
                    vertices_tmp.push_back(vert_id) 
                    vertices_tmp.push_back(side_coords[j][0])
                    vertices_tmp.push_back(side_coords[j][1])
                    vertices_tmp.push_back(side_coords[j][2])
                    regions_tmp.push_back(vert_id)

            if ((abs(side_coords[0][2]) < relative_difference) and (abs(side_coords[1][2]) < relative_difference)
                or (1-abs(side_coords[0][2]) < relative_difference) and (1-abs(side_coords[1][2]) < relative_difference)):
                continue

            pointA_id = i
            pointB_id = side.neighbour
            if pointA_id < pointB_id:
                point_con_tmp[0] = pointA_id
                point_con_tmp[1] = pointB_id
            else:
                point_con_tmp[0] = pointB_id
                point_con_tmp[1] = pointA_id
            if connection_list_tmp.size() == 0:
                for j in range(2):
                    connection_list_tmp.push_back(point_con_tmp[j])
                    ridge_vertices_tmp.push_back(side_indices[j] + vertex_start_num)
            else:
                in_list = False
                for coni in range(connection_list_tmp.size()//2):
                    if not ((connection_list_tmp[coni * 2] == point_con_tmp[0]) and (connection_list_tmp[coni*2+1] == point_con_tmp[1])):
                        pass
                    else:
                        in_list = True
                        continue
                if not in_list:
                    for j in range(2):
                        connection_list_tmp.push_back(point_con_tmp[j])
                        ridge_vertices_tmp.push_back(side_indices[j] + vertex_start_num)

        vertex_start_num += cell.vertices.size()

    for i in range(point_region_tmp.size()):
        point_region.append(point_region_tmp[i])

    for i in range(regions_tmp.size()):
        regions.append(regions_tmp[i])

    for i in range(ridge_vertices_tmp.size()):
        ridge_vertices.append(tuple([ridge_vertices_tmp[i*2], ridge_vertices_tmp[i*2+1]]))

    for i in range(connection_list_tmp.size()//2):
        connection_list.append(tuple([connection_list_tmp[i*2], connection_list_tmp[i*2+1]]))
    
    for i in range(vertices_tmp.size()//4):
        vertices.append(tuple([vertices_tmp[i*4], vertices_tmp[i*4+1], vertices_tmp[i*4+2], vertices_tmp[i*4+3]]))
    
    
