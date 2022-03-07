
from random import random
from random import seed

import numpy as np
from scipy.spatial import ConvexHull
import vtk
from vtk.util import numpy_support


def voronoi_plot_3d_vtk(vor):
    if vor.points.shape[1] != 3:
        raise ValueError("Voronoi diagram is not 3-D")

    colors = vtk.vtkNamedColors()

    renderer = vtk.vtkRenderer()
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.SetWindowName('Polyhedron')
    renderWindow.AddRenderer(renderer)
    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)

    # create polyhedron (cube)
    # The point Ids are: [0, 1, 2, 3, 4, 5, 6, 7]

    points = vtk.vtkPoints()

    for point in vor.vertices:
        points.InsertNextPoint(point)

    # points.SetData(vtk_np.numpy_to_vtk(points_arr))

    # These are the point ids corresponding to each face.
    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(points)
    #faces = [[0, 3, 2, 1], [0, 4, 7, 3], [4, 5, 6, 7], [5, 1, 2, 6], [0, 1, 5, 4], [2, 3, 7, 6]]
#    for i in range(len(vor.ridge_vertices)):
#        faceId = vtk.vtkIdList()
#        faceId.InsertNextId(1)  # Six faces make up the cell.

#        face = vor.ridge_vertices[i]
#        print(face)
#        faceId.InsertNextId(len(face))  # The number of points in the face.
#        [faceId.InsertNextId(i) for i in face]

#        ugrid.InsertNextCell(vtk.VTK_POLYHEDRON, faceId)
#    print(dir(ugrid))

    for pidx in range(vor.points.shape[0]):
        pridx = vor.point_region[pidx]
        if not pridx:
            print('XXXX', pidx + 1, pridx)
            continue
        if (-1 in vor.regions[pridx]):
            continue
        faceId = vtk.vtkIdList()
        faces = []
        for rpidx, rp in enumerate(vor.ridge_points):
            if pidx in rp:
                faces.append(vor.ridge_vertices[rpidx])
        faces = [i for i in faces if not -1 in i]
        faceId.InsertNextId(len(faces))
        for face in faces:
            faceId.InsertNextId(len(face))
            [faceId.InsertNextId(i) for i in face]
        ugrid.InsertNextCell(vtk.VTK_POLYHEDRON, faceId)

    scalar3 = np.arange(vor.points.shape[0])
    scalar3_array = numpy_support.numpy_to_vtk(scalar3)
    scalar3_array.SetName('scalar3')
    cell_data = ugrid.GetCellData()
    cell_data.AddArray(scalar3_array)

    # Here we write out the cube.
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetInputData(ugrid)
    writer.SetFileName('polyhedron.vtu')
    #writer.SetDataModeToAscii()
    writer.SetDataModeToBinary()
    #writer.Write()
    writer.Update()

    points2 = vtk.vtkPoints()
    for point in vor.points:
        points2.InsertNextPoint(point)
    ugrid2 = vtk.vtkUnstructuredGrid()
    ugrid2.SetPoints(points2)
    if hasattr(vor, 'weights'):
        scalar3 = vor.weights
    else:
        scalar3 = np.arange(vor.points.shape[0])
    scalar3_array = numpy_support.numpy_to_vtk(scalar3)
    scalar3_array.SetName('radii')
    point_data2 = ugrid2.GetPointData()
    point_data2.AddArray(scalar3_array)
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetInputData(ugrid2)
    writer.SetFileName('points.vtu')
    #writer.SetDataModeToAscii()
    writer.SetDataModeToBinary()
    #writer.Write()
    writer.Update()

    # Create a mapper and actor
    mapper = vtk.vtkDataSetMapper()
    mapper.SetInputData(ugrid)

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(colors.GetColor3d('Silver'))
    #actor.GetProperty().SetOpacity(.5)
    #actor.GetProperty().SetColor([1.0, 0.0, 0.0])
    #actor.GetProperty().SetLineWidth(3)
    actor.GetProperty().EdgeVisibilityOn()


    # Visualize
    renderer.AddActor(actor)
    renderer.SetBackground(colors.GetColor3d('Salmon'))
    renderer.ResetCamera()
    renderer.GetActiveCamera().Azimuth(30)
    renderer.GetActiveCamera().Elevation(30)
    renderWindow.Render()

    exporter = vtk.vtkVRMLExporter()
    exporter.SetRenderWindow(renderWindow)
    exporter.SetFileName("polyhedron.wrl")
    exporter.Write()
    exporter.Update()

    renderWindowInteractor.Start()


def voronoi_plot_3d_mlab(vor, fig=None, **kw):
    """
    Plot the given Voronoi diagram in 3-D
    Parameters
    ----------
    vor : scipy.spatial.Voronoi or power.PowerTesselation instance
        Diagram to plot
    ax : matplotlib.axes.Axes instance, optional
        Axes to plot on
    show_points: bool, optional
        Add the Voronoi points to the plot.
    show_vertices : bool, optional
        Add the Voronoi vertices to the plot.
    line_colors : string, optional
        Specifies the line color for polygon boundaries
    line_width : float, optional
        Specifies the line width for polygon boundaries
    line_alpha: float, optional
        Specifies the line alpha for polygon boundaries
    point_size: float, optional
        Specifies the size of points
    Returns
    -------
    fig : matplotlib.figure.Figure instance
        Figure for the plot
    See Also
    --------
    Voronoi
    Notes
    -----
    Requires Mayavi.
    Examples
    --------
    Set of point:
    >>> from mayavi import mlab
    >>> points = np.random.rand(10, 3) #random
    Voronoi diagram of the points:
    >>> from power import PowerTesselation
    >>> vor = PowerTesselation(points)
    using `voronoi_plot_3d` for visualisation:
    >>> fig = voronoi_plot_3d(vor)
    using `voronoi_plot_3d` for visualisation with enhancements:
    >>> fig = voronoi_plot_3d(vor, show_vertices=False, line_colors='orange',
    ...                 line_width=2, line_alpha=0.6, point_size=2)
    >>> mlab.show()
    """
    from mayavi import mlab
    if vor.points.shape[1] != 3:
        raise ValueError("Voronoi diagram is not 3-D")

    if kw.get('show_points', True):
        point_size = kw.get('point_size', None)
        #ax.plot(vor.points[:,0], vor.points[:,1], '.', markersize=point_size)
    if kw.get('show_vertices', True): pass
        #ax.plot(vor.vertices[:,0], vor.vertices[:,1], 'o')

    line_colors = kw.get('line_colors', 'k')
    line_width = kw.get('line_width', 1.0)
    line_alpha = kw.get('line_alpha', 1.0)

    if fig is None:
        fig = mlab.figure(1, size=(1000, 500), fgcolor=(0, 0, 0),
                          bgcolor=(1., 1., 1.))

    center = vor.points.mean(axis=0)
    ptp_bound = vor.points.ptp(axis=0)

    finite_segments = []
    infinite_segments = []
    for pointidx, simplex in zip(vor.ridge_points, vor.ridge_vertices):
        simplex = np.asarray(simplex)
        if np.all(simplex >= 0):
            finite_segments.append(vor.vertices[simplex])
        else:
            i = simplex[simplex >= 0][0]  # finite end Voronoi vertex

            t = vor.points[pointidx[1]] - vor.points[pointidx[0]]  # tangent
            t /= np.linalg.norm(t)
            n = np.array([-t[1], t[0]])  # normal

            midpoint = vor.points[pointidx].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, n)) * n
            far_point = vor.vertices[i] + direction * ptp_bound.max()

            infinite_segments.append([vor.vertices[i], far_point])

    fig.scene.parallel_projection = True

    fig.scene.disable_render = True
    mlab_points = []
    mlab_cells = []
    print(vor.points.shape, vor.weights.shape, len(vor.regions))
    for (x, y, z), r, region in zip(vor.points, vor.weights, vor.regions[1:]):
        # if prev_obj:
        #    prev_obj.actor.property.opacity = 1.0
        #    prev_obj.actor.property.color = (1,1,1)
        pts = mlab.quiver3d(x, y, z, r, r, r, scalars=r, mode='sphere',
                            opacity=1., scale_mode='none', resolution=30,
                            transparent=False)
        pts.glyph.color_mode = 'color_by_vector'
        pts.glyph.glyph.scaling = False
        pts.glyph.glyph_source.glyph_source.radius = r
        pts.glyph.glyph_source.glyph_source.center = [0, 0, 0]
        color = (255, 0, 0, 150)
        rgb_lut = np.array([color, color])
        scalar_lut_manager = pts.module_manager.scalar_lut_manager
        scalar_lut_manager.lut._vtk_obj.SetTableRange(0, rgb_lut.shape[0])
        scalar_lut_manager.lut.number_of_colors = rgb_lut.shape[0]
        scalar_lut_manager.lut.table = rgb_lut
        pts.module_manager.visible = True
        pts.actor.property.frontface_culling = True

        vertices = np.array(vor.vertices)[region]
        convexhull = ConvexHull(vertices)
        x = vertices[:, 0]
        y = vertices[:, 1]
        z = vertices[:, 2]
        obj1 = mlab.triangular_mesh(x, y, z, convexhull.simplices,
                                    color=(0, 1, 0), opacity=1.0,
                                    transparent=False)
        obj2 = mlab.triangular_mesh(x, y, z, convexhull.simplices,
                                    color=(0, 0, 0), line_width=1.,
                                    representation='wireframe')

        obj1.actor.property.lighting = False
        obj2.actor.property.lighting = False
        obj1.module_manager.visible = True
        obj2.module_manager.visible = True
        obj1.actor.property.frontface_culling = True
        obj2.actor.property.frontface_culling = True
        obj1.actor.property.backface_culling = True
        obj2.actor.property.backface_culling = True
        mlab_points.append(pts)
        mlab_cells.append((obj1, obj2))

    mlab.outline(extent=[0, 10, 0, 10, 0, 10], color=(0, 0, 0))
    fig.scene.isometric_view()
    fig.scene.disable_render = False

    return fig


def voronoi_plot_3d_mlab2(vor, fig=None, **kw):
    """
    Plot the given Voronoi diagram in 3-D
    Parameters
    ----------
    vor : scipy.spatial.Voronoi or power.PowerTesselation instance
        Diagram to plot
    ax : matplotlib.axes.Axes instance, optional
        Axes to plot on
    show_points: bool, optional
        Add the Voronoi points to the plot.
    show_vertices : bool, optional
        Add the Voronoi vertices to the plot.
    line_colors : string, optional
        Specifies the line color for polygon boundaries
    line_width : float, optional
        Specifies the line width for polygon boundaries
    line_alpha: float, optional
        Specifies the line alpha for polygon boundaries
    point_size: float, optional
        Specifies the size of points
    Returns
    -------
    fig : matplotlib.figure.Figure instance
        Figure for the plot
    See Also
    --------
    Voronoi
    Notes
    -----
    Requires Mayavi.
    Examples
    --------
    Set of point:
    >>> from mayavi import mlab
    >>> points = np.random.rand(10, 3) #random
    Voronoi diagram of the points:
    >>> from power import PowerTesselation
    >>> vor = PowerTesselation(points)
    using `voronoi_plot_3d` for visualisation:
    >>> fig = voronoi_plot_3d(vor)
    using `voronoi_plot_3d` for visualisation with enhancements:
    >>> fig = voronoi_plot_3d(vor, show_vertices=False, line_colors='orange',
    ...                 line_width=2, line_alpha=0.6, point_size=2)
    >>> mlab.show()
    """
    from mayavi import mlab
    if vor.points.shape[1] != 3:
        raise ValueError("Voronoi diagram is not 3-D")

    if fig is None:
        fig = mlab.figure(1, size=(1000, 500), fgcolor=(0, 0, 0),
                          bgcolor=(1., 1., 1.))

    center = vor.points.mean(axis=0)
    ptp_bound = vor.points.ptp(axis=0)

    finite_segments = []
    infinite_segments = []
    for pointidx, simplex in zip(vor.ridge_points, vor.ridge_vertices):
        simplex = np.asarray(simplex)
        if np.all(simplex >= 0):
            finite_segments.append(vor.vertices[simplex])
        else:
            i = simplex[simplex >= 0][0]  # finite end Voronoi vertex

            t = vor.points[pointidx[1]] - vor.points[pointidx[0]]  # tangent
            t /= np.linalg.norm(t)
            n = np.array([-t[1], t[0]])  # normal

            midpoint = vor.points[pointidx].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, n)) * n
            far_point = vor.vertices[i] + direction * ptp_bound.max()

            infinite_segments.append([vor.vertices[i], far_point])

    fig.scene.parallel_projection = True

    fig.scene.disable_render = True
    mlab_points = []
    mlab_cells = []
    print(vor.points.shape, vor.weights.shape, len(vor.regions))
    for (x, y, z), r, region in zip(vor.points, vor.weights, vor.regions[1:]):
        # if prev_obj:
        #    prev_obj.actor.property.opacity = 1.0
        #    prev_obj.actor.property.color = (1,1,1)
        pts = mlab.quiver3d(x, y, z, r, r, r, scalars=r, mode='sphere',
                            opacity=1., scale_mode='none', resolution=30,
                            transparent=False)
        pts.glyph.color_mode = 'color_by_vector'
        pts.glyph.glyph.scaling = False
        pts.glyph.glyph_source.glyph_source.radius = r
        pts.glyph.glyph_source.glyph_source.center = [0, 0, 0]
        color = (255, 0, 0, 150)
        rgb_lut = np.array([color, color])
        scalar_lut_manager = pts.module_manager.scalar_lut_manager
        scalar_lut_manager.lut._vtk_obj.SetTableRange(0, rgb_lut.shape[0])
        scalar_lut_manager.lut.number_of_colors = rgb_lut.shape[0]
        scalar_lut_manager.lut.table = rgb_lut
        pts.module_manager.visible = True
        pts.actor.property.frontface_culling = True

        vertices = np.array(vor.vertices)[region]
        convexhull = ConvexHull(vertices)
        x = vertices[:, 0]
        y = vertices[:, 1]
        z = vertices[:, 2]
        obj1 = mlab.triangular_mesh(x, y, z, convexhull.simplices,
                                    color=(0, 1, 0), opacity=1.0,
                                    transparent=False)
        obj2 = mlab.triangular_mesh(x, y, z, convexhull.simplices,
                                    color=(0, 0, 0), line_width=1.,
                                    representation='wireframe')

        obj1.actor.property.lighting = False
        obj2.actor.property.lighting = False
        obj1.module_manager.visible = True
        obj2.module_manager.visible = True
        obj1.actor.property.frontface_culling = True
        obj2.actor.property.frontface_culling = True
        obj1.actor.property.backface_culling = True
        obj2.actor.property.backface_culling = True
        mlab_points.append(pts)
        mlab_cells.append((obj1, obj2))

    mlab.outline(extent=[0, 10, 0, 10, 0, 10], color=(0, 0, 0))
    fig.scene.isometric_view()
    fig.scene.disable_render = False

    return fig


def tesselation_plot_3d_mlab(tes, animate=False):
    from mayavi import mlab
    # mlab.options.offscreen = True
    visible = True

    fig = mlab.figure(1, size=(1000, 500), fgcolor=(0, 0, 0),
                      bgcolor=(1., 1., 1.))
    # fig.scene.renderer.use_depth_peeling = True

    # x = points[:, 0]
    # y = points[:, 1]
    # z = points[:, 2]
    # obj = mlab.points3d(x, y, z, radii, color=(1,0,0), opacity=1.)

    # mlab.view(-60.0, 70.0, focalpoint = [0., 1., 1.])
    fig.scene.parallel_projection = True

    fig.scene.disable_render = True
    mlab_points = []
    mlab_cells = []
    for (x, y, z), r, region in zip(tes.points, tes.radii, tes.regions):
        # if prev_obj:
        #    prev_obj.actor.property.opacity = 1.0
        #    prev_obj.actor.property.color = (1,1,1)
        pts = mlab.quiver3d(x, y, z, r, r, r, scalars=r, mode='sphere',
                            opacity=1., scale_mode='none', resolution=30,
                            transparent=False)
        pts.glyph.color_mode = 'color_by_vector'
        pts.glyph.glyph.scaling = False
        pts.glyph.glyph_source.glyph_source.radius = r
        pts.glyph.glyph_source.glyph_source.center = [0, 0, 0]
        color = (255, 0, 0, 150)
        rgb_lut = np.array([color, color])
        scalar_lut_manager = pts.module_manager.scalar_lut_manager
        scalar_lut_manager.lut._vtk_obj.SetTableRange(0, rgb_lut.shape[0])
        scalar_lut_manager.lut.number_of_colors = rgb_lut.shape[0]
        scalar_lut_manager.lut.table = rgb_lut
        pts.module_manager.visible = visible
        pts.actor.property.frontface_culling = True

        vertices = np.array(tes.vertices)[:, 1:][region]
        convexhull = ConvexHull(vertices)
        x = vertices[:, 0]
        y = vertices[:, 1]
        z = vertices[:, 2]
        obj1 = mlab.triangular_mesh(x, y, z, convexhull.simplices,
                                    color=(0, 1, 0), opacity=1.0,
                                    transparent=False) # scalars=x, vmin=0, vmax=10)
        obj2 = mlab.triangular_mesh(x, y, z, convexhull.simplices,
                                    color=(0, 0, 0), line_width=1.,
                                    representation='wireframe')

        obj1.actor.property.lighting = False
        obj2.actor.property.lighting = False
        obj1.module_manager.visible = visible
        obj2.module_manager.visible = visible
        obj1.actor.property.frontface_culling = True
        obj2.actor.property.frontface_culling = True
        obj1.actor.property.backface_culling = True
        obj2.actor.property.backface_culling = True
        mlab_points.append(pts)
        mlab_cells.append((obj1, obj2))

    mlab.outline(extent=[0, 10, 0, 10, 0, 10], color=(0, 0, 0))
    fig.scene.isometric_view()
    fig.scene.disable_render = False
    # mlab.savefig('test.png')
    # mlab.show()
    if animate:
        @mlab.show
        @mlab.animate()
        def anim(mlab_points, mlab_cells):
            prev_obj1 = None
            prev_obj2 = None
            prev_pts = None
            # idx = 0
            for pts, (obj1, obj2) in zip(mlab_points, mlab_cells):
                if prev_obj1:
                    prev_obj1.actor.property.opacity = 1.0
                    prev_obj1.actor.property.lighting = False
                    # prev_obj1.actor.property.colormap = 'jet'
                    prev_obj1.module_manager.scalar_lut_manager.lut_mode = 'jet'
                    prev_obj1.actor.mapper.scalar_visibility = True
                    prev_obj1.actor.property.color = (.7, .7, .7)
                    prev_pts.module_manager.visible = False
                pts.module_manager.visible = True
                obj1.module_manager.visible = True
                obj2.module_manager.visible = True
                prev_obj1 = obj1
                prev_obj2 = obj2
                prev_pts = pts
                # mlab.savefig(filename='test%d.png' % idx)
                # idx += 1
                yield
        fig.scene.movie_maker.record = True
        anim(mlab_points, mlab_cells)
        # mlab.show()
    return fig
