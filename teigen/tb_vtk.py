#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging

logger = logging.getLogger(__name__)

import numpy as nm
import yaml
import argparse
import sys
import numpy as np
import vtk
from scipy.interpolate import InterpolatedUnivariateSpline

# new interface

class TBVTK:
    """
    This generator is called by generateTree() function as a general form.
    Other similar generator is used for generating LAR outputs.
    """

    def __init__(self, gtree, cylinder_resolution=30, sphere_resolution=30,
                 polygon_radius_selection_method="inscribed",
                 cylinder_radius_compensation_factor=1.0,
                 sphere_radius_compensation_factor=1.0,
                 tube_shape=True
                 ):
        """

        :param gtree:
        :param cylinder_resolution:
        :param sphere_resolution:
        :param polygon_radius_selection_method:
        :param cylinder_radius_compensation_factor:
        :param sphere_radius_compensation_factor:
        :param tube_shape: If true, the tube shape is generated.
                Otherwise the cylinder shape is used.
        """
        # self.shape = gtree.shape
        # self.data3d = np.zeros(gtree.shape, dtype=np.int)
        # self.voxelsize_mm = gtree.voxelsize_mm
        # make comapatible with old system
        self.polygon_radius_selection_method = polygon_radius_selection_method
        self.cylinder_radius_compensation_factor = cylinder_radius_compensation_factor
        self.sphere_radius_compensation_factor = sphere_radius_compensation_factor
        self.tree_data = gtree.tree_data

        self.tree_data_old = compatibility_processing(self.tree_data)
        self.cylinder_resolution = cylinder_resolution
        self.sphere_resolution = sphere_resolution
        self.tube_shape = tube_shape

    def add_cylinder(self, p1m, p2m, rad, id):
        """
        Funkce na vykresleni jednoho segmentu do 3D dat
        """
        pass

    def finish(self):
        # import ipdb; ipdb.set_trace()
        self.polyData = gen_tree(
            self.tree_data_old, self.cylinder_resolution, self.sphere_resolution,
            polygon_radius_selection_method=self.polygon_radius_selection_method,
            cylinder_radius_compensation_factor=self.cylinder_radius_compensation_factor,
            sphere_radius_compensation_factor=self.sphere_radius_compensation_factor,
            tube_shape=self.tube_shape
        )
        # import ipdb; ipdb.set_trace()

    def get_output(self):
        return self.polyData

    def save(self, outputfile, lc_all="C"):

        import vtk
        logger.debug("vtk version " + str(vtk.VTK_BUILD_VERSION))
        # import ipdb; ipdb.set_trace()
        if lc_all is not None:
            import locale
            locale.setlocale(locale.LC_ALL, lc_all)
        writer = vtk.vtkPolyDataWriter()
        writer.SetFileName(outputfile)
        try:
            writer.SetInputData(self.polyData)
        except:
            logger.warning("old vtk is used")
            writer.SetInput(self.polyData)
        writer.Write()

    def show(self):
        logger.info("there is no show implemented")


# old interface

def move_to_position(src, upper, direction, axis0=2, axis1=1, axis2=0):

    # along axis 1 axis1=2, axis2=1

    r1 = [0.0, 0.0, 0.0]
    r2 = [0.0, 0.0, 0.0]

    r1[axis0] = 1.0
    r2[axis1] = 1.0

    print r1, r2

    rot1 = vtk.vtkTransform()
    fi = nm.arccos(direction[axis1])

    rot1.RotateWXYZ(-nm.rad2deg(fi), r1[0], r1[1], r1[2])
    u = nm.abs(nm.sin(fi))
    rot2 = vtk.vtkTransform()
    if u > 1.0e-6:

        # sometimes d[0]/u little bit is over 1
        d0_over_u = direction[axis2] / u
        if d0_over_u > 1:
            psi = 0
        elif d0_over_u < -1:
            psi = 2 * nm.pi
        else:
            psi = nm.arccos(direction[axis2] / u)

        # logger.debug('d0 ' + str(direction[axis2]) + '  u ' + str(u) + ' psi ' + str(psi))
        if direction[axis0] < 0:
            psi = 2 * nm.pi - psi

        rot2.RotateWXYZ(-nm.rad2deg(psi), r2[0], r2[1], r2[2])

    tl = vtk.vtkTransform()
    tl.Translate(upper)

    tr1a = vtk.vtkTransformFilter()
    if "SetInputConnection" in dir(tr1a):
        tr1a.SetInputConnection(src.GetOutputPort())
    else:
        tr1a.SetInput(src.GetOutput())
    tr1a.SetTransform(rot1)

    tr1b = vtk.vtkTransformFilter()
    if "SetInputConnection" in dir(tr1b):
        tr1b.SetInputConnection(tr1a.GetOutputPort())
    else:
        tr1b.SetInput(tr1a.GetOutput())
    # tr1b.SetInput(tr1a.GetOutput())
    tr1b.SetTransform(rot2)

    tr2 = vtk.vtkTransformFilter()
    if "SetInputConnection" in dir(tr2):
        tr2.SetInputConnection(tr1b.GetOutputPort())
    else:
        tr2.SetInput(tr1b.GetOutput())
    # tr2.SetInput(tr1b.GetOutput())
    tr2.SetTransform(tl)

    tr2.Update()

    return tr2
    # return tr2.GetOutput()


def get_tube(radius=1.0, point=[0.0, 0.0, 0.0],
             direction=[0.0, 0.0, 1.0], length=1.0,
             sphere_resolution=10, cylinder_resolution=10,
             cylinder_radius_compensation_factor=1.0,
             sphere_radius_compensation_factor=1.0,
             tube_shape=True, axis=1
             ):
    point1 = [0.0, 0.0, 0.0]
    center = [0.0, 0.0, 0.0]
    point2 = [0.0, 0.0, 0.0]

    center[axis] = length / 2.0
    point2[axis] = length

    cylinder_radius = radius * cylinder_radius_compensation_factor
    sphere_radius = radius * sphere_radius_compensation_factor

    # print "radius ",  cylinder_radius_compensation_factor, sphere_radius_compensation_factor

    direction /= nm.linalg.norm(direction)
    lv = point + direction * length

    cylinderTri = vtk.vtkTriangleFilter()
    sphere1Tri = vtk.vtkTriangleFilter()
    sphere2Tri = vtk.vtkTriangleFilter()

    cylinder = vtk.vtkCylinderSource()
    cylinder.SetCenter(center)
    cylinder.SetHeight(length)
    cylinder.SetRadius(cylinder_radius)
    cylinder.SetResolution(cylinder_resolution)
    cylinder.Update()
    cylinderTri.SetInputData(cylinder.GetOutput())
    cylinderTri.Update()

    sphere1 = get_sphere(
        center=point1,
        radius=sphere_radius,
        resolution=sphere_resolution,
        start_phi=0,
        #end_phi=90,
        end_phi=180,
        axis=1
    )
    # sphere1.Update()

    sphere1Tri.SetInputData(sphere1)
    sphere1Tri.Update()

    sphere2 = get_sphere(
        center=point2,
        # radius= 1. - (cylinder_radius - sphere_radius),
        radius=sphere_radius,
        resolution=sphere_resolution,
        start_phi=0,
        end_phi=180,
        axis=1

    )
    sphere2Tri.SetInputData(sphere2)
    sphere2Tri.Update()

    boolean_operation1 = vtk.vtkBooleanOperationPolyDataFilter()
    boolean_operation2 = vtk.vtkBooleanOperationPolyDataFilter()
    boolean_operation1.SetOperationToUnion()
    boolean_operation2.SetOperationToUnion()

    # booleanOperation.SetInputData(0, cyl)
    boolean_operation1.SetInputData(0, cylinderTri.GetOutput())
    boolean_operation1.SetInputData(1, sphere1Tri.GetOutput())
    boolean_operation1.Update()
    boolean_operation2.SetInputData(0, boolean_operation1.GetOutput())
    boolean_operation2.SetInputData(1, sphere2Tri.GetOutput())
    # booleanOperation.SetInputData(2, sph2)
    boolean_operation2.Update()
    # tube_in_base_position = boolean_operation2.GetOutput()

    tube = move_to_position(boolean_operation2, point, direction, 1, 2)
    return tube.GetOutput()


def get_cylinder(upper, height, radius,
                 direction,
                 resolution=10):
    import vtk
    src = vtk.vtkCylinderSource()
    src.SetCenter((0, height / 2, 0))
    # src.SetHeight(height + radius/2.0)
    src.SetHeight(height)
    src.SetRadius(radius)
    src.SetResolution(resolution)
    return move_to_position(src, upper, direction).GetOutput()


def get_sphere(center, radius, resolution=10, start_phi=None, end_phi=None, axis=0):
    sph = get_sphere_source(center, radius, resolution, start_phi, end_phi, axis)
    #sph.Update()
    return sph

def get_sphere_source(center, radius, resolution=10, start_phi=None, end_phi=None, axis=0):
    # create source
    import vtk
    sphere = vtk.vtkSphereSource()
    sphere.SetPhiResolution(resolution)
    sphere.SetThetaResolution(resolution)
    # sphere.SetCenter(center[0], center[1], center[2])
    sphere.SetCenter(.0, .0, .0)
    sphere.SetRadius(radius)
    if start_phi is not None:
        sphere.SetStartPhi(start_phi)
    if end_phi is not None:
        sphere.SetEndPhi(end_phi)

    if axis == 0:
        translate = vtk.vtkTransform()
        translate.Translate(center)

        tr2 = vtk.vtkTransformFilter()
        tr2.SetInputConnection(sphere.GetOutputPort())
        tr2.SetTransform(translate)
        sphere = tr2

    if axis == 1:
        #sphere = move_to_position(sphere, center, [1., 1., 0.])

        rot1 = vtk.vtkTransform()
        rot1.RotateWXYZ(90, 1, 0, 0)
        translate = vtk.vtkTransform()
        translate.Translate(center)

        tr1 = vtk.vtkTransformFilter()
        # tr1a.SetInputConnection(src.GetOutputPort())
        tr1.SetInputConnection(sphere.GetOutputPort())
        tr1.SetTransform(rot1)
        tr1.Update()

        tr2 = vtk.vtkTransformFilter()
        tr2.SetInputConnection(tr1.GetOutputPort())
        tr2.SetTransform(translate)
        sphere = tr2
    sphere.Update()
    return sphere.GetOutput()

def polygon_radius_compensation_factos(
        polygon_radius_selection_method,
        cylinder_radius_compensation_factor,
        sphere_radius_compensation_factor,
        cylinder_resolution,
        sphere_resolution
):

    if polygon_radius_selection_method == "inscribed":
        cylinder_radius_compensation_factor = 1.0
        sphere_radius_compensation_factor = 1.0
        cylinder_radius_compensation_factor_long = cylinder_radius_compensation_factor
        sphere_radius_compensation_factor_long = sphere_radius_compensation_factor

    elif polygon_radius_selection_method == "circumscribed":
        # from .. import geometry3d as g3
        factor = circumscribed_polygon_radius(cylinder_resolution)
        cylinder_radius_compensation_factor = factor
        sphere_radius_compensation_factor = factor
        cylinder_radius_compensation_factor_long = cylinder_radius_compensation_factor
        sphere_radius_compensation_factor_long = sphere_radius_compensation_factor

    elif polygon_radius_selection_method == "compensation factors":
        cylinder_radius_compensation_factor_long = cylinder_radius_compensation_factor
        sphere_radius_compensation_factor_long = sphere_radius_compensation_factor
        pass

    elif polygon_radius_selection_method == "cylinder surface":
        # from .. import geometry3d as g3
        radius_compensation_factor =  regular_polygon_perimeter_equivalent_radius(cylinder_resolution)
        cylinder_radius_compensation_factor = radius_compensation_factor
        sphere_radius_compensation_factor = radius_compensation_factor
        cylinder_radius_compensation_factor_long = cylinder_radius_compensation_factor
        sphere_radius_compensation_factor_long = sphere_radius_compensation_factor

    elif polygon_radius_selection_method == "cylinder volume":
        # from .. import geometry3d as g3
        radius_compensation_factor =  regular_polygon_area_equivalent_radius(cylinder_resolution)
        cylinder_radius_compensation_factor = radius_compensation_factor
        sphere_radius_compensation_factor = radius_compensation_factor
        cylinder_radius_compensation_factor_long = cylinder_radius_compensation_factor
        sphere_radius_compensation_factor_long = sphere_radius_compensation_factor

    elif polygon_radius_selection_method == "average":
        # from .. import geometry3d as g3
        factor = circumscribed_polygon_radius(cylinder_resolution)
        factor = (factor + 1.0) / 2.0
        cylinder_radius_compensation_factor = factor
        sphere_radius_compensation_factor = factor
        cylinder_radius_compensation_factor_long = cylinder_radius_compensation_factor
        sphere_radius_compensation_factor_long = sphere_radius_compensation_factor


    elif polygon_radius_selection_method == "cylinder volume + sphere error":
        # analytically compensated cylinder + sphere compensate by measurement

        # x = [6, 7, 8, 10, 12, 16, 18, 20, 22, 24, 26, 28, 30, 34, 38, 42, 46, 50, 100, 200]
        # y = [0.907761069989, 0.933196213162, 0.949394505649, 0.968085936432, 0.978051435286,
        #     0.987801062653, 0.990399429011, 0.99224803231, 0.993609871041, 0.994641988211,
        #     0.995442810978, 0.996076647961, 0.996586890825, 0.997348551083, 0.997881034775,
        #     0.998267843958, 0.998557648375, 0.998780367897, 1.0, 1.0]
        x = [6, 8, 10, 12, 14, 17, 21, 23, 29, 31, 39, 46, 50, 60, 70, 80, 90, 100, 150, 200]
        y = [0.907761087455, 0.949394472294, 0.968085949491, 0.978051451209, 0.983984947219, 0.989216674476, 0.992978216638, 0.994159974337, 0.996345087831, 0.996805451193, 0.997989063564, 0.998557649514, 0.998780374797, 0.999154602003, 0.999379702886, 0.999525551526, 0.999625411186, 0.999696768773, 0.999865478984, 1.0]
        # x = [6, 7, 8, 10, 12, 16, 20, 25, 30, 40, 50, 1000]
        # y = [0.99820681, 0.99990171, 1.00057384, 1.00090875,
        #      1.00086617, 1.00064401, 1.00046984, 1.00032942,
        #      1.00024186, 1.00014509, 1.00009627, 1.]
        spl1 = InterpolatedUnivariateSpline(x, y)
        # radius_compensation_factor *= 1. / spl1(cylinder_resolution)
        # cylinder_radius_compensation_factor = radius_compensation_factor
        sphere_radius_compensation_factor = 1. / spl1(cylinder_resolution)
        cylinder_radius_compensation_factor =  regular_polygon_area_equivalent_radius(cylinder_resolution)
        cylinder_radius_compensation_factor_long = cylinder_radius_compensation_factor
        sphere_radius_compensation_factor_long = sphere_radius_compensation_factor
        # cylinder_radius_compensation_factor = 1.0
        # sphere_radius_compensation_factor = 1.0

    elif polygon_radius_selection_method == "cylinder surface + sphere error polygon perimeter equivalent":
        # analytically compensated cylinder + sphere compensate by measurement
        radius_compensation_factor =  regular_polygon_perimeter_equivalent_radius(cylinder_resolution)
        x = [6, 8, 10, 12, 14, 17, 21, 23, 29, 31, 39, 46, 50, 60, 70, 80, 90, 100, 150, 200]
        y = [0.975228602567, 0.987423714247, 0.992397574304, 0.994910516827, 0.996355306056, 0.997593221093, 0.998459193066, 0.998726453568, 0.999213600518, 0.999314918173, 0.999572948543, 0.999695445089, 0.999743138504, 0.99982282109, 0.999870450728, 0.999901168474, 0.999922129169, 0.999937063217, 0.999972213153, 1.0]

        #x = [6, 7, 8, 10, 12, 16, 20, 25, 30, 40, 50, 100, 200]
        #y = [0.97522857799, 0.982858482408, 0.987423696432, 0.99239757445,
        #     0.994910515802, 0.997261880581, 0.998292863128, 0.998929760493,
        #     0.999266886531, 0.999594545757, 0.999743129653, 1.0, 1.0]
        spl1 = InterpolatedUnivariateSpline(x, y)
        radius_compensation_factor *= 1. / spl1(cylinder_resolution)
        cylinder_radius_compensation_factor = radius_compensation_factor
        sphere_radius_compensation_factor = radius_compensation_factor
        cylinder_radius_compensation_factor_long = cylinder_radius_compensation_factor
        sphere_radius_compensation_factor_long = sphere_radius_compensation_factor

    elif polygon_radius_selection_method == "cylinder surface + sphere error":
        # analytically compensated cylinder + sphere compensate by measurement
        radius_compensation_factor =  regular_polygon_perimeter_equivalent_radius(cylinder_resolution)
        x = [6, 7, 8, 10, 12, 16, 20, 25, 30, 40, 50, 100, 200]
        y = [0.97522857799, 0.982858482408, 0.987423696432, 0.99239757445,
             0.994910515802, 0.997261880581, 0.998292863128, 0.998929760493,
             0.999266886531, 0.999594545757, 0.999743129653, 1.0, 1.0]
        spl1 = InterpolatedUnivariateSpline(x, y)
        radius_compensation_factor *= 1. / spl1(cylinder_resolution)
        cylinder_radius_compensation_factor = radius_compensation_factor
        sphere_radius_compensation_factor = radius_compensation_factor
        cylinder_radius_compensation_factor_long = cylinder_radius_compensation_factor
        sphere_radius_compensation_factor_long = sphere_radius_compensation_factor

    elif polygon_radius_selection_method == "cylinder surface + sphere error + join error":
        # sphere like objects
        # analytically compensated cylinder + sphere compensate by measurement
        radius_compensation_factor =  regular_polygon_perimeter_equivalent_radius(cylinder_resolution)
        x = [6, 7, 8, 10, 12, 16, 20, 25, 30, 40, 50, 100, 200]
        y = [0.97522857799, 0.982858482408, 0.987423696432, 0.99239757445,
             0.994910515802, 0.997261880581, 0.998292863128, 0.998929760493,
             0.999266886531, 0.999594545757, 0.999743129653, 1.0, 1.0]
        spl1 = InterpolatedUnivariateSpline(x, y)
        radius_compensation_factor *= 1. / spl1(cylinder_resolution)
        cylinder_radius_compensation_factor = radius_compensation_factor
        sphere_radius_compensation_factor = radius_compensation_factor

        # long objects
        # analytically compensated cylinder + sphere compensate by measurement
        radius_compensation_factor =  regular_polygon_perimeter_equivalent_radius(cylinder_resolution)
        x = [6, 8, 10, 12, 16, 18, 20, 22, 26, 29, 33, 37, 42, 100, 200]
        y = [0.868491038036, 0.930734501628, 0.957118271516, 0.969229258189, 0.983939714778, 0.986159525166, 0.990680841143, 0.991986941042, 0.994813132295, 0.995235630724, 0.99633325877, 0.997090867513, 0.997631260904, 1.0, 1.0]
        spl1 = InterpolatedUnivariateSpline(x, y)
        radius_compensation_factor *= 1. / spl1(cylinder_resolution)
        cylinder_radius_compensation_factor_long = radius_compensation_factor
        sphere_radius_compensation_factor_long = radius_compensation_factor
    else:
        logger.error("Unknown compensation method '{}'".format(polygon_radius_selection_method))

    return cylinder_radius_compensation_factor, sphere_radius_compensation_factor, cylinder_radius_compensation_factor_long, sphere_radius_compensation_factor_long


def gen_tree(tree_data, cylinder_resolution=10, sphere_resolution=10,
             polygon_radius_selection_method="inscribed",
             cylinder_radius_compensation_factor=1.0,
             sphere_radius_compensation_factor=1.0,
             tube_shape=True
             ):
    """
    
    :param polygon_radius_selection_method:
        "inscribed":
        "compensation factors"
        "cylinder volume"
        "cylinder surface"
    :param tree_data:
    :param cylinder_resolution: 
    :param sphere_resolution: 
    :param cylinder_radius_compensation_factor: is used to change radius of cylinder and spheres
    :return: 
    """
    import vtk
    # appendFilter = vtk.vtkAppendPolyData()
    appended_data = None
    if vtk.VTK_MAJOR_VERSION <= 5:
        logger.error("VTK 6 required")
    factors = polygon_radius_compensation_factos(
        polygon_radius_selection_method,
        cylinder_radius_compensation_factor,
        sphere_radius_compensation_factor,
        cylinder_resolution,
        sphere_resolution
    )

    cylinder_radius_compensation_factor, sphere_radius_compensation_factor,\
    cylinder_radius_compensation_factor_long, sphere_radius_compensation_factor_long = factors
    print "factors ", factors


    # import ipdb; ipdb.set_trace()
    for br in tree_data:
        # import ipdb;
        # ipdb.set_trace()
        something_to_add = True
        radius = br['radius']
        length = br["length"]
        direction = br["direction"]
        uv = br['upperVertex']

        dbg_msg = "generating edge with length: " + str(br["length"])
        logger.debug(dbg_msg)

        # tube = get_tube_old(radius, uv, direction, length,
        if length == 0:
            tube = get_sphere(uv, radius * sphere_radius_compensation_factor, sphere_resolution)
        else:
            tube = get_tube(radius, uv, direction, length,
                            sphere_resolution, cylinder_resolution,
                            cylinder_radius_compensation_factor=cylinder_radius_compensation_factor_long,
                            sphere_radius_compensation_factor=sphere_radius_compensation_factor_long,
                            tube_shape=tube_shape)
        # this is simple version
        # appendFilter.AddInputData(boolean_operation2.GetOutput())
        # print "object connected, starting addind to general space " + str(br["length"])
        if something_to_add:
            if appended_data is None:
                #appended_data = boolean_operation2.GetOutput()
                appended_data = tube
            else:
                boolean_operation3 = vtk.vtkBooleanOperationPolyDataFilter()
                boolean_operation3.SetOperationToUnion()
                boolean_operation3.SetInputData(0, appended_data)
                boolean_operation3.SetInputData(1, tube)
                boolean_operation3.Update()
                appended_data = boolean_operation3.GetOutput()

    # import ipdb; ipdb.set_trace()


    # del (cylinderTri)
    # del (sphere1Tri)
    # del (sphere2Tri)
    # del (boolean_operation1)
    # del (boolean_operation2)
    # del (boolean_operation3)
    logger.debug("konec gen_tree()")
    # appendFilter.Update()
    # appended_data = appendFilter.GetOutput()
    return appended_data

def get_tube_old(radius, point, direction, length,
                 sphere_resolution, cylinder_resolution,
                 cylinder_radius_compensation_factor=1.0,
                 sphere_radius_compensation_factor=1.0,
                 tube_shape=True
                 ):
    """ Create a tube with ending half-spherese in Z axis.

    :param radius:
    :param point:
    :param direction:
    :param length:
    :param sphere_resolution:
    :param cylinder_resolution:
    :param cylinder_radius_compensation_factor:
    :param sphere_radius_compensation_factor:
    :param tube_shape:
    :return:
    """
    cylinderTri = vtk.vtkTriangleFilter()
    sphere1Tri = vtk.vtkTriangleFilter()
    sphere2Tri = vtk.vtkTriangleFilter()
    boolean_operation1 = vtk.vtkBooleanOperationPolyDataFilter()
    boolean_operation2 = vtk.vtkBooleanOperationPolyDataFilter()
    boolean_operation1.SetOperationToUnion()
    boolean_operation2.SetOperationToUnion()

    # print(dbg_msg)
    cylinder_radius = radius * cylinder_radius_compensation_factor
    sphere_radius = radius * sphere_radius_compensation_factor
    retval = None

    if tube_shape:
        sphere1 = get_sphere(point, sphere_radius, resolution=sphere_resolution)
        sphere1Tri.SetInputData(sphere1)
        sphere1Tri.Update()
    if length > 0:
        cylinder = get_cylinder(point,
                                length,
                                cylinder_radius,
                                direction,
                                resolution=cylinder_resolution)

        cylinderTri.SetInputData(cylinder)
        cylinderTri.Update()
        direction /= nm.linalg.norm(direction)

        lv = point + direction * length
        if tube_shape:
            sphere2 = get_sphere(lv, sphere_radius, resolution=sphere_resolution)
            sphere2Tri.SetInputData(sphere2)
            sphere2Tri.Update()

            # booleanOperation.SetInputData(0, cyl)
            boolean_operation1.SetInputData(0, cylinderTri.GetOutput())
            boolean_operation1.SetInputData(1, sphere1Tri.GetOutput())
            boolean_operation1.Update()
            boolean_operation2.SetInputData(0, boolean_operation1.GetOutput())
            boolean_operation2.SetInputData(1, sphere2Tri.GetOutput())
            # booleanOperation.SetInputData(2, sph2)
            boolean_operation2.Update()
        else:
            boolean_operation2 = cylinderTri

        retval = boolean_operation2.GetOutput()
    else:
        if tube_shape:
            boolean_operation2 = sphere1Tri
            retval = boolean_operation2.GetOutput()
        else:
            # length == 0 but no spheres
            # so we are generating just flat shape
            # boolean_operation2 = cylinderTri

            # return empty space
            retval = vtk.vtkTriangleFilter().GetOutput()
            # something_to_add = False
    return retval


def gen_tree_old(tree_data):
    import vtk
    points = vtk.vtkPoints()
    polyData = vtk.vtkPolyData()
    polyData.Allocate(1000, 1)
    polyData.SetPoints(points)
    poffset = 0

    for br in tree_data:
        logger.debug("generating edge " + str(br["length"]))
        cyl = get_cylinder(br['upperVertex'],
                           br['length'],
                           br['radius'],
                           br['direction'],
                           resolution=16)

        for ii in xrange(cyl.GetNumberOfPoints()):
            points.InsertPoint(poffset + ii, cyl.GetPoint(ii))

        for ii in xrange(cyl.GetNumberOfCells()):
            cell = cyl.GetCell(ii)
            cellIds = cell.GetPointIds()
            for jj in xrange(cellIds.GetNumberOfIds()):
                oldId = cellIds.GetId(jj)
                cellIds.SetId(jj, oldId + poffset)

            polyData.InsertNextCell(cell.GetCellType(),
                                    cell.GetPointIds())
        # spheres part
        # cyl = get_sphere(br['upperVertex'],
        #                    br['radius']
        #                    )
        # for ii in xrange(cyl.GetNumberOfPoints()):
        #     points.InsertPoint(poffset + ii, cyl.GetPoint(ii))
        #
        # for ii in xrange(cyl.GetNumberOfCells()):
        #     cell = cyl.GetCell(ii)
        #     cellIds = cell.GetPointIds()
        #     for jj in xrange(cellIds.GetNumberOfIds()):
        #         oldId = cellIds.GetId(jj)
        #         cellIds.SetId(jj, oldId + poffset)
        #
        #     polyData.InsertNextCell(cell.GetCellType(),
        #                             cell.GetPointIds())

        poffset += cyl.GetNumberOfPoints()

    return polyData


def compatibility_processing(indata):
    scale = 1e-3
    scale = 1

    outdata = []
    for key in indata:
        ii = indata[key]
        # logger.debug(ii)
        br = {}

        lengthEstimation = None
        try:
            # old version of yaml tree
            vA = ii['upperVertexXYZmm']
            vB = ii['lowerVertexXYZmm']
            radi = ii['radius']
            lengthEstimation = ii['length']
        except:
            # new version of yaml
            try:
                vA = ii['nodeA_ZYX_mm']
                vB = ii['nodeB_ZYX_mm']
                radi = ii['radius_mm']
                if "lengthEstimation" in ii.keys():
                    lengthEstimation = ii['lengthEstimation']
            except:
                import traceback
                logger.debug(traceback.format_exc())
                continue

        br['upperVertex'] = nm.array(vA) * scale
        br['radius'] = radi * scale
        if lengthEstimation is None:

            br['real_length'] = None
        else:
            br['real_length'] = lengthEstimation * scale

        vv = nm.array(vB) * scale - br['upperVertex']
        br['direction'] = vv / nm.linalg.norm(vv)
        br['length'] = nm.linalg.norm(vv)
        outdata.append(br)

    return outdata


def fix_tree_structure(tree_raw_data):
    if 'graph' in tree_raw_data:
        trees = tree_raw_data['graph']
    else:
        trees = tree_raw_data['Graph']
    return trees


def vt_file_2_vtk_file(infile, outfile, text_label=None):
    """
    From vessel_tree.yaml to output.vtk

    :param vessel_tree:  vt structure
    :param outfile: filename with .vtk extension
    :param text_label: text label like 'porta' or 'hepatic_veins'
    :return:

    """
    yaml_file = open(infile, 'r')
    tree_raw_data = yaml.load(yaml_file)
    vt2vtk_file(tree_raw_data, outfile, text_label)


def vt2vtk_file(vessel_tree, outfile, text_label=None):
    """
    vessel_tree structure
    :param vessel_tree:  vt structure
    :param outfile: filename with .vtk extension
    :param text_label: text label like 'porta' or 'hepatic_veins'
    :return:
    """
    import vtk
    trees = fix_tree_structure(vessel_tree)

    tkeys = trees.keys()
    if text_label is None:
        text_label = tkeys[0]

    tree_data = compatibility_processing(trees[text_label])
    polyData = gen_tree(tree_data)

    import ipdb;
    ipdb.set_trace()
    writer = vtk.vtkPolyDataWriter()
    writer.SetFileName(outfile)
    try:
        writer.SetInputData(polyData)
    except:
        logger.warning("old vtk is used")
        writer.SetInput(polyData)
    writer.Write()


# From  geometry3.py

def circumscribed_polygon_radius(n, radius=1.0):
    """ Get circumscribed polygon radius.

    :param n: number of polygon elements
    :param radius: radius of inscribed circle
    :return: radius (distance from center to the corner) of polygon circumscribed to the
     circle
    """

    theta = 2 * np.pi / n
    radius_out = radius / np.cos(theta/2)

    return radius_out

def inscribed_polygon_radius(radius, n):
    """

    :param radius:
    :return:
    """
    pass

def regular_polygon_area_equivalent_radius(n, radius=1.0):
    """ Compute equivalent radius to obtain same surface as circle.

    \theta = \frac{2 \pi}{n}

    r_{eqs} = \sqrt{\frac{\theta r^2}{\sin{\theta}}}

    :param radius: circle radius
    :param n:  number of regular polygon segments
    :return:  equivalent regular polynom surface
    """

    theta = 2 * np.pi / n

    r = np.sqrt((theta * radius**2) / np.sin(theta))
    return r

def regular_polygon_perimeter_equivalent_radius(n, radius=1.0):
    """ Compute equivalent radius to obtain same perimeter as circle.

    \theta = \frac{2 \pi}{n}

    r_{eqp} = \frac{\theta r}{2 \sin{\frac{\theta}}{2}}

    :param radius: circle radius
    :param n:  number of regular polygon segments
    :return:  equivalent regular polynom surface
    """

    theta = 2 * np.pi / n

    r = (theta * radius) / (2 * np.sin(theta/2.0))
    return r

def main():
    logger = logging.getLogger()

    logger.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    logger.addHandler(ch)

    # create file handler which logs even debug messages
    # fh = logging.FileHandler('log.txt')
    # fh.setLevel(logging.DEBUG)
    # formatter = logging.Formatter(
    #     '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    # fh.setFormatter(formatter)
    # logger.addHandler(fh)
    # logger.debug('start')

    # input parser
    parser = argparse.ArgumentParser(
        description=__doc__
    )
    parser.add_argument(
        'inputfile',
        default=None,
        help='input file'
    )
    parser.add_argument(
        'outputfile',
        default='output.vtk',
        nargs='?',
        help='output file'
    )
    parser.add_argument(
        '-l', '--label',
        default=None,
        help='text label of vessel tree. f.e. "porta" or "hepatic_veins". \
        First label is used if it is set to None'
    )
    parser.add_argument(
        '-d', '--debug', action='store_true',
        help='Debug mode')
    args = parser.parse_args()

    if args.debug:
        ch.setLevel(logging.DEBUG)

    vt_file_2_vtk_file(args.inputfile, args.outputfile, args.label)


if __name__ == "__main__":
    main()
