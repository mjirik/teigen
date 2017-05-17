#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging

logger = logging.getLogger(__name__)

import numpy as nm
import yaml
import argparse
import sys
import numpy as np
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
                 sphere_radius_compensation_factor=1.0
                 ):
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

    def add_cylinder(self, p1m, p2m, rad, id):
        """
        Funkce na vykresleni jednoho segmentu do 3D dat
        """
        pass

    def finish(self):
        # import ipdb; ipdb.set_trace()
        self.polyData = gen_tree(self.tree_data_old, self.cylinder_resolution, self.sphere_resolution,
                                 polygon_radius_selection_method=self.polygon_radius_selection_method,
                                 cylinder_radius_compensation_factor=self.cylinder_radius_compensation_factor,
                                 sphere_radius_compensation_factor=self.sphere_radius_compensation_factor)
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

    rot1 = vtk.vtkTransform()
    fi = nm.arccos(direction[1])

    rot1.RotateWXYZ(-nm.rad2deg(fi), 0.0, 0.0, 1.0)
    u = nm.abs(nm.sin(fi))
    rot2 = vtk.vtkTransform()
    if u > 1.0e-6:

        # sometimes d[0]/u little bit is over 1
        d0_over_u = direction[0] / u
        if d0_over_u > 1:
            psi = 0
        elif d0_over_u < -1:
            psi = 2 * nm.pi
        else:
            psi = nm.arccos(direction[0] / u)

        logger.debug('d0 ' + str(direction[0]) + '  u ' + str(u) + ' psi ' + str(psi))
        if direction[2] < 0:
            psi = 2 * nm.pi - psi

        rot2.RotateWXYZ(-nm.rad2deg(psi), 0.0, 1.0, 0.0)

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

    return tr2.GetOutput()


def get_sphere(center, radius, resolution=10):
    # create source
    import vtk
    source = vtk.vtkSphereSource()
    source.SetPhiResolution(resolution)
    source.SetThetaResolution(resolution)
    source.SetCenter(center[0], center[1], center[2])
    source.SetRadius(radius)
    source.Update()
    return source.GetOutput()

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

    elif polygon_radius_selection_method == "circumscribed":
        # from .. import geometry3d as g3
        factor = circumscribed_polygon_radius(cylinder_resolution)
        cylinder_radius_compensation_factor = factor
        sphere_radius_compensation_factor = factor

    elif polygon_radius_selection_method == "compensation factors":
        pass

    elif polygon_radius_selection_method == "cylinder surface":
        # from .. import geometry3d as g3
        radius_compensation_factor =  regular_polygon_perimeter_equivalent_radius(cylinder_resolution)
        cylinder_radius_compensation_factor = radius_compensation_factor
        sphere_radius_compensation_factor = radius_compensation_factor

    elif polygon_radius_selection_method == "cylinder volume":
        # from .. import geometry3d as g3
        radius_compensation_factor =  regular_polygon_surface_equivalent_radius(cylinder_resolution)
        cylinder_radius_compensation_factor = radius_compensation_factor
        sphere_radius_compensation_factor = radius_compensation_factor

    elif polygon_radius_selection_method == "average":
        # from .. import geometry3d as g3
        factor = circumscribed_polygon_radius(cylinder_resolution)
        factor = (factor + 1.0) / 2.0
        cylinder_radius_compensation_factor = factor
        sphere_radius_compensation_factor = factor


    elif polygon_radius_selection_method == "cylinder volume + sphere compensation":
        # analytically compensated cylinder + sphere compensate by measurement
        radius_compensation_factor =  regular_polygon_surface_equivalent_radius(cylinder_resolution)

        x = [6, 7, 8, 10, 12, 16, 20, 25, 30, 40, 50, 1000]
        y = [0.99820681, 0.99990171, 1.00057384, 1.00090875,
             1.00086617, 1.00064401, 1.00046984, 1.00032942,
             1.00024186, 1.00014509, 1.00009627, 1.]
        spl1 = InterpolatedUnivariateSpline(x, y)
        radius_compensation_factor *= 1. / spl1(cylinder_resolution)
        cylinder_radius_compensation_factor = radius_compensation_factor
        sphere_radius_compensation_factor = radius_compensation_factor

    elif polygon_radius_selection_method == "cylinder surface + sphere compensation":
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

    return cylinder_radius_compensation_factor, sphere_radius_compensation_factor


def gen_tree(tree_data, cylinder_resolution=10, sphere_resolution=10,
             polygon_radius_selection_method="inscribed",
             cylinder_radius_compensation_factor=1.0,
             sphere_radius_compensation_factor=1.0):
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
    if vtk.VTK_MAJOR_VERSION > 5:
        pass
    factors = polygon_radius_compensation_factos(
        polygon_radius_selection_method,
        cylinder_radius_compensation_factor,
        sphere_radius_compensation_factor,
        cylinder_resolution,
        sphere_resolution
    )

    cylinder_radius_compensation_factor, sphere_radius_compensation_factor = factors


    # import ipdb; ipdb.set_trace()
    for br in tree_data:
        # import ipdb;
        # ipdb.set_trace()

        dbg_msg = "generating edge " + str(br["length"])
        logger.debug(dbg_msg)
        # print(dbg_msg)
        radius = br['radius']
        cylinder_radius = radius * cylinder_radius_compensation_factor
        sphere_radius = radius * sphere_radius_compensation_factor
        cylinder = get_cylinder(br['upperVertex'],
                                br['length'],
                                cylinder_radius,
                                br['direction'],
                                resolution=cylinder_resolution)
        sphere1 = get_sphere(br['upperVertex'], sphere_radius, resolution=sphere_resolution)
        uv = br['upperVertex']
        length = br["length"]
        direction = br["direction"]
        # length = nm.linalg.norm(direction)
        # print "obj ", uv, length
        if length > 0:
            direction /= nm.linalg.norm(direction)

            lv = uv + direction * length
            sphere2 = get_sphere(lv, radius, resolution=sphere_resolution)

        if vtk.VTK_MAJOR_VERSION <= 5:
            appendFilter.AddInputConnection(cylinder.GetProducerPort())
            appendFilter.AddInputConnection(sphere1.GetProducerPort())
            appendFilter.AddInputConnection(sphere2.GetProducerPort())
        else:
            cylinderTri = vtk.vtkTriangleFilter()
            sphere1Tri = vtk.vtkTriangleFilter()
            sphere2Tri = vtk.vtkTriangleFilter()
            boolean_operation1 = vtk.vtkBooleanOperationPolyDataFilter()
            boolean_operation2 = vtk.vtkBooleanOperationPolyDataFilter()
            boolean_operation3 = vtk.vtkBooleanOperationPolyDataFilter()
            boolean_operation1.SetOperationToUnion()
            boolean_operation2.SetOperationToUnion()
            boolean_operation3.SetOperationToUnion()

            sphere1Tri.SetInputData(sphere1)
            sphere1Tri.Update()
            if length > 0:
                cylinderTri.SetInputData(cylinder)
                cylinderTri.Update()
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
                boolean_operation2 = sphere1Tri

            # this is simple version
            # appendFilter.AddInputData(boolean_operation2.GetOutput())
            # print "object connected, starting addind to general space " + str(br["length"])
            if appended_data is None:
                appended_data = boolean_operation2.GetOutput()
            else:
                boolean_operation3.SetInputData(0, appended_data)
                boolean_operation3.SetInputData(1, boolean_operation2.GetOutput())
                boolean_operation3.Update()
                appended_data = boolean_operation3.GetOutput()

    # import ipdb; ipdb.set_trace()

    if vtk.VTK_MAJOR_VERSION > 5:
        del (cylinderTri)
        del (sphere1Tri)
        del (sphere2Tri)
        del (boolean_operation1)
        del (boolean_operation2)
        del (boolean_operation3)
        pass
    logger.debug("konec gen_tree()")
    # appendFilter.Update()
    # appended_data = appendFilter.GetOutput()
    return appended_data


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

def regular_polygon_surface_equivalent_radius(n, radius=1.0):
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
