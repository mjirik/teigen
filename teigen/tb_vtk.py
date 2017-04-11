#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
logger = logging.getLogger(__name__)

import numpy as nm
import yaml
import argparse
import sys
import numpy as np


# new interface

class TBVTK:
    """
    This generator is called by generateTree() function as a general form.
    Other similar generator is used for generating LAR outputs.
    """
    def __init__(self, gtree, cylinder_resolution=50, sphere_resolution=50):
        # self.shape = gtree.shape
        # self.data3d = np.zeros(gtree.shape, dtype=np.int)
        # self.voxelsize_mm = gtree.voxelsize_mm
        # make comapatible with old system
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
        self.polyData = gen_tree(self.tree_data_old, self.cylinder_resolution, self.sphere_resolution)
        # import ipdb; ipdb.set_trace()

    def get_output(self):
        return self.polyData

    def save(self, outputfile):

        import vtk
        logger.debug("vtk version " + str(vtk.VTK_BUILD_VERSION))
        # import ipdb; ipdb.set_trace()
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
    src.SetCenter((0, height/2, 0))
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

        logger.debug('d0 '+str(direction[0])+'  u '+str(u)+' psi '+str(psi))
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

def gen_tree(tree_data, cylinder_resolution=10, sphere_resolution=10):
    import vtk
    # appendFilter = vtk.vtkAppendPolyData()
    appendedData = None

    for br in tree_data:

        logger.debug("generating edge " + str(br["length"]))
        cylinder = get_cylinder(br['upperVertex'],
                           br['length'],
                           br['radius'],
                           br['direction'],
                           resolution=cylinder_resolution)
        sphere1 = get_sphere(br['upperVertex'], br['radius'], resolution=sphere_resolution )
        uv = br['upperVertex']
        # length = br["length"]
        direction = br["direction"]
        length = nm.linalg.norm(direction)
        # print "obj ", uv, length
        if length > 0:

            direction = direction / length

            lv = uv + direction * length
            sphere2 = get_sphere(lv, br['radius'], resolution=sphere_resolution )

        if vtk.VTK_MAJOR_VERSION <= 5:
            appendFilter.AddInputConnection(cylinder.GetProducerPort())
            appendFilter.AddInputConnection(sphere1.GetProducerPort())
            appendFilter.AddInputConnection(sphere2.GetProducerPort())
        else:
            sphere1Tri = vtk.vtkTriangleFilter()
            sphere1Tri.SetInputData(sphere1)
            sphere1Tri.Update()
            if length > 0:
                cylinderTri = vtk.vtkTriangleFilter()
                cylinderTri.SetInputData(cylinder)
                cylinderTri.Update()
                sphere2Tri = vtk.vtkTriangleFilter()
                sphere2Tri.SetInputData(sphere2)
                sphere2Tri.Update()

                booleanOperation1 = vtk.vtkBooleanOperationPolyDataFilter()
                booleanOperation1.SetOperationToUnion()
                booleanOperation2 = vtk.vtkBooleanOperationPolyDataFilter()
                booleanOperation2.SetOperationToUnion()
                # booleanOperation.SetInputData(0, cyl)
                booleanOperation1.SetInputData(0, cylinderTri.GetOutput())
                booleanOperation1.SetInputData(1, sphere1Tri.GetOutput())
                booleanOperation1.Update()
                booleanOperation2.SetInputData(0, booleanOperation1.GetOutput())
                booleanOperation2.SetInputData(1, sphere2Tri.GetOutput())
                # booleanOperation.SetInputData(2, sph2)
                booleanOperation2.Update()
            else:
                booleanOperation2 = sphere1Tri

            # this is simple version
            # appendFilter.AddInputData(booleanOperation2.GetOutput())
            # print "object connected, starting addind to general space"
            if appendedData is None:
                appendedData = booleanOperation2.GetOutput()
            else:
                booleanOperation3 = vtk.vtkBooleanOperationPolyDataFilter()
                booleanOperation3.SetOperationToUnion()
                booleanOperation3.SetInputData(0, appendedData)
                booleanOperation3.SetInputData(1, booleanOperation2.GetOutput())
                booleanOperation3.Update()
                appendedData = booleanOperation3.GetOutput()

            # import ipdb; ipdb.set_trace()
    # import ipdb; ipdb.set_trace()

    # appendFilter.Update()
    # appendedData = appendFilter.GetOutput()
    return appendedData

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

    import ipdb; ipdb.set_trace()
    writer = vtk.vtkPolyDataWriter()
    writer.SetFileName(outfile)
    try:
        writer.SetInputData(polyData)
    except:
        logger.warning("old vtk is used")
        writer.SetInput(polyData)
    writer.Write()


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
        '-l','--label',
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
