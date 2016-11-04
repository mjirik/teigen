#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © %YEAR% %USER% <%MAIL%>
#
# Distributed under terms of the %LICENSE% license.

"""
%HERE%
"""

import logging

logger = logging.getLogger(__name__)
import argparse
import numpy as np
# import ..geometry3 as g3
# from ..geometry3d import plane_fit
from .. import geometry3d as g3
import os.path


def __half_plane(self, perp, plane_point, point):
    cdf = (np.array(point) - np.array(plane_point))
    out = perp[0] * cdf[0] + \
          perp[1] * cdf[1] + \
          perp[2] * cdf[2]
    return  out > 0

class CylinderGenerator:

    def __init__(self,
                 build=True,
                 gtree=None,
                 endDistMultiplicator=1,
                 use_joints=True,
                 element_number=60
                 ):
        """
        gtree is information about input data structure.
        endDistMultiplicator: make cylinder shorter by multiplication of radius
        """
        self.build = build
        self.area_shape = np.ones([3]) * 100
        self.max_radius = 3
        self.element_number = element_number
        self._cylinder_nodes = []
        # input of geometry and topology
        # self.V = []
        # self.CV = []
        # self.joints = {}
        # self.joints_lar = []
        # self.gtree = gtree
        # self.endDistMultiplicator = endDistMultiplicator
        # self.use_joints = use_joints
        self.surface = 0


    def _check_cylinder_position(self, pt1, pt2, step):

        if pt1 is not None \
            and self._is_in_area(pt1) \
            and self._is_in_area(pt2):

            if len(self._cylinder_nodes) == 0:
                return True
            else:
                line_nodes = g3.get_points_in_line_segment(pt1, pt2, step)
                safe_dist2 = (self.max_radius * 1.5)**2
                for node in line_nodes:
                    dist_closest = g3.closest_node_square_dist(node, self._cylinder_nodes)
                    if dist_closest < safe_dist2:
                        return False
                return True
        return False

    def run(self):


        tree_data = {

        }
        np.random.seed(0)
        pts = np.random.random([self.element_number, 3]) * 100

        # construct voronoi
        import scipy.spatial
        import itertools
        vor3 = scipy.spatial.Voronoi(pts)

        radius = self.max_radius

        # for i, two_points in enumerate(vor3.ridge_points):
        for i, simplex in enumerate(vor3.ridge_vertices):
            simplex = np.asarray(simplex)
            # fallowing line removes all ridges with oulayers
            simplex = simplex[simplex > 0]
            if np.all(simplex >= 0):

                x = vor3.vertices[simplex, 0]
                y = vor3.vertices[simplex, 1]
                z = vor3.vertices[simplex, 2]
                for two_points_id in itertools.combinations(simplex, 2):
                    pt1 = vor3.vertices[two_points_id[0]]
                    pt2 = vor3.vertices[two_points_id[1]]
                    pt1, pt2 = self._make_cylinder_shorter(pt1, pt2, radius*2)
                    pt1 = np.asarray(pt1)
                    pt2 = np.asarray(pt2)
                    if self._check_cylinder_position(pt1, pt2, radius):


                        edge = {
                            "nodeA_ZYX_mm": pt1,
                            "nodeB_ZYX_mm": pt2,
                            "radius_mm": radius
                        }
                        tree_data[i] = edge
                        line_nodes = g3.get_points_in_line_segment(pt1, pt2, radius)
                        self._cylinder_nodes.extend(line_nodes)
                        surf = 2 * np.pi * (radius + np.linalg.norm(pt1 - pt2))
                        self.surface += surf
            else:
                pass

        show_input_points = False
        if show_input_points:
            length = len(tree_data)
            for i in range(self.element_number):
                edge = {
                    #         #"nodeA_ZYX_mm": np.random.random(3) * 100,
                    "nodeA_ZYX_mm": pts[i-1],
                    "nodeB_ZYX_mm": pts[i],
                    #         "nodeB_ZYX_mm": np.random.random(3) * 100,
                    "radius_mm": 1
                }
                tree_data[i+length] = edge
        if self.build:
            from ..tree import TreeBuilder

            tvg = TreeBuilder('vtk')
            # yaml_path = os.path.join(path_to_script, "./hist_stats_test.yaml")
            # tvg.importFromYaml(yaml_path)
            tvg.voxelsize_mm = [1, 1, 1]
            tvg.shape = [100, 100, 100]
            tvg.tree_data = tree_data
            output = tvg.buildTree() # noqa
            # tvg.show()
            tvg.saveToFile("tree_output.vtk")


            tvgvol = TreeBuilder('vol')
            tvgvol.voxelsize_mm = [1, 1, 1]
            tvgvol.shape = [100, 100, 100]
            tvgvol.tree_data = tree_data
            outputvol = tvgvol.buildTree()
            tvgvol.saveToFile("tree_volume.pklz")
        # self.assertTrue(False)

        print "Surface: ", self.surface

    def _make_cylinder_shorter(self, nodeA, nodeB, radius): #, radius, cylinder_id):
        vector = (np.asarray(nodeA) - np.asarray(nodeB)).tolist()
        if np.linalg.norm(vector) < 2*radius:
            return None, None

        # mov circles to center of cylinder by size of radius because of joint
        nodeA = g3.translate(nodeA, vector,
                             -radius) # * self.endDistMultiplicator)
        nodeB = g3.translate(nodeB, vector,
                             radius) #  * self.endDistMultiplicator)
        return nodeA, nodeB

    def _is_in_area(self, node, radius=None):
        """
        check if point is in area with considering eventual maximum radius
        :param node:
        :param radius:
        :return:
        """
        node = np.asarray(node)
        if radius is None:
            radius = self.max_radius

        if np.all(node > (0 + radius)) and np.all(node < (self.area_shape - radius)):
            return  True
        else:
            return False



    def add_cylinder(self, nodeA, nodeB, radius, cylinder_id):

        try:
            idA = tuple(nodeA)  # self.gtree.tree_data[cylinder_id]['nodeIdA']
            idB = tuple(nodeB)  # self.gtree.tree_data[cylinder_id]['nodeIdB']
        except:
            idA = 0
            idB = 0
            self.use_joints = False

        # vect = nodeA - nodeB
        # self.__draw_circle(nodeB, vect, radius)

        vector = (np.asarray(nodeA) - np.asarray(nodeB)).tolist()

        # mov circles to center of cylinder by size of radius because of joint
        nodeA = g3.translate(nodeA, vector,
                             -radius * self.endDistMultiplicator)
        nodeB = g3.translate(nodeB, vector,
                             radius * self.endDistMultiplicator)

        if all(nodeA == nodeB):
            logger.error("End points are on same place")

        ptsA, ptsB = g3.cylinder_circles(nodeA, nodeB, radius,
                                         element_number=30)
        CVlistA = self.__construct_cylinder_end(ptsA, idA)
        CVlistB = self.__construct_cylinder_end(ptsB, idB)

        CVlist = CVlistA + CVlistB

        self.CV.append(CVlist)

# lar add ball
#         ball0 = mapper.larBall(radius, angle1=PI, angle2=2*PI)([10, 16])
#         V, CV = ball0
#         # mapper.T
#         # ball = STRUCT(MKPOLS(ball0))
#
#         # mapper.T(1)(nodeA[0])(mapper.T(2)(nodeA[1])(mapper.T(3)(nodeA[1])(ball)))
#
#         lenV = len(self.V)
#
#         self.V = self.V + (np.array(V) + np.array(nodeA)).tolist()
#         self.CV = self.CV + (np.array(CV) + lenV).tolist()


def main():
    logging.basicConfig()
    logger = logging.getLogger()
    logger.setLevel(logging.WARNING)

    # input parser
    parser = argparse.ArgumentParser(
        description='Histology analyser reporter. Try: \
python src/tb_volume.py -i ./tests/hist_stats_test.yaml'
    )
    parser.add_argument(
        '-i', '--inputfile',
        default=None,
        # required=True,
        help='input file, yaml file'
    )
    parser.add_argument(
        '-o', '--outputfile',
        default=None,
        help='output file, .raw, .dcm, .tiff, given by extension '
    )
    parser.add_argument(
        '-ot', '--outputfiletype',
        default='pkl',
        help='output file type.  raw, dcm, tiff, or pkl,   default is pkl, '
    )
    parser.add_argument(
        '-vs', '--voxelsize',
        default=[1.0, 1.0, 1.0],
        type=float,
        metavar='N',
        nargs='+',
        help='size of voxel (ZYX)'
    )
    parser.add_argument(
        '-ds', '--datashape',
        default=[200, 200, 200],
        type=int,
        metavar='N',
        nargs='+',
        help='size of output data in pixels for each axis (ZYX)'
    )
    parser.add_argument(
        '-g', '--generator',
        default='vol',
        type=str,
        help='Volume or surface model can be generated by use this option. \
                Use "vol", "volume" for volumetric model. For LAR surface model\
                use "lar". For VTK file use "vtk".'
    )
    parser.add_argument(
        '-d', '--debug', action='store_true',
        help='Debug mode')
    parser.add_argument(
        '-l', '--useLar', action='store_true',
        help='Use LAR')
    args = parser.parse_args()

    if args.debug:
        logger.setLevel(logging.DEBUG)

    # startTime = datetime.now()

    generator_params = None
    generator_class = args.generator

    # if args.generator == "vtk":
    #     import gen_vtk_tree
    #     gen_vtk_tree.vt2vtk_file(args.inputfile, args.outputfile)
    #     return
    cylgen = CylinderGenerator()
    cylgen.run()

    # tg = TreeBuilder(generator_class, generator_params)
    # tg.importFromYaml(args.inputfile)
    # tg.voxelsize_mm = args.voxelsize
    # tg.shape = args.datashape
    # tg.use_lar = args.useLar
    # data3d = tg.buildTree()
    #
    # logger.info("TimeUsed:" + str(datetime.now() - startTime))
    # # volume_px = sum(sum(sum(data3d)))
    # # volume_mm3 = volume_px * \
    # #     (tg.voxelsize_mm[0] * tg.voxelsize_mm[1] * tg.voxelsize_mm[2])
    # # logger.info("Volume px:" + str(volume_px))
    # # logger.info("Volume mm3:" + str(volume_mm3))
    #
    # # vizualizace
    # logger.debug("before visualization")
    # tg.show()
    # logger.debug("after visualization")

    # ukládání do souboru
    if args.outputfile is not None:
        tg.saveToFile(args.outputfile, args.outputfiletype)


# class TreeGenerator(TreeConstructor):
#     """
#     back compatibility
#     """
#     pass

if __name__ == "__main__":
    main()

