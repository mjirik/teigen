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
                 # endDistMultiplicator=1,
                 # use_joints=True,
                 voxelsize_mm=[1.0, 1.0, 1.0],
                 # voxelsize_mm_z=1.0,
                 # voxelsize_mm_x=1.0,
                 # voxelsize_mm_y=1.0,
                 areasize_px=[100, 100, 100],
                 # area_shape_z=100,
                 # area_shape_x=100,
                 # area_shape_y=100,
                 element_number=30,
                 uniform_radius_distribution=True,
                 normal_radius_distribution=False,
                 fixed_radius_distribution=False,
                 radius_distribution_minimum=2.0,
                 radius_distribution_maximum=10.0,
                 radius_distribution_mean=5.0,
                 radius_distribution_standard_deviation=5.0,
                 # intensity_profile=None
                 intensity_profile_radius=[0.4, 0.7, 1.0, 1.3],
                 intensity_profile_intensity=[195, 190, 200, 30],
                 random_generator_seed=0
                 ):
        """
        gtree is information about input data structure.
        endDistMultiplicator: make cylinder shorter by multiplication of radius
        intensity_profile: Dictionary type. Key is radius and value is required intensity.
        """
        # area_shape = [area_shape_z,area_shape_x, area_shape_y]
        # voxelsize_mm = [
        #     voxelsize_mm_z,
        #     voxelsize_mm_x,
        #     voxelsize_mm_y
        # ]
        self.build = build
        # self.filename = "output{:05d}.jpg"
        self.areasize_px = np.asarray(areasize_px)
        self.voxelsize_mm = np.asarray(voxelsize_mm)
        self.element_number = element_number
        self.radius_maximum = radius_distribution_maximum
        self.radius_minimum = radius_distribution_minimum
        # self.intensity_profile = intensity_profile
        self.intensity_profile = dict(zip(intensity_profile_radius, intensity_profile_intensity))
        self._cylinder_nodes = []
        self._cylinder_nodes_radiuses = []
        self.random_generator_seed = random_generator_seed
        self.radius_generator = self._const
        self.radius_generator_args=[radius_distribution_mean]
        self.area_volume = np.prod(self.areasize_px * self.voxelsize_mm)
        if uniform_radius_distribution:
            self.radius_generator = np.random.uniform
            self.radius_generator_args = [radius_distribution_minimum, radius_distribution_maximum]
        if normal_radius_distribution:
            self.radius_generator = np.random.normal
            self.radius_generator_args = [radius_distribution_mean, radius_distribution_standard_deviation]
        # import ipdb; ipdb.set_trace()
        # input of geometry and topology
        # self.V = []
        # self.CV = []
        # self.joints = {}
        # self.joints_lar = []
        # self.gtree = gtree
        # self.endDistMultiplicator = endDistMultiplicator
        # self.use_joints = use_joints
        self.surface = 0
        # self.LEN_STEP_CONSTANT = 0.1
        self.MAKE_IT_SHORTER_CONSTANT = 2.0
        # self.DIST_MAX_RADIUS_MULTIPLICATOR = 3.0
        self.OVERLAPS_ALOWED = False
        self.tree_data = {}
        self.data3d = None
        self.progress_callback = None

    def _const(self, value):
        return value


    def _cylinder_collision(self, pt1, pt2, radius):
        # TODO use geometry3.check_collision_along_line
        collision, new_nodes = g3.cylinder_collision(
            pt1,
            pt2,
            radius,
            other_points=self._cylinder_nodes,
            other_points_radiuses=self._cylinder_nodes_radiuses,
            areasize_mm=self.areasize_px,
            # DIST_MAX_RADIUS_MULTIPLICATOR=self.DIST_MAX_RADIUS_MULTIPLICATOR,
            OVERLAPS_ALOWED=self.OVERLAPS_ALOWED
        )

        if not collision:
            self._cylinder_nodes.extend(new_nodes)
            self._cylinder_nodes_radiuses.extend([radius] * len(new_nodes))

        return collision



        # if pt1 is not None \
        #     and self._is_in_area(pt1) \
        #     and self._is_in_area(pt2):
        #
        #     if self.OVERLAPS_ALOWED:
        #         return True
        #
        #     if len(self._cylinder_nodes) == 0:
        #         return True
        #     else:
        #         line_nodes = g3.get_points_in_line_segment(pt1, pt2, step)
        #         safe_dist2 = (self.radius_maximum * self.DIST_MAX_RADIUS_MULTIPLICATOR) ** 2
        #         for node in line_nodes:
        #             dist_closest = g3.closest_node_square_dist(node, self._cylinder_nodes)
        #             if dist_closest < safe_dist2:
        #                 return False
        #         return True
        # return False

    def run(self):


        tree_data = {

        }
        np.random.seed(self.random_generator_seed)
        pts = np.random.random([self.element_number, 3]) * self.areasize_px * self.voxelsize_mm

        # construct voronoi
        import scipy.spatial
        import itertools
        vor3 = scipy.spatial.Voronoi(pts)
        self.geometry_data = {
            "length":[],
            "radius":[],
            "volume":[],
            "surface":[],
            "vector":[]
        }

        # radius = self.radius_maximum
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
                    radius = self.radius_generator(*self.radius_generator_args)
                    if radius > self.radius_maximum:
                        continue
                    if radius < self.radius_minimum:
                        continue
                    pt1 = vor3.vertices[two_points_id[0]]
                    pt2 = vor3.vertices[two_points_id[1]]
                    pt1, pt2 = self._make_cylinder_shorter(pt1, pt2, radius*self.MAKE_IT_SHORTER_CONSTANT)
                    pt1 = np.asarray(pt1)
                    pt2 = np.asarray(pt2)
                    if not self._cylinder_collision(pt1, pt2, radius):


                        edge = {
                            "nodeA_ZYX_mm": pt1,
                            "nodeB_ZYX_mm": pt2,
                            # "radius_mm": radius
                            # "radius_mm": 1 + np.random.rand() * (self.max_radius -1 )
                            "radius_mm": radius,
                        }
                        tree_data[i] = edge
                        # line_nodes = g3.get_points_in_line_segment(pt1, pt2, radius)
                        # self._cylinder_nodes.extend(line_nodes)
                        length = np.linalg.norm(pt1 - pt2)
                        surf = 2 * np.pi * (radius + length)
                        volume =  np.pi * radius**2 * length
                        vector = pt1 - pt2

                        self.geometry_data["length"].append(length)
                        self.geometry_data["surface"].append(surf)
                        self.geometry_data["radius"].append(radius)
                        self.geometry_data["volume"].append(volume)
                        self.geometry_data["vector"].append(vector)
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

        self.tree_data = tree_data
        if self.build:
            from ..tree import TreeBuilder

            tvg = TreeBuilder('vtk')
            # yaml_path = os.path.join(path_to_script, "./hist_stats_test.yaml")
            # tvg.importFromYaml(yaml_path)
            tvg.voxelsize_mm = self.voxelsize_mm
            tvg.shape = self.areasize_px
            tvg.tree_data = tree_data
            output = tvg.buildTree() # noqa
            # tvg.show()
            tvg.saveToFile("tree_output.vtk")

        self.getStats()
        self.data3d = None


    def getStats(self):
        # self.assertTrue(False)
        print "Surface: ", self.surface
        import pandas as pd
        df = pd.DataFrame(self.geometry_data)
        desc = df.describe()

        print desc
        return df

    def generate_volume(self):
        from ..tree import TreeBuilder

        self.tvgvol = TreeBuilder('vol')
        self.tvgvol.voxelsize_mm = self.voxelsize_mm # [1, 1, 1]
        self.tvgvol.shape = self.areasize_px # [100, 100, 100]
        self.tvgvol.tree_data = self.tree_data
        self.tvgvol.finish_progress_callback = self.progress_callback
        if self.intensity_profile is not None:
            self.tvgvol.intensity_profile = self.intensity_profile
        self.data3d = self.tvgvol.buildTree()
        return self.data3d


    def saveVolumeToFile(self, filename="output{:06d}.jpg"):
        if self.data3d is None:
            self.generate_volume()

        # self.tvgvol.saveToFile(filename)
        import io3d
        import io3d.misc
        import numpy as np
        data = {
            'data3d': self.data3d.astype(np.uint8), #* self.output_intensity,
            'voxelsize_mm': self.voxelsize_mm,
            # 'segmentation': np.zeros_like(self.data3d, dtype=np.int8)
        }
        io3d.write(data, filename)

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
            radius = self.radius_maximum
        return g3.is_in_area(node, self.areasize_px, radius=radius)

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
    # if args.outputfile is not None:
    #     tg.saveToFile(args.outputfile, args.outputfiletype)


# class TreeGenerator(TreeConstructor):
#     """
#     back compatibility
#     """
#     pass

if __name__ == "__main__":
    main()

